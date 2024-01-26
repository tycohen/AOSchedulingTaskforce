import os
from os import path
import sys
import re
import cPickle
import warnings
import argparse
import pandas as pd
import numpy as np
import pypulse
from glob import glob
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import pyne2001
from PTAOptimizer.foppulsar import FOpPulsar
from resultsreader import ResultsReader
from pulsar import Pulsar
from pta import PTA

SHORT_DESCRIPTION = """
Collect pulsar attributes to frequencyoptimizer.PulsarNoise
from various sources. Distances are all DM distances even if
pulsar has measured parallax.
"""

BACKENDS = ['puppi', 'guppi', 'yuppi', 'VEGAS'] # in order of precedence 
LBAND_RCVRS = {'puppi': 'L-wide',
               'guppi': 'Rcvr1_2',
               'yuppi': '1.5GHz',
               'VEGAS': 'Rcvr1_2'}
DATADIR = '/lustre/aoc/users/pdemores/timing/toagen/data'
TIMINGDIR = '/lustre/aoc/users/pdemores/timing'
BACKUPDATADIR = '/lustre/aoc/users/pdemores/timing/15yr_prelim/working'
LOCALDATADIR = '.'
TEMPLATEDIR = path.join(DATADIR, 'templates')
LOCALTEMPLATEDIR = path.join(LOCALDATADIR, 'templates')
VALID_COLNAMES = ['name',
                  'period',
                  'DM',
                  'DEC',
                  'RA',
                  'dtd',
                  'dnud',
                  'taud',
                  'dist',
                  'w50',
                  'weff',
                  'uscale',
                  'S_1000',
                  'spindex',
                  'sig_jitter']

def _write_bool(filepath, overwrite=None):
    """Return True if file doesnt exist, otherwise ask me"""
    while path.isfile(filepath):
        while overwrite not in ('y', 'n', 'yes', 'no'):
            overwrite = raw_input("File {} exists. "
                                  "Overwrite (y/n)?\n".format(filepath)).lower()
            if overwrite in ("yes", 'y'):
                return True
            elif overwrite in ("no", 'n'):
                return False
            else:
                break
    return True

def read_txt_to_df(inpath):
    df = pd.read_csv(inpath,
                     header=0,
                     delimiter='\s+')
    df.columns.values[0] = df.columns.values[0].strip('#')
    return df

def write_df_to_txt(df, outpath):
    if _write_bool(outpath):
        df.rename(columns={'name': '#name'}, inplace=True)
        df.to_csv(outpath,
                  sep='\t',
                  index=False)
        df.rename(columns={'#name': 'name'}, inplace=True)
    return

def get_period(df):
    periods = []
    for name in df['name']:
        parfile = get_parfile(name)
        if parfile is None:
            p0 = np.nan
        with open(parfile, 'r') as f:
            for line in f:
                matches = re.match('F0\s+', line)
                if matches is not None:
                    p0 = 1/float(line.split()[1])
                    break
                p0 = np.nan
        periods.append(p0)
    return periods

def get_DM(df):
    dms = []
    for name in df['name']:
        parfile = get_parfile(name)
        if parfile is None:
            dms.append(np.nan)
            continue
        with open(parfile, 'r') as f:
            for line in f:
                matches = re.match('DM\s+', line)
                if matches is not None:
                    dm = float(line.split()[1])
                    break
                dm = np.nan
        dms.append(dm)
    return dms

def get_ra(df):
    '''
    Search for parfiles for list of pulsar names, then search parfiles
    for sky coordinates and convert to right ascension
    '''
    ras = []
    for name in df['name']:
        parfile = get_parfile(name)
        if parfile is None:
            ra = np.nan
        with open(parfile, 'r') as f:
            for line in f:
                if re.match('RAJ\s+', line) is not None:
                    valstr = line.split()[1].split(':')
                    rastr = "".join(i for j in zip(valstr,
                                                    ['h', 'm', 's']) for i in j)
                    ra = Angle(rastr).degree
                    break
                elif re.match('ELAT\s+', line) is not None:
                    lat = float(line.split()[1])
                    for line in f:
                        if re.match('ELONG\s+', line) is not None:
                            lon = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            ra = ecoord.icrs.ra.degree
                            break
                elif re.match('ELONG\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('ELAT\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            ra = ecoord.icrs.ra.degree
                            break
                elif re.match('BETA\s+', line) is not None:
                    lat = float(line.split()[1])
                    for line in f:
                        if re.match('LAMBDA\s+', line) is not None:
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            ra = ecoord.icrs.ra.degree
                            lon = float(line.split()[1])
                            break
                elif re.match('LAMBDA\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('BETA\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            ra = ecoord.icrs.ra.degree
                            break
        ras.append(ra)
    return ras

def get_dec(df):
    '''
    Search for parfiles for list of pulsar names, then search parfiles
    for sky coordinates and convert to declination
    '''
    decs = []
    for name in df['name']:
        parfile = get_parfile(name)
        if parfile is None:
            dec = np.nan
        with open(parfile, 'r') as f:
            for line in f:
                if re.match('DECJ\s+', line) is not None:
                    valstr = line.split()[1].split(':')
                    decstr = "".join(i for j in zip(valstr,
                                                    ['d', 'm', 's']) for i in j)
                    dec = Angle(decstr).value
                    break
                elif re.match('ELAT\s+', line) is not None:
                    lat = float(line.split()[1])
                    for line in f:
                        if re.match('ELONG\s+', line) is not None:
                            lon = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
                elif re.match('ELONG\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('ELAT\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
                elif re.match('BETA\s+', line) is not None:
                    lat = float(line.split()[1])
                    for line in f:
                        if re.match('LAMBDA\s+', line) is not None:
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            dec = ecoord.icrs.dec.degree
                            lon = float(line.split()[1])
                            break
                elif re.match('LAMBDA\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('BETA\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentrictrueecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
        decs.append(dec)
    return decs

def get_ne2001_pars(ras, decs, dms):
    """
    Read NE2001 params from pyne2001 dict
    Parameters
    ----------
    ras : list or np.ndarray of right ascensions in deg
    decs : list or np.ndarray of declinations in deg

    Returns
    -------
    dtds, dnuds, tauds, dists
    """
    dtds = []
    dnuds = []
    tauds = []
    dists = []
    eq_coords = SkyCoord(ra=ras*u.degree,
                         dec=decs*u.degree,
                         frame="icrs")
    longs = eq_coords.galactic.l.degree
    lats = eq_coords.galactic.b.degree
    for l, b, dm in zip(longs, lats, dms):
        ne2001_dict = pyne2001.get_dist_full(l, b, dm)
        dtds.append(ne2001_dict['SCINTIME'])
        dnuds.append(ne2001_dict['SBW'] / 1000.) # MHz to GHz
        tauds.append(ne2001_dict['TAU'] * 1000.) # ms to us
        dists.append(ne2001_dict['DIST'])
    return dtds, dnuds, tauds, dists

def update_flux_spindex_file(newfile,
                             oldfile="psr_spect_index_flux_stat.txt"):
    """
    Replace old flux entries for a given pulsar-rcvr pair with a new entry
    and add any entries unique to either file
    """
    updated_lines = []
    with open(newfile, "r") as newf:
        newlines = newf.readlines()
    with open(oldfile, "r") as oldf:
        for line in oldf:
            oldstart = re.search("^([^\s]+)\s+([^\s]+)\s+([^\s]+)",
                                 line).group(0)
            for newl in newlines:
                if newl.startswith(oldstart):
                    updated_lines.append(newl)
                    newlines.remove(newl)
                    break
            else:
                updated_lines.append(line)
    updated_lines += newlines
    with open(oldfile, "w") as oldf:
        oldf.write("".join([l.replace("\r", "") for l in updated_lines]))
    return

def get_flux_spindex(df):
    """
    Computes fluxes and spectral indices from David's file
    Spectral indices computed as in 12.5 yr alpha = log(S2 / S1)/log(nu2 / nu1)
    where 800/1400 used for GBT pulsars, 430/1400 for AO when available and
    1400/2000 for AO when 430 unavailable
    Returns flux at 1GHz (mJy) and spectral index
    """
    with open('psr_spect_index_flux_stat.txt', 'r') as f:
        flux_1GHz = []
        spindices = []
        lines = f.readlines()
        names = [l.split()[0] for l in lines]
    for n in df['name']:
        if n == "J1910+1256":
            flux = 0.98
            spindex = -1.66
        elif n == "J1911-1114":
            flux = 1.3
            spindex = -2.9
        elif n == "J0437-4715":
            flux = 210.
            spindex = -1.0
        elif n not in names:
            flux = np.nan
            spindex = np.nan
            print('No fluxes for {}'.format(n))
        else:
            rcvrs = [l.split()[1] for l in lines if l.startswith(n)]
            fluxes = [float(l.split()[7]) for l in lines if l.startswith(n)]
            fluxdict = dict(zip(rcvrs, fluxes))
            if "430_PUPPI" in fluxdict:
                spindex = np.log10(fluxdict['L-wide_PUPPI'] / fluxdict['430_PUPPI']) / np.log10(1400. / 430.)
                flux = fluxdict['L-wide_PUPPI'] * (1. / 1.4) ** spindex
            elif "S-wide_PUPPI" in fluxdict:
                spindex = np.log10(fluxdict['S-wide_PUPPI'] / fluxdict['L-wide_PUPPI']) / np.log10(2000. / 1400.)
                flux = fluxdict['L-wide_PUPPI'] * (1. / 1.4) ** spindex
            else:
                spindex = np.log10(fluxdict['Rcvr1_2_GUPPI'] / fluxdict['Rcvr_800_GUPPI']) / np.log10(1400. / 800.)
                flux = fluxdict['Rcvr1_2_GUPPI'] * (1. / 1.4) ** spindex
        flux_1GHz.append(flux)
        if spindex > 0.:
            spindices.append(0.)
        else:
            spindices.append(spindex)
    return flux_1GHz, spindices
        

def get_template_pars(df):
    """
    Returns
    -------
    weffs, uscales, w50s
    """
    weffs = []
    uscales = []
    w50s = []
    for name, p0 in zip(df['name'], df['period']):
        tf = get_template(name)
        if tf is None:
            weffs.append(np.nan)
            uscales.append(np.nan)
            w50s.append(np.nan)
        else:
            ar = pypulse.archive.Archive(tf, verbose=False)
            temp_units = ar.getTimeUnit()
            if temp_units == 'SEC':
                sp = pypulse.singlepulse.SinglePulse(ar.getData(),
                                                     period=p0)
                sp.normalize()
            else:
                raise ValueError("Template has units of {}, "
                                 " need time units.".format(temp_units))
            w50s.append(sp.getFWHM() * 1e6) # s to us
            weffs.append(sp.getWeff() * 1e6) # s to us
            uscales.append(sp.getUscale())
    return weffs, uscales, w50s

def get_jitter(df):
    """Get single pulse RMS jitter from Michael's 12.5 yr paper"""
    rr = ResultsReader()
    sigma_js = []
    for n, w50 in zip(df['name'], df['w50']):
        try:
            jitmodel = rr.get_pulsar_model(n, "constant")
            if jitmodel[1][0]:
                s_j = jitmodel[0][0][0]
            else:
                s_j = jitmodel[0][0][3]
        except KeyError:
            # not in 12.5 yr jitter paper, use scaling from Fig.7 Lam et al. 2018
            # Freq-dependence of jitter
            s_j = 0.63 * w50
        sigma_js.append(s_j)
    return sigma_js

def make_pta(df, ptaname):
    pta = PTA(name=ptaname)
    for idx, psr in df.iterrows():
        pulsar = Pulsar(name=psr['name'],
                        period=psr['period'],
                        dm=psr['DM'],
                        dec=psr['DEC'],
                        dtd=psr['dtd'],
                        dnud=psr['dnud'],
                        taud=psr['taud'],
                        dist=psr['dist'],
                        w50=psr['w50'],
                        weff=psr['weff'],
                        uscale=psr['uscale'],
                        s_1000=psr['S_1000'],
                        spindex=psr['spindex'],
                        sig_j_single=psr['sig_jitter'])
        pta.psrlist.append(pulsar)
    ptapath = '{}.pta'.format(ptaname)
    if _write_bool(ptapath):
        with open(ptapath, 'wb') as f:
            cPickle.dump(pta, f)
    return pta

def get_parfile(name):
    for b in BACKENDS:
        try:
            for f in os.listdir(path.join(DATADIR, b, name)):
                if '.par' in f:
                    return path.join(DATADIR, b, name, f)
        except OSError:
            pass
    try:
        for f in os.listdir(path.join(BACKUPDATADIR, name)):
            if '.par' in f:
                return path.join(BACKUPDATADIR, name, f)
    except OSError:
        pass
    local_parfiles = glob(path.join(LOCALDATADIR,
                                       "parfiles",
                                       "*{}*.par".format(name)))
    try:
        return local_parfiles[0]
    except IndexError:
        pass
    warnings.warn("No parfile for {}".format(name))
    return None

def get_template(name):
    timing_yr = ["2013", "2015", "2016", "2017", "2019"]
    for b in BACKENDS:
        r = LBAND_RCVRS[b]
        # temppaths in order of search precedence
        temppath = glob(path.join(TEMPLATEDIR,
                                  "{}.{}.{}.*.x.sum.sm".format(name,
                                                               r,
                                                               b.upper())))
        temppath += glob(path.join(LOCALTEMPLATEDIR,
                                   "{}.{}.{}.*.x.sum.sm".format(name,
                                                                r,
                                                                b.upper())))
        # this catches some VLA L-band templates
        temppath += glob(path.join(DATADIR, b, name,
                                   "{}.{}.{}.*.P.sum.sm".format(name,
                                                                r,
                                                                b.upper())))
        for y in timing_yr:
            temppath += glob(path.join(TIMINGDIR,
                                       "nanograv_timing_{}".format(y),
                                       "templates",
                                       "{}.{}.{}.*.x.sum.sm".format(name,
                                                                    r,
                                                                    b.upper())))
        if len(temppath) > 0:
            return temppath[0]
    print("No template for {}".format(name))
    return None

def update_column(df, psrname, col_name, col_val):
    """
    Update (in-place) column value for a pulsar 

    Parameters
    ----------
    df : pandas.DataFrame
    psrname : str
              pulsar name
    col_name : str
               name of column to update with value
    col_val : float
              value to update column with
    """
    df.loc[df["name"] == psrname, [col_name]] = col_val
    return

def drop_pulsar(name, df):
    return df[~df['name'].isin([name])]

def parse_args(cmd_line_args):
    """
    Parse command line args
    """
    parser = argparse.ArgumentParser(description=SHORT_DESCRIPTION,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--infile", type=str, dest="infile",
                        default="15yr_psrs.txt",
                        help="ASCII text file to read into dataframe.\n" +\
                        "Must have the (order-agnostic) header:\n" +\
                        "#{}\n".format(" ".join(VALID_COLNAMES)) +\
                        "and at least the 'name' column filled in.")
    parser.add_argument("-o", "--outfile", type=str, dest="outfile",
                        default="15yr_psrs.txt",
                        help="ASCII text file to write dataframe to.\n" +\
                        "Default is to overwrite the input file.")
    parser.add_argument("-pta", "--ptaname", type=str, dest="ptaname",
                        default=None,
                        help="If supplied, creates a pickled pta.PTA object" +\
                        "in a file named <ptaname>.pta.")
    parser.add_argument("--debug", dest="debug",
                        action="store_true", default=False,
                        help="Toggle on the debug message")
    args = parser.parse_args(cmd_line_args)
    if args.debug:
        print(args)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = parse_args(sys.argv[1:])

    if not path.isfile(args.infile):
        raise IOError("'infile' not found: {}".format(args.infile))
    
    df = read_txt_to_df(args.infile)
    df["period"] = get_period(df)
    df["DM"] = get_DM(df)
    df["DEC"] = get_dec(df)
    df["RA"] = get_ra(df)
    dtds, dnuds, tauds, dists = get_ne2001_pars(list(df["RA"]),
                                                list(df["DEC"]),
                                                list(df["DM"]))
    df["dtd"] = dtds
    df["dnud"] = dnuds
    df["taud"] = tauds
    df["dist"] = dists
    weffs, uscales, w50s = get_template_pars(df)
    df["weff"] = weffs
    df["uscale"] = uscales
    df["w50"] = w50s
    flux_1GHz, spindices = get_flux_spindex(df)
    df["S_1000"] = flux_1GHz
    df["spindex"] = spindices
    df["sig_jitter"] = get_jitter(df)
    write_df_to_txt(df, args.outfile)
    with open("testdf.dataframe", "wb") as f:
        cPickle.dump(df, f)
    if args.ptaname:
        make_pta(df, args.ptaname)
    
if __name__ == '__main__':
    sys.exit(main())
