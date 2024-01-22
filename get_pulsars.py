import os
import sys
import re
import cPickle
import warnings
import pandas as pd
import numpy as np
import pypulse
from glob import glob
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from PTAOptimizer.foppulsar import FOpPulsar
from resultsreader import ResultsReader
import ne2001_15y
from pulsar import Pulsar
from pta import PTA
"""
Get 15-year pulsar inputs to frequencyoptimizer.PulsarNoise
from various sources. Distances are all DM distances even if
pulsar has measured parallax.
"""

BACKENDS = ['puppi', 'guppi', 'yuppi', 'VEGAS'] # in order of precedence 
LBAND_RCVRS = {'puppi': 'L-wide',
               'guppi': 'Rcvr1_2',
               'yuppi': 'L-Band',
               'VEGAS': 'Rcvr1_2'} 
DATADIR = '/lustre/aoc/users/pdemores/timing/toagen/data'
BACKUPDATADIR = '/lustre/aoc/users/pdemores/timing/15yr_prelim/working'
LOCALDATADIR = '.'
TEMPLATEDIR = os.path.join(DATADIR, 'templates')
LOCALTEMPLATEDIR = os.path.join(LOCALDATADIR, 'templates')

def read_txt_to_df():
    df = pd.read_csv('15yr_psrs.txt',
                     header=0,
                     delimiter='\s+')
    df.columns.values[0] = df.columns.values[0].strip('#')
    return df

def write_df_to_txt(df):
    df.rename(columns={'name': '#name'}, inplace=True)
    df.to_csv('15yr_psrs.txt',
              sep='\t',
              index=False)
    df.rename(columns={'#name': 'name'}, inplace=True)
    return

def get_period():
    df = read_txt_to_df()
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

def get_DM():
    df = read_txt_to_df()
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

def get_dec():
    '''
    Search for parfiles for list of pulsar names, then search parfiles
    for sky coordinates and convert to declination
    '''
    df = read_txt_to_df()
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
                                              frame='barycentricmeanecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
                elif re.match('ELONG\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('ELAT\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentricmeanecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
                elif re.match('BETA\s+', line) is not None:
                    lat = float(line.split()[1])
                    for line in f:
                        if re.match('LAMBDA\s+', line) is not None:
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentricmeanecliptic')
                            dec = ecoord.icrs.dec.degree
                            lon = float(line.split()[1])
                            break
                elif re.match('LAMBDA\s+', line) is not None:
                    lon = float(line.split()[1])
                    for line in f:
                        if re.match('BETA\s+', line) is not None:
                            lat = float(line.split()[1])
                            ecoord = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, 
                                              frame='barycentricmeanecliptic')
                            dec = ecoord.icrs.dec.degree
                            break
        decs.append(dec)
    return decs

def get_ne2001_pars(df):
    """
    Read NE2001 params from pyne2001 dict
    Parameters
    ----------
    df : pandas.DataFrame

    Returns
    -------
    dtds, dnuds, tauds, dists
    """
    ne2001_dict = ne2001_15y.ne2001_results
    dtds = []
    dnuds = []
    tauds = []
    dists = []
    for n in df['name']:
        try:
            dtds.append(ne2001_dict[n]['SCINTIME'])
            dnuds.append(ne2001_dict[n]['SBW'] / 1000.) # MHz to GHz
            tauds.append(ne2001_dict[n]['TAU'] * 1000.) # ms to us
            dists.append(ne2001_dict[n]['DIST'])
        except KeyError:
            dtds.append(np.nan)
            dnuds.append(np.nan)
            tauds.append(np.nan)
            dists.append(np.nan)
    return dtds, dnuds, tauds, dists

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
            ar = pypulse.archive.Archive(tf)
            temp_units = ar.getTimeUnit()
            if temp_units == 'SEC':
                sp = pypulse.singlepulse.SinglePulse(ar.getData(),
                                                     period=p0)
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

def make_pta(df, write=False):
    pta = PTA(name='NG15yr')
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
    if write:
        with open('NG15yr.pta', 'wb') as f:
            cPickle.dump(pta, f)
    return pta


def get_parfile(name):
    for b in BACKENDS:
        try:
            for f in os.listdir(os.path.join(DATADIR, b, name)):
                if '.par' in f:
                    return os.path.join(DATADIR, b, name, f)
        except OSError:
            pass
    try:
        for f in os.listdir(os.path.join(BACKUPDATADIR, name)):
            if '.par' in f:
                return os.path.join(BACKUPDATADIR, name, f)
    except OSError:
        pass
    local_parfiles = glob(os.path.join(LOCALDATADIR,
                                       "parfiles",
                                       "*{}*.par".format(name)))
    try:
        return local_parfiles[0]
    except IndexError:
            pass
    warnings.warn("No parfile for {}".format(name))
    return None

def get_template(name):
    for b in BACKENDS:
        r = LBAND_RCVRS[b]
        temppath = glob(os.path.join(TEMPLATEDIR,
                                     "{}.{}.{}.*.x.sum.sm".format(name,
                                                                  r,
                                                                  b.upper())))
        temppath += glob(os.path.join(LOCALTEMPLATEDIR,
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
