import cPickle
import numpy as np
from os import path
from astropy.coordinates import SkyCoord
import astropy.units as u
import frequencyoptimizer as fop
from PTAOptimizer.telescope import Telescope
import PTAOptimizer.observatory_ops as oops


def calc_timing(pta,
                nus,
                rxspecfile=None,
                t_int=1800.,
                dec_lim=None,
                lat=None,
                gainmodel=None,
                gainexp=None,
                timefac=0.):
    if rxspecfile is None:
        raise ValueError('rxspecfile must be defined')
    for p in pta.psrlist:
        scope = Telescope(name=path.basename(rxspecfile).strip(".txt"),
                          dec_lim=dec_lim,
                          lat=lat,
                          gainmodel=gainmodel,
                          gainexp=gainexp)
        scope.timefac = timefac
        ra_str = p.name[1:3] + 'h' + p.name[3:5] + 'm' # get RA from Jname
        j2k_coords = SkyCoord(ra=ra_str, dec=p.dec*u.deg, frame='icrs')
        # initial scope noise to get the rx specs
        scope_noise_init = fop.TelescopeNoise(1.,
                                              1.,
                                              T=t_int,
                                              rxspecfile=rxspecfile)
        scope_noise_init.gain = oops.get_gains(scope,
                                               p.dec,
                                               scope_noise_init.get_gain(nus))
        if 0. in scope_noise_init.gain:
            # if any gains are zero, psr below at least 1 scopes horizon
            p.add_sigmas(scope.name, (-2, -2, -2, -2, -2))
            continue
        else:
            if isinstance(timefac, np.ndarray):
                scope_noise_init.T = get_tobs(scope_noise_init.get_T(nus),
                                         scope,
                                         p.dec)
            else:
                scope_noise_init.T = scope_noise_init.get_T(nus)
            # if 0. in scope_noise_init.T:
            #     print('Zero in scope_noise_init.T : {}'.format(scope_noise_init.T))
            #     print('Gain = {}'.format(scope_noise_init.gain))

            # re-initialize telescope noise with dec-dependent gains, int time
            scope_noise = fop.TelescopeNoise(rx_nu=nus,
                                             gain=scope_noise_init.gain,
                                             T_rx=scope_noise_init.get_T_rx(nus),
                                             epsilon=scope_noise_init.get_epsilon(nus),
                                             T=scope_noise_init.T)
            p.telescope_noise.update({scope.name : scope_noise})
            pulsar_noise = fop.PulsarNoise('', 
                                           alpha=-1 * p.spindex,
                                           dtd=p.dtd,
                                           dnud=p.dnud,
                                           taud=p.taud,
                                           C1=1.16,
                                           I_0=p.s_1000,
                                           DM=p.dm,
                                           D=p.dist,
                                           tauvar=0.5 * p.taud,
                                           Weffs=p.weff,
                                           W50s=p.w50,
                                           Uscale=p.uscale,
                                           sigma_Js=p.sigma_jitter(scope_noise.T),
                                           glon=j2k_coords.galactic.b.degree,
                                           glat=j2k_coords.galactic.l.degree)
            gal_noise = fop.GalacticNoise()
            fop_inst = fop.FrequencyOptimizer(pulsar_noise,
                                              gal_noise,
                                              scope_noise,
                                              nchan=len(nus),
                                              numax=get_ctrfreq(nus),
                                              numin=get_ctrfreq(nus),
                                              vverbose=False)
            sigma_tup = fop_inst.calc_single(nus)
            p.add_sigmas(scope.name, sigma_tup)
    return

def get_tobs(t0, scope, psr_dec, horiz=0., cutoff=1.08e5):
    if abs(psr_dec - scope.lat) >= 90. - horiz:
        # source never rises
        if isinstance(t0, (list, np.ndarray)): 
            t_obs = np.zeros(len(t0))
        elif isinstance(t0, (int, float)):
            t_obs = 0.
    elif abs(psr_dec + scope.lat) >= 90. + horiz:
        # source never sets
        t_obs = t0 * (2 * np.cos(np.radians(psr_dec)) ** -1) ** scope.timefac
    else:
        t_obs = t0 * (np.cos(np.radians(psr_dec)) ** -1) ** scope.timefac
    return np.clip(t_obs, 0., cutoff)

def get_ctrfreq(nus):
    mid = float(len(nus)) / 2
    if len(nus) % 2 != 0:
        return nus[int(mid - .5)]
    else:
        return (nus[int(mid)] + nus[int(mid - 1)]) / 2.

def chime_only(write=False):
    """Returns a CHIME-only pta.PTA object"""
    with open('NG15yr.pta', 'rb') as ptaf:
        pta = cPickle.load(ptaf)
    for p in pta.psrlist:   # reset the sigma dicts to being empty
        p.sigmas = {}
    print('Timing CHIME')
    chime_nus = np.arange(0.6 - 0.4 / 2., 0.6 + 0.4 / 2., 0.005)[:-1]
    chime_gainexp = np.full(len(chime_nus), 1.)
    chime_timefac = np.full(len(chime_nus), 1.)
    calc_timing(pta,
                chime_nus,
                rxspecfile="CHIME.txt",
                dec_lim=(90., -35.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_gainexp,
                timefac=chime_timefac)
    if write:
        with open('NG15yr_CHIMEonly.pta', 'wb') as ptaf:
            cPickle.dump(pta, ptaf)
    return pta
    
if __name__ == '__main__':
    """
    'Main' function; calculate sigmas for multiple telescope configs
    and writes out to .pta file. File is overwritten each time.
    """
    with open('NG15yr.pta', 'rb') as ptaf:
        pta = cPickle.load(ptaf)
    print('Timing AO L-S')
    LbandSlo_nus = np.arange(1.44 - .618 / 2, 1.868, 0.011)
    Shi_nus = np.arange(2.227 - .354 / 2, 2.227 + .354 / 2, 0.01)[:-1]
    aoLS_nus = np.sort(np.concatenate([LbandSlo_nus, Shi_nus]))
    calc_timing(pta,
                aoLS_nus,
                rxspecfile="AO_Lwide_Swide_logain.txt",
                dec_lim=(39., 0.),
                t_int=1800.,
                lat=18.44,
                gainmodel=None,
                gainexp=None)

    print('Timing AO 430-L')
    nus_ao430 = np.arange(.432 - .02 / 2, .432 + .02 / 2, 0.00125)[:-1]
    nus_aoL = np.arange(1.44 - .58 / 2, 1.44 + .58 / 2, 0.00125)[:-1]
    ao430L_nus = np.concatenate([nus_ao430, nus_aoL])
    calc_timing(pta,
                ao430L_nus,
                rxspecfile="AO_430_Lwide_logain.txt",
                t_int=1800.,
                dec_lim=(39., 0.),
                lat=18.44,
                gainmodel=None,
                gainexp=None)
    
    print('Timing GB 800-1200')
    nus_gb800 = np.arange(.820 - .200 / 2, .820 + .200 / 2, 0.009)
    nus_gb1_2 = np.arange(1.510 - .800 / 2, 1.510 + .800 / 2, 0.009)[:-1]
    gbt80012_nus = np.concatenate([nus_gb800, nus_gb1_2])
    calc_timing(pta,
                gbt80012_nus,
                rxspecfile="GBT_Rcvr_800-Rcvr_1_2_logain.txt",
                t_int=1800.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GBT-L + VLA-S')
    vlaS_nus = np.arange(3. - 2. / 2., 3. + 2. / 2., 0.009)
    gbL_vlaS_nus = np.concatenate([nus_gb1_2, vlaS_nus])
    calc_timing(pta,
                gbL_vlaS_nus,
                rxspecfile="GBT_Rcvr_1_2_VLAS_logain.txt",
                t_int=1800.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing CHIME + GBT-L')
    chime_nus = np.arange(0.6 - 0.4 / 2., 0.6 + 0.4 / 2., 0.009)[:-1]
    chime_gbtL_nus = np.concatenate([chime_nus, nus_gb1_2])
    chime_gbtL_gainexp = np.concatenate([np.full(len(chime_nus), 1.),
                                         np.full(len(nus_gb1_2), 0.)])
    chime_gbtL_timefac = np.concatenate([np.full(len(chime_nus), 1.),
                                         np.full(len(nus_gb1_2), 0.)])
    calc_timing(pta,
                chime_gbtL_nus,
                rxspecfile="CHIME-GBTL_logain.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_gbtL_gainexp,
                timefac=chime_gbtL_timefac)
    
    print('Timing CHIME + UWBR')
    gbuwb_ctrfreq = 2.35 #GHz
    gbuwb_bw = 3.3 # GHz
    gbuwb_nus = np.arange(gbuwb_ctrfreq - gbuwb_bw / 2,
                        gbuwb_ctrfreq + gbuwb_bw / 2,
                        0.009)[:-1]
    chime_lt_uwbr_nus = chime_nus[chime_nus < gbuwb_nus[0]]
    chime_uwbr_nus = np.concatenate([chime_lt_uwbr_nus,
                                     gbuwb_nus])
    chime_uwbr_gainexp = np.concatenate([np.full(len(chime_lt_uwbr_nus), 1.),
                                         np.full(len(gbuwb_nus), 0.)])
    chime_uwbr_timefac = np.concatenate([np.full(len(chime_lt_uwbr_nus), 1.),
                                         np.full(len(gbuwb_nus), 0.)])
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    print('Timing CHIME + UWBR improved (no excess noise)')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_NoExcNoise.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    print('Timing CHIME + UWBR pessimistic (with excess noise)')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_WithExcNoise.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    # print('Timing CHIME + UWBR May 2023 rcvr specs (no RFI)')
    # calc_timing(pta,
    #             chime_uwbr_nus,
    #             rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_NoRFI_May23.txt",
    #             dec_lim=(90., -20.),
    #             lat=49.32,
    #             gainmodel='cos',
    #             gainexp=chime_uwbr_gainexp,
    #             timefac=chime_uwbr_timefac)

    # print('Timing CHIME + UWBR May 2023 rcvr specs (no RFI, improved gain)')
    # calc_timing(pta,
    #             chime_uwbr_nus,
    #             rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_NoRFI_May23_ImprovedGain.txt",
    #             dec_lim=(90., -20.),
    #             lat=49.32,
    #             gainmodel='cos',
    #             gainexp=chime_uwbr_gainexp,
    #             timefac=chime_uwbr_timefac)
    
    # print('Timing CHIME + UWBR May 2023 rcvr specs (no RFI, improved gain, improved Tsys)')
    # calc_timing(pta,
    #             chime_uwbr_nus,
    #             rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_NoRFI_May23_ImprovedGainTsys.txt",
    #             dec_lim=(90., -20.),
    #             lat=49.32,
    #             gainmodel='cos',
    #             gainexp=chime_uwbr_gainexp,
    #             timefac=chime_uwbr_timefac)
    
    print('Timing UWBR')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing UWBR improved (no excess noise)')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_NoExcNoise.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing UWBR pessimistic (with excess noise)')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_WithExcNoise.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    # print('Timing UWBR using May 2023 rcvr specs (no RFI)')
    # calc_timing(pta,
    #             gbuwb_nus,
    #             rxspecfile="./uwbr_rxspecs/GBTUWBR_NoRFI_May23.txt",
    #             t_int=3600.,
    #             dec_lim=(90., -46.),
    #             lat=38.42,
    #             gainmodel=None,
    #             gainexp=None)

    # print('Timing UWBR using May 2023 rcvr specs (no RFI, improved gain)')
    # calc_timing(pta,
    #             gbuwb_nus,
    #             rxspecfile="./uwbr_rxspecs/GBTUWBR_NoRFI_May23_ImprovedGain.txt",
    #             t_int=3600.,
    #             dec_lim=(90., -46.),
    #             lat=38.42,
    #             gainmodel=None,
    #             gainexp=None)
    
    # print('Timing UWBR using May 2023 rcvr specs (no RFI, improved gain, improved Tsys)')
    # calc_timing(pta,
    #             gbuwb_nus,
    #             rxspecfile="./uwbr_rxspecs/GBTUWBR_NoRFI_May23_ImprovedGainTsys.txt",
    #             t_int=3600.,
    #             dec_lim=(90., -46.),
    #             lat=38.42,
    #             gainmodel=None,
    #             gainexp=None)

    print('Timing GBTUWBR_May_2023_Improved_Gain_No_RFI')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_May_2023_Improved_Gain_No_RFI.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GBTUWBR_May_2023_OnSky_No_RFI')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_May_2023_OnSky_No_RFI.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GBTUWBR_May_2023_Improved_Gain_and_Tsys_No_RFI')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_May_2023_Improved_Gain_and_Tsys_No_RFI.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GBTUWBR_Original_Specs')
    calc_timing(pta,
                gbuwb_nus,
                rxspecfile="./uwbr_rxspecs/GBTUWBR_Original_Specs.txt",
                t_int=3600.,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing CHIME-GBTUWBR_May_2023_Improved_Gain_No_RFI')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_May_2023_Improved_Gain_No_RFI.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    print('Timing CHIME-GBTUWBR_May_2023_OnSky_No_RFI')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_May_2023_OnSky_No_RFI.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    print('Timing CHIME-GBTUWBR_May_2023_Improved_Gain_and_Tsys_No_RFI')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_May_2023_Improved_Gain_and_Tsys_No_RFI.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)

    print('Timing CHIME-GBTUWBR_Original_Specs')
    calc_timing(pta,
                chime_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/CHIME-GBTUWBR_Original_Specs.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_uwbr_gainexp,
                timefac=chime_uwbr_timefac)
    
    
    print('Timing GB140ft 400-800 MHz')
    gb140lo_nus = np.arange(.4, .8, .009)
    calc_timing(pta,
                gb140lo_nus,
                rxspecfile="GB140_400-800.txt",
                t_int=1.08e5,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GB140ft 1.3-1.8 GHz')
    gb140L_nus = np.arange(1.3, 1.8, .009)
    calc_timing(pta,
                gb140L_nus,
                rxspecfile="GB140_LBand.txt",
                t_int=1.08e5,
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing CHIME + GB140ft L-band')
    chime_gb140L_nus = np.concatenate([chime_nus,
                                       gb140L_nus])
    chime_gb140L_gainexp = np.concatenate([np.full(len(chime_nus), 1.),
                                         np.full(len(gb140L_nus), 0.)])
    chime_gb140L_timefac = np.concatenate([np.full(len(chime_nus), 1.),
                                         np.full(len(gb140L_nus), 0.)])
    calc_timing(pta,
                chime_gb140L_nus,
                rxspecfile="CHIME-GB140Lband_logain.txt",
                dec_lim=(90., -20.),
                lat=49.32,
                gainmodel='cos',
                gainexp=chime_gb140L_gainexp,
                timefac=chime_gb140L_timefac)

    print('Timing GB140ft 400-800 MHz + GBT-L')
    gb140lo_gbtL_nus = np.concatenate([gb140lo_nus, nus_gb1_2])
    calc_timing(pta,
                gb140lo_gbtL_nus,
                rxspecfile="GB140_400-800_GBTLband.txt",
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    print('Timing GB140ft 400-800 MHz + UWBR')
    gb140lo_uwbr_nus = np.concatenate([gb140lo_nus[gb140lo_nus < gbuwb_nus[0]],
                                       gbuwb_nus])
    calc_timing(pta,
                gb140lo_uwbr_nus,
                rxspecfile="./uwbr_rxspecs/GB140_400-800_GBTUWBR.txt",
                dec_lim=(90., -46.),
                lat=38.42,
                gainmodel=None,
                gainexp=None)

    with open('NG15yr.pta', 'wb') as ptaf:
        cPickle.dump(pta, ptaf)
