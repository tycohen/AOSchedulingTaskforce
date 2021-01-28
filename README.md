
# AO Scheduling Taskforce &sigma;TOA Estimates


Code used to aggregate pulsar parameters and estimate the TOA uncertainty referenced to the infinite frequency TOA (as described in [Lam et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...861...12L/abstract)) for NANOGrav 15-yr pulsars with a variety of telescopes including:
* Arecibo 430 & L-band
* Arecibo L-band and S-band
* GBT 800 & L-band
* GBT L-band & VLA S-band
* CHIME & GBT L-band
* CHIME & GBT UWBR


NG15yr_totRMS.txt contains the total RMS for each pulsar at each telescope, assuming 30 min integration time per band per epoch, except at CHIME where the integration time is dec-dependent and at most 30 hours / epoch. Total noise can be scaled to different integration time using the individual noise components. See **Example Usage** to get individual noise components.

## Summary

Best telescope (excluding Arecibo and GBT UWBR) for each pulsar based on RMS estimates taking into consideration declination-limits.

|  Pulsar |   Best Telescope(s)     |Total RMS (&mu;s)|   2nd Best Telescope(s)     |Total RMS (&mu;s)|
|---------|-------------------------|-----------------|-------------------------|-----------------|
|J0605+3757 | CHIME-GBTL | 1.1552 | GBT_Rcvr_800-Rcvr_1_2 | 2.0358|
|J1312+0051 | CHIME-GBTL | 2.5590 | GBT_Rcvr_800-Rcvr_1_2 | 3.7960|
|J2022+2534 | CHIME-GBTL | 1.0419 | GBT_Rcvr_800-Rcvr_1_2 | 1.1488|
|J1853+1303 | CHIME-GBTL | 1.2095 | GBT_Rcvr_800-Rcvr_1_2 | 1.9298|
|J0709+0458 | CHIME-GBTL | 8.7690 | GBT_Rcvr_1_2_VLAS | 12.1082|
|J2010-1323 | CHIME-GBTL | 0.7913 | GBT_Rcvr_800-Rcvr_1_2 | 1.1157|
|J1024-0719 | CHIME-GBTL | 0.6790 | GBT_Rcvr_800-Rcvr_1_2 | 0.9731|
|J1630+3550 | CHIME-GBTL | 3.6222 | GBT_Rcvr_800-Rcvr_1_2 | 9.3529|
|J1719-1438 | CHIME-GBTL | 1.4084 | GBT_Rcvr_800-Rcvr_1_2 | 1.8039|
|J2150-0326 | CHIME-GBTL | 1.0480 | GBT_Rcvr_800-Rcvr_1_2 | 1.6443|
|J0740+6620 | CHIME-GBTL | 0.5070 | GBT_Rcvr_800-Rcvr_1_2 | 1.1347|
|J0931-1902 | CHIME-GBTL | 1.1971 | GBT_Rcvr_800-Rcvr_1_2 | 1.5985|
|J1944+0907 | CHIME-GBTL | 1.6993 | GBT_Rcvr_800-Rcvr_1_2 | 2.7582|
|J1744-1134 | CHIME-GBTL | 0.3753 | GBT_Rcvr_800-Rcvr_1_2 | 0.5293|
|J0406+3039 | GBT_Rcvr_800-Rcvr_1_2 | 1.2742 | CHIME-GBTL | 1.4583|
|J1741+1351 | CHIME-GBTL | 1.6546 | GBT_Rcvr_800-Rcvr_1_2 | 2.9330|
|J0751+1807 | CHIME-GBTL | 0.8289 | GBT_Rcvr_800-Rcvr_1_2 | 1.1757|
|J0732+2314 | CHIME-GBTL | 2.0457 | GBT_Rcvr_800-Rcvr_1_2 | 2.7689|
|J2124-3358 | GBT_Rcvr_800-Rcvr_1_2 | 1.0778 | GBT_Rcvr_1_2_VLAS | 1.2549|
|J1918-0642 | CHIME-GBTL | 0.9180 | GBT_Rcvr_800-Rcvr_1_2 | 1.3219|
|J1453+1902 | CHIME-GBTL | 3.6239 | GBT_Rcvr_800-Rcvr_1_2 | 5.6557|
|B1855+09 | CHIME-GBTL | 0.4579 | GBT_Rcvr_800-Rcvr_1_2 | 0.6280|
|J2039-3616 | GBT_Rcvr_800-Rcvr_1_2 | 1.4707 | GBT_Rcvr_1_2_VLAS | 1.7814|
|J2043+1711 | CHIME-GBTL | 1.3399 | GBT_Rcvr_800-Rcvr_1_2 | 2.2708|
|J1713+0747 | CHIME-GBTL | 0.1602 | GBT_Rcvr_1_2_VLAS | 0.1969|
|J2145-0750 | CHIME-GBTL | 0.6875 | GBT_Rcvr_800-Rcvr_1_2 | 1.0324|
|J2322+2057 | CHIME-GBTL | 2.5864 | GBT_Rcvr_800-Rcvr_1_2 | 4.3082|
|J0636+5128 | CHIME-GBTL | 0.3880 | GBT_Rcvr_800-Rcvr_1_2 | 0.9103|
|J0509+0856 | CHIME-GBTL | 1.9387 | GBT_Rcvr_800-Rcvr_1_2 | 2.7891|
|J0154+1833 | CHIME-GBTL | 1.3340 | GBT_Rcvr_800-Rcvr_1_2 | 2.8961|
|J1903+0327 | GBT_Rcvr_1_2_VLAS | 48.3718 | CHIME-GBTL | 245.6239|
|J0645+5158 | CHIME-GBTL | 1.4743 | GBT_Rcvr_800-Rcvr_1_2 | 2.9557|
|J2017+0603 | GBT_Rcvr_1_2_VLAS | 1.6264 | CHIME-GBTL | 2.0987|
|J1923+2515 | CHIME-GBTL | 2.0951 | GBT_Rcvr_800-Rcvr_1_2 | 3.4809|
|J0557+1551 | GBT_Rcvr_1_2_VLAS | 5.3571 | CHIME-GBTL | 9.4253|
|J1705-1903 | GBT_Rcvr_1_2_VLAS | 0.5512 | GBT_Rcvr_800-Rcvr_1_2 | 0.6364|
|J0030+0451 | CHIME-GBTL | 0.8110 | GBT_Rcvr_800-Rcvr_1_2 | 1.2065|
|J2317+1439 | CHIME-GBTL | 0.7490 | GBT_Rcvr_800-Rcvr_1_2 | 1.5689|
|J1327+3423 | CHIME-GBTL | 7.5790 | GBT_Rcvr_800-Rcvr_1_2 | 13.8936|
|J2234+0944 | CHIME-GBTL | 1.1274 | GBT_Rcvr_1_2_VLAS | 1.5500|
|J1745+1017 | CHIME-GBTL | 0.7703 | GBT_Rcvr_800-Rcvr_1_2 | 1.0975|
|J1640+2224 | CHIME-GBTL | 0.7894 | GBT_Rcvr_800-Rcvr_1_2 | 1.4290|
|J1911+1347 | CHIME-GBTL | 0.9963 | GBT_Rcvr_800-Rcvr_1_2 | 1.4606|
|J0610-2100 | GBT_Rcvr_800-Rcvr_1_2 | 1.6001 | GBT_Rcvr_1_2_VLAS | 2.2534|
|J1832-0836 | CHIME-GBTL | 0.7827 | GBT_Rcvr_800-Rcvr_1_2 | 1.0559|
|J1803+1358 | GBT_Rcvr_800-Rcvr_1_2 | 1.6802 | CHIME-GBTL | 4.0567|
|J0340+4130 | GBT_Rcvr_800-Rcvr_1_2 | 2.0188 | CHIME-GBTL | 2.6424|
|J1730-2304 | GBT_Rcvr_800-Rcvr_1_2 | 0.8273 | GBT_Rcvr_1_2_VLAS | 1.1629|
|J0023+0923 | CHIME-GBTL | 1.4126 | GBT_Rcvr_800-Rcvr_1_2 | 2.2837|
|J1012-4235 | GBT_Rcvr_1_2_VLAS | 1.1461 | GBT_Rcvr_800-Rcvr_1_2 | 1.6197|
|J0613-0200 | GBT_Rcvr_800-Rcvr_1_2 | 0.3512 | CHIME-GBTL | 0.5378|
|J1747-4036 | GBT_Rcvr_1_2_VLAS | 2.7376 | GBT_Rcvr_800-Rcvr_1_2 | 8.1213|
|J1012+5307 | CHIME-GBTL | 0.2555 | GBT_Rcvr_800-Rcvr_1_2 | 0.5858|
|J1125+7819 | CHIME-GBTL | 0.3854 | GBT_Rcvr_800-Rcvr_1_2 | 1.5188|
|B1953+29 | GBT_Rcvr_800-Rcvr_1_2 | 3.4353 | GBT_Rcvr_1_2_VLAS | 5.6623|
|J1811-2405 | GBT_Rcvr_800-Rcvr_1_2 | 1.5610 | GBT_Rcvr_1_2_VLAS | 1.6729|
|J1751-2857 | GBT_Rcvr_800-Rcvr_1_2 | 2.3538 | GBT_Rcvr_1_2_VLAS | 3.3697|
|J2033+1734 | CHIME-GBTL | 2.9154 | GBT_Rcvr_800-Rcvr_1_2 | 5.2802|
|J1843-1113 | GBT_Rcvr_800-Rcvr_1_2 | 1.0656 | GBT_Rcvr_1_2_VLAS | 1.3051|
|J2234+0611 | CHIME-GBTL | 1.3866 | GBT_Rcvr_800-Rcvr_1_2 | 2.0389|
|J2302+4442 | CHIME-GBTL | 1.0490 | GBT_Rcvr_800-Rcvr_1_2 | 2.1463|
|J1909-3744 | GBT_Rcvr_800-Rcvr_1_2 | 0.3791 | GBT_Rcvr_1_2_VLAS | 0.6362|
|J1802-2124 | GBT_Rcvr_1_2_VLAS | 7.9382 | GBT_Rcvr_800-Rcvr_1_2 | 38.9509|
|J1630+3734 | CHIME-GBTL | 1.0158 | GBT_Rcvr_800-Rcvr_1_2 | 1.9963|
|B1937+21 | GBT_Rcvr_1_2_VLAS | 0.1365 | GBT_Rcvr_800-Rcvr_1_2 | 0.3972|
|J0125-2327 | GBT_Rcvr_1_2_VLAS | 0.2866 | GBT_Rcvr_800-Rcvr_1_2 | 0.3791|
|J1910+1256 | CHIME-GBTL | 1.5494 | GBT_Rcvr_800-Rcvr_1_2 | 2.3068|
|J1946+3417 | GBT_Rcvr_800-Rcvr_1_2 | 1.5064 | GBT_Rcvr_1_2_VLAS | 1.7415|
|J1643-1224 | GBT_Rcvr_800-Rcvr_1_2 | 0.5250 | GBT_Rcvr_1_2_VLAS | 0.7314|
|J1455-3330 | GBT_Rcvr_800-Rcvr_1_2 | 2.4215 | GBT_Rcvr_1_2_VLAS | 3.8961|
|J0614-3329 | GBT_Rcvr_800-Rcvr_1_2 | 0.6971 | GBT_Rcvr_1_2_VLAS | 0.9803|
|J2229+2643 | CHIME-GBTL | 1.2753 | GBT_Rcvr_800-Rcvr_1_2 | 2.6486|
|J2214+3000 | GBT_Rcvr_1_2_VLAS | 2.0539 | CHIME-GBTL | 2.4623|
|J1614-2230 | GBT_Rcvr_800-Rcvr_1_2 | 0.9900 | GBT_Rcvr_1_2_VLAS | 1.4292|
|J1022+1001 | CHIME-GBTL | 0.5194 | GBT_Rcvr_800-Rcvr_1_2 | 0.6706|
|J1738+0333 | GBT_Rcvr_1_2_VLAS | 1.7060 | CHIME-GBTL | 1.7663|
|J1600-3053 | GBT_Rcvr_1_2_VLAS | 0.2576 | GBT_Rcvr_800-Rcvr_1_2 | 0.3965|


## 15yr_psrs.txt

Contains the parameters of each 15yr pulsar. J0437 not included b/c only visible with VLA.

### Columns

* name
* period (s)
* DM (pc cm^-3)
* DEC: declination (deg)
* dtd: scintillation timescale from NE2001 (s)
* dnud: scintillation bandwidth from NE2001 (GHz)
* taud: scattering timescale from NE2001 (us)
* dist: DM distance from NE2001 (kpc)
* w50: FWHM of the L-band template (us)
* weff: effective width of the L-band template (us)
* uscale: scales intensity across pulse phase
* S_1000: flux density at 1 GHz (mJy)
* spindex: spectral index computed from log(flux ratio)/log(freq ratio)
using 800/1400 for GBT pulsars, 430/1400 for AO when available, otherwise
1400/2000. Cutoff at zero
* sig_jitter: single pulse RMS jitter (us) from [Lam et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...872..193L/abstract) or 0.63 * W50 if not in 12.5-yr (from Fig. 7 of same paper) 


## The Pulsar class


Class to store an individual 15-yr pulsar

### Attributes: 

* `name` : string

pulsar name
* `period` : float

pulse period in seconds
* `dm` : float

dispersion measure in pc cm^-3	   
* `dec` : float

declination in degrees
* `dtd` : float

scintillation timescale in seconds
* `dnud` : float

scintillation bandwidth in GHz
* `taud` : float

scattering timescale in us
* `dist` : float

DM distance in kpc	    
* `w50` : float

FWHM of L-band pulse profile in us
* `weff` : float

Effective width of the L-band profile in us
* `uscale` : float

scaling factor to distribute intensity across
                pulse profile
* `s_1000` : float

flux density at 1 GHz in mJy
* `spindex` : float

spectral index
* `sig_j_single` : float

single-pulse RMS jitter in us
* `sigmas` : dict

dictionary of dictionaries of RMS components
                for each instrument


### Methods

`sigma_jitter(self, t_int)`

Return intrinsic jitter noise (in us)
for given integration time in seconds

`add_sigmas(self, instr_name, sigma_tup)`

instr_name : str

name of timing instrument

sigma_tup : tuple

tuple of RMS components

`get_instr_keys(self)`

returns instrument key names for sigmas dict


## The PTA class


Class to store timed pulsars

### Attributes: 

* `name` : string (optional)

PTA name
* `psrlist` : list

list of pulsar.Pulsar objects


### Methods

`get_single_pulsar(self, psr_name)`

return pulsar.Pulsar object whose name matches 'psr_name'

`sigma_best(self, exclude=[])`

Get the best instrument for each pulsar
and return list of tuples of (pulsar name, instrument, sigma_tot)

<ins>Parameters</ins>

`exclude` : list (optional, default = [])

list of telescope name substrings to exclude when sorting RMS's




`write_to_text(self, filename)`

Write total RMS for each pulsar at each instrument to file


## Example Usage

To get a particular pulsar's noise estimates at a particular telescope:

```
import cPickle
with open('NG15yr.pta', 'rb') as f:
    pta = cPickle.load(f)
j1713 = pta.get_single_pulsar("J1713+0747")
print(j1713.get_instr_keys()) # get keys for sigmas dict
print(j1713.sigmas["AO_430_Lwide_logain"])
```

To see which non-Arecibo telescope is best for each pulsar

```
import cPickle
with open('NG15yr.pta', 'rb') as f:
    pta = cPickle.load(f)
print(pta.best_sigma(exclude=["AO"]))
```

