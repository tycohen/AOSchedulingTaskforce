
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

Best telescope (excluding Arecibo) for each pulsar based on RMS estimates

|  Pulsar |   Best Telescope(s)     |Total RMS (&mu;s)|   2nd Best Telescope(s)     |Total RMS (&mu;s)|
|---------|-------------------------|-----------------|-------------------------|-----------------|
|J0605+3757 | CHIME-GBTUWBR | 1.0799 | CHIME-GBTL | 1.1473|
|J1312+0051 | CHIME-GBTUWBR | 2.2679 | CHIME-GBTL | 2.4773|
|J2022+2534 | CHIME-GBTUWBR | 1.0320 | CHIME-GBTL | 1.0636|
|J1853+1303 | CHIME-GBTUWBR | 1.0487 | CHIME-GBTL | 1.1687|
|J0709+0458 | CHIME-GBTUWBR | 7.9440 | CHIME-GBTL | 8.6264|
|J2010-1323 | CHIME-GBTUWBR | 0.6917 | CHIME-GBTL | 0.7625|
|J1024-0719 | CHIME-GBTUWBR | 0.6070 | CHIME-GBTL | 0.6614|
|J1630+3550 | CHIME-GBTUWBR | 3.2956 | CHIME-GBTL | 3.5511|
|J1719-1438 | CHIME-GBTUWBR | 1.3918 | CHIME-GBTL | 1.4206|
|J2150-0326 | CHIME-GBTUWBR | 0.8654 | CHIME-GBTL | 0.9844|
|J0740+6620 | CHIME-GBTL | 0.5293 | CHIME-GBTUWBR | 0.5401|
|J0931-1902 | CHIME-GBTUWBR | 1.1245 | CHIME-GBTL | 1.1845|
|J1944+0907 | CHIME-GBTUWBR | 1.4603 | CHIME-GBTL | 1.6198|
|J1744-1134 | CHIME-GBTUWBR | 0.3236 | CHIME-GBTL | 0.3608|
|J0406+3039 | GBT_Rcvr_800-Rcvr_1_2 | 1.2742 | CHIME-GBTL | 1.4940|
|J1741+1351 | CHIME-GBTUWBR | 1.3960 | CHIME-GBTL | 1.5731|
|J0751+1807 | CHIME-GBTUWBR | 0.7779 | CHIME-GBTL | 0.8236|
|J0732+2314 | CHIME-GBTUWBR | 2.0237 | CHIME-GBTL | 2.0428|
|J2124-3358 | CHIME-GBTUWBR | 0.8718 | CHIME-GBTL | 0.8875|
|J1918-0642 | CHIME-GBTUWBR | 0.7934 | CHIME-GBTL | 0.8865|
|J1453+1902 | CHIME-GBTUWBR | 3.2892 | CHIME-GBTL | 3.5434|
|B1855+09 | CHIME-GBTUWBR | 0.4293 | CHIME-GBTL | 0.4507|
|J2039-3616 | CHIME-GBTUWBR | 1.0548 | CHIME-GBTL | 1.0986|
|J2043+1711 | CHIME-GBTUWBR | 1.1719 | CHIME-GBTL | 1.2911|
|J1713+0747 | CHIME-GBTUWBR | 0.1461 | CHIME-GBTL | 0.1591|
|J2145-0750 | CHIME-GBTUWBR | 0.5921 | CHIME-GBTL | 0.6329|
|J2322+2057 | CHIME-GBTUWBR | 2.3165 | CHIME-GBTL | 2.5141|
|J0636+5128 | CHIME-GBTUWBR | 0.4552 | CHIME-GBTL | 0.4771|
|J0509+0856 | CHIME-GBTUWBR | 1.8154 | CHIME-GBTL | 1.9033|
|J0154+1833 | CHIME-GBTUWBR | 1.1320 | CHIME-GBTL | 1.2629|
|J1903+0327 | CHIME-GBTUWBR | 126.5030 | GBT_Rcvr_800-Rcvr_1_2 | 257.6498|
|J0645+5158 | CHIME-GBTUWBR | 1.4924 | CHIME-GBTL | 1.4945|
|J2017+0603 | CHIME-GBTUWBR | 1.5020 | CHIME-GBTL | 2.0490|
|J1923+2515 | CHIME-GBTUWBR | 1.8922 | CHIME-GBTL | 2.0482|
|J0557+1551 | CHIME-GBTUWBR | 6.2106 | CHIME-GBTL | 9.4408|
|J1705-1903 | GBT_Rcvr_800-Rcvr_1_2 | 0.6364 | CHIME-GBTL | 1.4167|
|J0030+0451 | CHIME-GBTUWBR | 0.7301 | CHIME-GBTL | 0.7848|
|J2317+1439 | CHIME-GBTUWBR | 0.6241 | CHIME-GBTL | 0.7017|
|J1327+3423 | CHIME-GBTUWBR | 7.0927 | CHIME-GBTL | 7.4655|
|J2234+0944 | CHIME-GBTUWBR | 1.0049 | CHIME-GBTL | 1.1113|
|J1745+1017 | CHIME-GBTUWBR | 0.7148 | CHIME-GBTL | 0.7589|
|J1640+2224 | CHIME-GBTUWBR | 0.6984 | CHIME-GBTL | 0.7659|
|J1911+1347 | CHIME-GBTUWBR | 0.9063 | CHIME-GBTL | 0.9820|
|J0610-2100 | CHIME-GBTL | 1.3959 | CHIME-GBTUWBR | 1.4633|
|J1832-0836 | CHIME-GBTUWBR | 0.7109 | CHIME-GBTL | 0.7748|
|J1803+1358 | GBT_Rcvr_800-Rcvr_1_2 | 1.6802 | CHIME-GBTUWBR | 3.7014|
|J0340+4130 | GBT_Rcvr_800-Rcvr_1_2 | 2.0188 | CHIME-GBTUWBR | 2.2533|
|J1730-2304 | CHIME-GBTL | 0.6779 | CHIME-GBTUWBR | 0.6806|
|J0023+0923 | CHIME-GBTUWBR | 1.2225 | CHIME-GBTL | 1.3534|
|J1012-4235 | CHIME-GBTUWBR | 1.1273 | GBT_Rcvr_800-Rcvr_1_2 | 1.6197|
|J0613-0200 | GBT_Rcvr_800-Rcvr_1_2 | 0.3512 | CHIME-GBTL | 0.6102|
|J1747-4036 | GBT_Rcvr_800-Rcvr_1_2 | 8.1213 | CHIME-GBTL | 34.6205|
|J1012+5307 | CHIME-GBTL | 0.2603 | CHIME-GBTUWBR | 0.2635|
|J1125+7819 | CHIME-GBTL | 0.4302 | CHIME-GBTUWBR | 0.4622|
|B1953+29 | GBT_Rcvr_800-Rcvr_1_2 | 3.4353 | CHIME-GBTUWBR | 8.4933|
|J1811-2405 | GBT_Rcvr_800-Rcvr_1_2 | 1.5610 | CHIME-GBTL | 3.9791|
|J1751-2857 | CHIME-GBTUWBR | 1.8187 | CHIME-GBTL | 1.8704|
|J2033+1734 | CHIME-GBTUWBR | 2.5483 | CHIME-GBTL | 2.7857|
|J1843-1113 | GBT_Rcvr_800-Rcvr_1_2 | 1.0656 | CHIME-GBTL | 2.6128|
|J2234+0611 | CHIME-GBTUWBR | 1.2458 | CHIME-GBTL | 1.3622|
|J2302+4442 | CHIME-GBTUWBR | 1.2123 | CHIME-GBTL | 1.2493|
|J1909-3744 | CHIME-GBTUWBR | 0.2481 | CHIME-GBTL | 0.2789|
|J1802-2124 | GBT_Rcvr_800-Rcvr_1_2 | 38.9509 | CHIME-GBTL | 127.7888|
|J1630+3734 | CHIME-GBTUWBR | 0.9388 | CHIME-GBTL | 1.0046|
|B1937+21 | GBT_Rcvr_800-Rcvr_1_2 | 0.3972 | CHIME-GBTL | 2.0839|
|J0125-2327 | CHIME-GBTUWBR | 0.2452 | CHIME-GBTL | 0.2826|
|J1910+1256 | CHIME-GBTUWBR | 1.3944 | CHIME-GBTL | 1.5184|
|J1946+3417 | GBT_Rcvr_800-Rcvr_1_2 | 1.5064 | CHIME-GBTL | 4.8968|
|J1643-1224 | GBT_Rcvr_800-Rcvr_1_2 | 0.5250 | CHIME-GBTL | 2.0532|
|J1455-3330 | CHIME-GBTUWBR | 1.5372 | CHIME-GBTL | 1.7465|
|J0614-3329 | CHIME-GBTUWBR | 0.5088 | CHIME-GBTL | 0.5277|
|J2229+2643 | CHIME-GBTUWBR | 1.1285 | CHIME-GBTL | 1.2324|
|J2214+3000 | CHIME-GBTUWBR | 1.8393 | CHIME-GBTL | 2.4495|
|J1614-2230 | CHIME-GBTL | 0.8259 | CHIME-GBTUWBR | 0.8546|
|J1022+1001 | CHIME-GBTL | 0.5063 | CHIME-GBTUWBR | 0.5121|
|J1738+0333 | CHIME-GBTUWBR | 1.4134 | CHIME-GBTL | 1.7284|
|J1600-3053 | GBT_Rcvr_800-Rcvr_1_2 | 0.3965 | CHIME-GBTUWBR | 0.5824|


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

