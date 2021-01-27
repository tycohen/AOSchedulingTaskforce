
# AO Scheduling Taskforce &sigma;TOA Estimates


Code used to aggregate pulsar parameters and estimate the TOA uncertainty referenced to the infinite frequency TOA (as described in [Lam et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...861...12L/abstract)) for NANOGrav 15-yr pulsars with a variety of telescopes including:
* Arecibo 430 & L-band
* Arecibo L-band and S-band
* GBT 800 & L-band
* GBT L-band & VLA S-band
* CHIME & GBT L-band
* CHIME & GBT UWBR


NG15yr_totRMS.txt contains the total RMS for each pulsar at each telescope, assuming 30 min integration time per band per epoch, except at CHIME where the integration time is dec-dependent and at most 30 hours / epoch. Total noise can be scaled to different integration time using the individual noise components. See **Example Usage** to get individual noise components.

## tl;dr

|    Pulsar       |   Best Telescope(s)     |    Total RMS (&mu;s)   |
|-----------------|-------------------------|------------------------|
|J0605+3757  |  CHIME-GBTL  |  1.2063  |
|J1312+0051  |  CHIME-GBTL  |  2.6885  |
|J2022+2534  |  CHIME-GBTL  |  1.1379  |
|J1853+1303  |  CHIME-GBTL  |  1.2463  |
|J0709+0458  |  CHIME-GBTL  |  9.4069  |
|J2010-1323  |  CHIME-GBTL  |  0.8306  |
|J1024-0719  |  CHIME-GBTL  |  0.7260  |
|J1630+3550  |  CHIME-GBTL  |  3.6182  |
|J1719-1438  |  CHIME-GBTL  |  1.5405  |
|J2150-0326  |  CHIME-GBTL  |  1.0460  |
|J0740+6620  |  CHIME-GBTL  |  0.5465  |
|J0931-1902  |  CHIME-GBTL  |  1.2916  |
|J1944+0907  |  CHIME-GBTL  |  1.7100  |
|J1744-1134  |  CHIME-GBTL  |  0.3929  |
|J0406+3039  |  GBT_Rcvr_800-Rcvr_1_2  |  1.2704  |
|J1741+1351  |  CHIME-GBTL  |  1.6506  |
|J0751+1807  |  CHIME-GBTL  |  0.8939  |
|J0732+2314  |  CHIME-GBTL  |  2.1619  |
|J2124-3358  |  CHIME-GBTL  |  0.9422  |
|J1918-0642  |  CHIME-GBTL  |  0.9583  |
|J1453+1902  |  CHIME-GBTL  |  3.8104  |
|B1855+09  |  CHIME-GBTL  |  0.4843  |
|J2039-3616  |  CHIME-GBTL  |  1.2106  |
|J2043+1711  |  CHIME-GBTL  |  1.3668  |
|J1713+0747  |  CHIME-GBTL  |  0.1726  |
|J2145-0750  |  CHIME-GBTL  |  0.6594  |
|J2322+2057  |  CHIME-GBTL  |  2.6748  |
|J0636+5128  |  CHIME-GBTL  |  0.4985  |
|J0509+0856  |  CHIME-GBTL  |  2.0392  |
|J0154+1833  |  CHIME-GBTL  |  1.2971  |
|J1903+0327  |  GBT_Rcvr_1_2_VLAS  |  70.5327  |
|J0645+5158  |  CHIME-GBTL  |  1.5587  |
|J2017+0603  |  CHIME-GBTL  |  2.2493  |
|J1923+2515  |  CHIME-GBTL  |  2.1760  |
|J0557+1551  |  GBT_Rcvr_1_2_VLAS  |  9.5636  |
|J1705-1903  |  GBT_Rcvr_800-Rcvr_1_2  |  0.6350  |
|J0030+0451  |  CHIME-GBTL  |  0.8435  |
|J2317+1439  |  CHIME-GBTL  |  0.7227  |
|J1327+3423  |  CHIME-GBTL  |  7.7983  |
|J2234+0944  |  CHIME-GBTL  |  1.2209  |
|J1745+1017  |  CHIME-GBTL  |  0.8236  |
|J1640+2224  |  CHIME-GBTL  |  0.8007  |
|J1911+1347  |  CHIME-GBTL  |  1.0752  |
|J0610-2100  |  CHIME-GBTL  |  1.5220  |
|J1832-0836  |  CHIME-GBTL  |  0.8467  |
|J1803+1358  |  GBT_Rcvr_800-Rcvr_1_2  |  1.6841  |
|J0340+4130  |  GBT_Rcvr_800-Rcvr_1_2  |  2.0146  |
|J1730-2304  |  CHIME-GBTL  |  0.7036  |
|J0023+0923  |  CHIME-GBTL  |  1.4426  |
|J1012-4235  |  GBT_Rcvr_800-Rcvr_1_2  |  1.6141  |
|J0613-0200  |  GBT_Rcvr_800-Rcvr_1_2  |  0.3506  |
|J1747-4036  |  GBT_Rcvr_1_2_VLAS  |  3.1925  |
|J1012+5307  |  CHIME-GBTL  |  0.2674  |
|J1125+7819  |  CHIME-GBTL  |  0.4325  |
|B1953+29  |  GBT_Rcvr_800-Rcvr_1_2  |  3.4281  |
|J1811-2405  |  GBT_Rcvr_800-Rcvr_1_2  |  1.5581  |
|J1751-2857  |  CHIME-GBTL  |  2.0446  |
|J2033+1734  |  CHIME-GBTL  |  2.9056  |
|J1843-1113  |  GBT_Rcvr_800-Rcvr_1_2  |  1.0631  |
|J2234+0611  |  CHIME-GBTL  |  1.4945  |
|J2302+4442  |  CHIME-GBTL  |  1.3122  |
|J1909-3744  |  CHIME-GBTL  |  0.3082  |
|J1802-2124  |  GBT_Rcvr_1_2_VLAS  |  10.7155  |
|J1630+3734  |  CHIME-GBTL  |  1.0463  |
|B1937+21  |  GBT_Rcvr_1_2_VLAS  |  0.1572  |
|J0125-2327  |  CHIME-GBTL  |  0.3062  |
|J1910+1256  |  CHIME-GBTL  |  1.6490  |
|J1946+3417  |  GBT_Rcvr_800-Rcvr_1_2  |  1.5023  |
|J1643-1224  |  GBT_Rcvr_800-Rcvr_1_2  |  0.5238  |
|J1455-3330  |  CHIME-GBTL  |  1.9242  |
|J0614-3329  |  CHIME-GBTL  |  0.5798  |
|J2229+2643  |  CHIME-GBTL  |  1.2725  |
|J2214+3000  |  CHIME-GBTL  |  2.6973  |
|J1614-2230  |  CHIME-GBTL  |  0.9013  |
|J1022+1001  |  CHIME-GBTL  |  0.5227  |
|J1738+0333  |  CHIME-GBTL  |  1.9000  |
|J1600-3053  |  GBT_Rcvr_800-Rcvr_1_2  |  0.3959  |

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

`sigma_best(self, exclude="*")`

Get the best instrument for each pulsar
and return list of tuples of (pulsar name, instrument, sigma_tot)
Set exclude = string to ignore a particular telescope

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
print(pta.best_sigma(exclude="AO"))
```

