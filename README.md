
# AO Scheduling Taskforce &sigma;TOA Estimates


Code used to aggregate pulsar parameters and estimate the TOA uncertainty referenced to the infinite frequency TOA (as described in [Lam et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...861...12L/abstract)) for NANOGrav 15-yr pulsars with a variety of telescopes including:
* Arecibo 430 & L-band
* Arecibo L-band and S-band
* GBT 800 & L-band
* GBT L-band & VLA S-band
* CHIME & GBT L-band
* CHIME & GBT UWBR

NG15yr_totRMS.txt contains the total RMS for each pulsar at each telescope. See Example Usage to get individual noise components.

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

