
# AO Scheduling Taskforce &sigma;TOA Estimates


Code used to aggregate pulsar parameters and estimate the TOA uncertainty referenced to the infinite frequency TOA (as described in [Lam et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...861...12L/abstract)) for NANOGrav 15-yr pulsars with a variety of telescopes including:
* Arecibo 430 & L-band
* Arecibo L-band and S-band
* GBT 800 & L-band
* GBT L-band & VLA S-band
* CHIME & GBT L-band
* CHIME & GBT UWBR



## The Pulsar class


Class to store an individual 15-yr pulsar

### Attributes: 

* name : string
            pulsar name
* period : float
                pulse period in seconds
* dm : float
            dispersion measure in pc cm^-3
* dec : float
            declination in degrees
* dtd : float
            scintillation timescale in seconds
* dnud : float
            scintillation bandwidth in GHz
* taud : float
            scattering timescale in us
* dist : float
            DM distance in kpc
* w50 : float
            FWHM of L-band pulse profile in us
* weff : float
            Effective width of the L-band profile in us
* uscale : float
                scaling factor to distribute intensity across
                pulse profile
* s_1000 : float
                flux density at 1 GHz in mJy
* spindex : float
                spectral index
* sig_j_single : float
                    single-pulse RMS jitter in us
* sigmas : dict
                dictionary of dictionaries of RMS components
                for each instrument


### Methods

sigma_jitter(self, t_int):
    Return intrinsic jitter noise (in us)
    for given integration time in seconds

add_sigmas(self, instr_name, sigma_tup):
    instr_name : str
                 name of timing instrument
    sigma_tup : tuple
                tuple of RMS components"""

get_instr_keys(self):
    returns instrument key names for sigmas dict