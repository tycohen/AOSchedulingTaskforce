import numpy as np

class Pulsar(object):
    """
    Class to store an individual 15-yr pulsar

    Attributes:
    __________
    name : string
            pulsar name
    period : float
                pulse period in seconds
    dm : float
            dispersion measure in pc cm^-3
    dec : float
            declination in degrees
    dtd : float
            scintillation timescale in seconds
    dnud : float
            scintillation bandwidth in GHz
    taud : float
            scattering timescale in us
    dist : float
            DM distance in kpc
    w50 : float
            FWHM of L-band pulse profile in us
    weff : float
            Effective width of the L-band profile in us
    uscale : float
                scaling factor to distribute intensity across
                pulse profile
    s_1000 : float
                flux density at 1 GHz in mJy
    spindex : float
                spectral index
    sig_j_single : float
                    single-pulse RMS jitter
    sigmas : dict
                dictionary of dictionaries of RMS components
                for each instrument
        """
    def __init__(self,
                 name=None,
                 period=None,
                 dm=None,
                 dec=None,
                 dtd=None,
                 dnud=None,
                 taud=None,
                 dist=None,
                 w50=None,
                 weff=None,
                 uscale=None,
                 s_1000=None,
                 spindex=None,
                 sig_j_single=None,
                 sigmas=None,
                 *args,
                 **kwargs):
        """
        ___init___ function for the Pulsar class
        """

        self.name = name
        self.period = period
        self.dm = dm
        self.dec = dec
        self.dtd = dtd
        self.dnud = dnud
        self.taud = taud
        self.dist = dist
        self.w50 = w50
        self.weff = weff
        self.uscale = uscale
        self.s_1000 = s_1000
        self.spindex = spindex
        self.sig_j_single = sig_j_single
        self.sigmas = {}

    def sigma_jitter(self, t_int):
        """Return intrinsic jitter noise (in us)
        for given integration time in seconds"""
        n_pulses = t_int / self.period
        return self.sig_j_single / np.sqrt(n_pulses)            

    def add_sigmas(self, instr_name, sigma_tup):
        """instr_name : str
                        name of timing instrument
        sigma_tup : tuple
                        tuple of RMS components"""
        keys = ['sigma_tot',
                'sigma_white',
                'sigma_dm',
                'sigma_tel',
                'sigma_rn']
        self.sigmas.update({instr_name: dict(zip(keys, sigma_tup))})

    def get_instr_keys(self):
        return [k for k in self.sigmas]

    def detected(self, only=None):
        """
        Is the pulsar is detected by every telescope?
        
        Parameters
        ----------
        only : None or list (optional)
                list of receiver keys to only check if not None

        Returns
        -------
        bool
        """
        if only is None:
            rcvrs = self.get_instr_keys()
        else:
            rcvrs = only
        return all([self.sigmas[r]["sigma_tot"] > 0 for r in rcvrs])
        
