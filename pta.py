class PTA(object):
    """
    Class to store all timed pulsars

    Attributes:
    __________
    name : string (optional)
           PTA name
    psrlist : list
              list of pulsar.Pulsar objects
        """
    def __init__(self,
                 name=None,
                 psrlist=None,
                 *args,
                 **kwargs):
        """
        ___init___ function for the Pulsar class
        """

        self.name = name
        self.psrlist = []

    def get_single_pulsar(self, psr_name):
        """return pulsar.Pulsar object whose name matches psr_name
        if no match, returns None"""
        for p in self.psrlist:
            if p.name == psr_name:
                return p
        else:
            print("No pulsar named {} in PTA".format(psr_name))
            return

