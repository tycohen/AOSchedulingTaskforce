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
        """return pulsar.Pulsar object whose name matches 'psr_name'
        if no match, returns None"""
        for p in self.psrlist:
            if p.name == psr_name:
                return p
        else:
            raise ValueError("No pulsar named {} in PTA".format(psr_name))

    def write_to_txt(self, filename):
        """Write total RMS for each pulsar at each instrument to file"""
        key_names = sorted([k for k in self.psrlist[0].get_instr_keys()])
        instr_names = sorted([k.strip("_logain") for k in key_names])
        lines = []
        header = "\t".join(["# name"] + instr_names)
        lines.append(header)
        for p in self.psrlist:
             s = "\t".join([p.name] + [str(p.sigmas[k]['sigma_tot']) for k in key_names])
             lines.append(s)
        all_lines = "\n".join(lines)
        with open(filename, 'wb') as f:
            f.write(all_lines)
        return
