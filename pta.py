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
                 psrlist=[],
                 *args,
                 **kwargs):
        """
        ___init___ function for the Pulsar class
        """

        self.name = name
        self.psrlist = psrlist

    def get_single_pulsar(self, psr_name):
        """return pulsar.Pulsar object whose name matches 'psr_name'
        """
        for p in self.psrlist:
            if p.name == psr_name:
                return p
        else:
            raise ValueError("No pulsar named {} in PTA".format(psr_name))

    def sigma_best(self, exclude=[]):
        """
        Get the best instrument for each pulsar
        and return list of tuples of (pulsar name, instrument, sigma_tot)
        
        Parameters
        ----------
        exclude : list
                  list of telescope name substrings to exclude
        """
        if not isinstance(exclude, list):
            raise TypeError("'exclude' must be a list of substrings not "
                            "{}".format(type(exclude)))
        best_instr_list = []
        for p in self.psrlist:
            best_tup = sorted([(p.name, k.strip('_logain'), v['sigma_tot'])
                               for k, v in p.sigmas.iteritems()
                               if not any([e in k for e in exclude])],
                              key=lambda t: (t[2] < 0., t[2]))[0]
            best_instr_list.append(best_tup)
        return best_instr_list

    # def sigma_2ndbest(self, exclude='*'):
    #     """
    #     Get the best and 2nd instrument for each pulsar
    #     and return list of tuples of (pulsar name, instrument, % diff)
    #     Set exclude = string to ignore a particular telescope
    #     """
    #     best_instr_list = []
    #     for p in self.psrlist:
    #         sigmas_sorted = sorted([(p.name, k.strip('_logain'), v['sigma_tot'])
    #                                 for k, v in p.sigmas.iteritems()
    #                                 if not exclude in k],
    #                                key=lambda t: (t[2] < 0., t[2]))
    #         best = sigmas_sorted[0]
    #         secbest = sigmas_sorted[1]
            
    #         best_instr_list.append(best_tup)
    #     return best_instr_list
    

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
