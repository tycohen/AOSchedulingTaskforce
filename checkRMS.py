import os
import cPickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class RMSPlotter(object):
    """
    Main Plotting Class
    
    Attributes
    ----------
    ptafile : str
            path to pickled .pta file (Default ./NG15yr.pta)
    """
    
    def __init__(self,
                 ptafile='NG15yr.pta',
                 psrnames=None,
                 fig=None,
                 ax=None,
                 annot=None,
                 scatter=None):
        """__init__ function for plotter"""
        try:
            with open(ptafile, "rb") as f:
                self.pta = cPickle.load(f)
        except IOError as e:
            print("ptafile does not exist.")
            raise
        self.ptafile = ptafile
        self.psrnames = psrnames
        self.fig = fig
        self.ax = ax
        self.annot = annot
        self.scatter = scatter

    def check_rcvr_keys(self, x_rcvr, y_rcvr):
        """Check if receiver/telescope system in .pta file"""
        instr_keys = self.pta.psrlist[0].get_instr_keys()
        if not x_rcvr in instr_keys:
            raise ValueError("x_rcvr={} not a valid instr_key in {}.\n"
                             "Valid keys are:\n{}.".format(x_rcvr,
                                                         self.ptafile,
                                                           "\n".join(instr_keys)))
        if not y_rcvr in instr_keys:
            raise ValueError("y_rcvr={} not a valid instr_key in {}.\n"
                             "Valid keys are:\n{}".format(y_rcvr,
                                                          self.ptafile,
                                                          "\n".join(instr_keys)))
        return

    def plot_vs(self, x_rcvr, y_rcvr, sigma_key="sigma_tot", save=False):
        self.check_rcvr_keys(x_rcvr, y_rcvr)
        detected = [p for p in self.pta.psrlist if p.detected(only=[x_rcvr, y_rcvr])]
        self.psrnames = [p.name for p in detected]
        x_tot_rms = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected])
        y_tot_rms = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected])
        j1713x = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])
        j1713y = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])

        self.fig = plt.figure(figsize=(12, 8))
        self.ax = plt.subplot(111)
        plt.ion()
        plt.show()
        self.ax.set_title("{} vs {}".format(x_rcvr.strip("_logain"),
                                       y_rcvr.strip("_logain")))
        self.ax.set_ylabel("{} / {}".format(x_rcvr.strip("_logain"),
                                       y_rcvr.strip("_logain")))
        self.ax.set_xlabel("{} (us)".format(x_rcvr.strip("_logain")))
        self.ax.axhline(1., color="black", ls=":")
        self.ax.set_xscale('log')
        self.ax.set_ylim([0, 5])
        self.scatter = self.ax.scatter(x_tot_rms, x_tot_rms / y_tot_rms,
                                       marker='x')
        self.ax.plot(j1713x, j1713x / j1713y,
                     'o', color="orange",
                     markersize=12, fillstyle="none",
                     label="J1713+0747")
        self.ax.legend(loc="upper right")
        self.ax.grid()
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20, 20),
                            textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"))
        self.annot.set_visible(False)

        print("Above dashed line is better with {}".format(y_rcvr))

        if save:
            filename = "{}_{}_vs_{}_{}.png".format(self.ptafile.strip(".pta"),
                                                   x_rcvr.strip("_logain"),
                                                   y_rcvr.strip("_logain"),
                                                   sigma_key)
            if not os.path.exists("plots"):
                warnings.warn("./plots dir does not exist. Creating.")
                os.mkdir("./plots")
            savepath = os.path.join("./plots", filename)
            print("Saved ./plots/{}".format(filename))
            plt.savefig(savepath)
            
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        plt.draw()
        plt.pause(0.001)
        return
    
    def plot_vs_dec(self, x_rcvr, y_rcvr, sigma_key='sigma_tot',
                    save=False):
        self.check_rcvr_keys(x_rcvr, y_rcvr)
        detected = [p for p in self.pta.psrlist if p.detected(only=[x_rcvr, y_rcvr])]
        self.psrnames = [p.name for p in detected]
        x_tot_rms = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected])
        y_tot_rms = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected])
        decs = np.array([p.dec for p in detected])
        j1713x = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])
        j1713y = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])
        j1713_dec = np.array([p.dec for p in detected
                                  if p.name.startswith("J1713+0747")])

        self.fig = plt.figure(figsize=(12, 8))
        self.ax = plt.subplot(111)
        plt.ion()
        plt.show()
        self.ax.set_title("{} & {} ({}) vs dec".format(x_rcvr.strip("_logain"),
                                                       y_rcvr.strip("_logain"),
                                                       sigma_key))
        self.ax.set_ylabel(r"$\sigma$ ({}) / $\sigma$ ({})".format(x_rcvr.strip("_logain"),
                                                              y_rcvr.strip("_logain")))
        self.ax.set_xlabel(r"$\mathrm{dec\ (deg)}$")
        self.ax.axhline(1., color="black", ls=":")
#        self.ax.set_ylim([0, 5])
        self.scatter = self.ax.scatter(decs, x_tot_rms / y_tot_rms,
                                       marker='x')
        self.ax.plot(j1713_dec, j1713x / j1713y,
                     'o', markersize=12, color="orange",
                     fillstyle="none", label="J1713+0747")
        self.ax.legend(loc="upper right")
        self.ax.grid()
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20, 20),
                                      textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"))
        self.annot.set_visible(False)

        print("Above dashed line is better with {}".format(y_rcvr))

        if save:
            if not os.path.exists("plots"):
                warnings.warn("./plots dir does not exist. Creating.")
                os.mkdir("./plots")
            filename = "{}_{}_and_{}_{}_vs_dec.png".format(self.ptafile.strip(".pta"),
                                                           x_rcvr.strip("_logain"),
                                                           y_rcvr.strip("_logain"),
                                                           sigma_key)
            savepath = os.path.join("./plots", filename)
            plt.savefig(savepath)
            print("Saved ./plots/{}".format(filename))
            
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        plt.draw()
        plt.pause(0.001)
        return        

    def plot_vs_scattering(self, x_rcvr, y_rcvr, sigma_key='sigma_tot',
                           save=False):
        self.check_rcvr_keys(x_rcvr, y_rcvr)
        detected = [p for p in self.pta.psrlist if p.detected(only=[x_rcvr, y_rcvr])]
        self.psrnames = [p.name for p in detected]
        x_tot_rms = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected])
        y_tot_rms = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected])
        t_scatter = np.array([p.taud for p in detected])
        j1713x = np.array([p.sigmas[x_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])
        j1713y = np.array([p.sigmas[y_rcvr][sigma_key] for p in detected
                           if p.name.startswith("J1713+0747")])
        j1713_scatter = np.array([p.taud for p in detected
                                  if p.name.startswith("J1713+0747")])

        self.fig = plt.figure(figsize=(12, 8))
        self.ax = plt.subplot(111)
        plt.ion()
        plt.show()
        self.ax.set_title("{} & {} vs scattering time".format(x_rcvr.strip("_logain"),
                                                         y_rcvr.strip("_logain")))
        self.ax.set_ylabel(r"$\sigma$ ({}) / $\sigma$ ({})".format(x_rcvr.strip("_logain"),
                                                              y_rcvr.strip("_logain")))
        self.ax.set_xlabel(r"$\tau_\mathrm{scatter}\ \mathrm{(us)}$")
        self.ax.axhline(1., color="black", ls=":")
        self.ax.set_xscale('log')
        self.ax.set_ylim([0, 5])
        self.scatter = self.ax.scatter(t_scatter, x_tot_rms / y_tot_rms,
                                       marker='x')
        self.ax.plot(j1713_scatter, j1713x / j1713y,
                     'o', markersize=12, color="orange",
                     fillstyle="none", label="J1713+0747")
        self.ax.legend(loc="upper right")
        self.ax.grid()
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20, 20),
                                      textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"))
        self.annot.set_visible(False)

        print("Above dashed line is better with {}".format(y_rcvr))

        if save:
            if not os.path.exists("plots"):
                warnings.warn("./plots dir does not exist. Creating.")
                os.mkdir("./plots")
            filename = "{}_{}_and_{}_{}_vs_scatter.png".format(self.ptafile.strip(".pta"),
                                                               x_rcvr.strip("_logain"),
                                                               y_rcvr.strip("_logain"),
                                                               sigma_key)
            savepath = os.path.join("./plots", filename)
            plt.savefig(savepath)
            
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        plt.draw()
        plt.pause(0.001)
        return        

    
    def update_annot(self, ind):
        pos = self.scatter.get_offsets()[ind["ind"][0]]
        self.annot.xy = pos
        text = "{}".format(", ".join([self.psrnames[n] for n in ind["ind"]]))
        self.annot.set_text(text)
        #self.annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        self.annot.get_bbox_patch().set_alpha(0.7)


    def hover(self, event):
        """from https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib"""
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            cont, ind = self.scatter.contains(event)
            if cont:
                self.update_annot(ind)
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()

def check_best_low_dec_CHIME():
    """
    Check which telescope low-declination pulsars moved to
    after correcting CHIME latitude from 39 to 49 deg
    """
    with open('NG15yr.pta', 'rb') as ptaf:
   	pta = cPickle.load(ptaf)
    with open('NG15yr_fixCHIMElat.pta', 'rb') as ptaf:
        fixCHIMEpta = cPickle.load(ptaf)
    print("Name     Best After Fix     RMS(us)        Best Before      RMS(us)")
    for p,t in zip(pta.sigma_2best(exclude=["AO", "UWBR"]),
                   fixCHIMEpta.sigma_2best(exclude=["AO", "UWBR"])):
        if p[1] != t[1]:
            print(t[0],
                  t[1],
                  "{:.4f}".format(t[2]),
                  p[1], "{:.4f}".format(p[2]),
                  pta.get_single_pulsar(t[0]).dec)
    return
