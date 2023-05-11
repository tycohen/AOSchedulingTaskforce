import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

freq, trx, gain = np.loadtxt("uwbr_rxspecs/November_2023_specs_nozap.txt", unpack=True) 
xfreq, xtrx, xgain = np.loadtxt("uwbr_rxspecs/Predicted_Specs_With_Excess_Improved_Gain.txt", unpack=True)
noxfreq, noxtrx, noxgain = np.loadtxt("uwbr_rxspecs/Predicted_Specs_No_Excess_Improved_Gain.txt", unpack=True)

def plot_current_v_target(save=False):
    fig, ax = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    plt.subplots_adjust(hspace=0)
    mpl.rcParams.update({'font.size': 14})
    ax[0].set_title("UWBR rxspecs current v target")
    ax[0].set_ylabel("gain (K/Jy)")
    ax[0].grid(ls=":")
    ax[0].plot(freq, gain, label="current")
    ax[0].plot(xfreq, xgain, ls="--", label="predicted, excess noise")
    ax[0].plot(noxfreq, noxgain, ls=":", label="target, no excess noise")  
    ax[0].legend(loc="best")

    ax[1].set_xlabel("freq (GHz)")
    ax[1].set_ylabel("T_rx (K)")
    ax[1].grid(ls=":")
    ax[1].plot(freq, trx, label="current")                 
    ax[1].plot(xfreq, xtrx, ls="--", label="predicted, excess noise")
    ax[1].plot(noxfreq, noxtrx, ls=":", label="target, no excess noise")
    if save:
        plt.savefig("./plots/UWBR_rxspecs_current_v_target.pdf")
        plt.savefig("./plots/UWBR_rxspecs_current_v_target.png")
    plt.show()
    
    return

def prepare_rxspecfiles():
    chime_str = """
# CHIME Gains multiplied by sqrt(.7) to account for BW loss to RFI
# Bands overlap from ~700 - 800 MHz, UWBR figures chosen in range for lower SEFD
# Integration times (s) for pulsar with dec=0
#freq	Trx	G	eps	t_int
.400	30.	0.97	0.01	1.8e4
0.699	30.	0.97	0.01	1.03e4
"""
    
    fstr = ""
    with open("uwbr_rxspecs/CHIME-GBTUWBR_WithExcNoise.txt", "wb") as outf:
        outf.write("# Hacky CHIME F-Engine + GBT UWBR pessimistic"
                   " with excess rcvr noise" + chime_str)
        for f,t,g in zip(xfreq, xtrx, xgain):
            fstr += "{}\t{}\t{}\t0.01\t3600.\n".format(f / 1000.,t,g)
        outf.write(fstr)

    fstr = ""
    with open("uwbr_rxspecs/CHIME-GBTUWBR_NoExcNoise.txt", "wb") as outf:
        outf.write("# Hacky CHIME F-Engine + GBT UWBR with no "
                   "excess receiver noise" + chime_str)
        for f,t,g in zip(noxfreq, noxtrx, noxgain):
            fstr += "{}\t{}\t{}\t0.01\t3600.\n".format(f / 1000.,t,g)
        outf.write(fstr)
        
    fstr = ""
    with open("uwbr_rxspecs/GBTUWBR_WithExcNoise.txt", "wb") as outf:
        outf.write("""# GBT UWBR pessimistic with excess receiver noise - Ryan Lynch 4/19/23
#freq	Trx	G	eps\n""")
        for f,t,g in zip(xfreq, xtrx, xgain):
            fstr += "{}\t{}\t{}\t0.01\n".format(f / 1000.,t,g)       
        outf.write(fstr)
        
    fstr = ""
    with open("uwbr_rxspecs/GBTUWBR_NoExcNoise.txt", "wb") as outf:
        outf.write("""# GBT UWBR with no excess receiver noise - Ryan Lynch 4/19/23
#freq	Trx	G	eps\n""")
        for f,t,g in zip(noxfreq, noxtrx, noxgain):
            fstr += "{}\t{}\t{}\t0.01\n".format(f / 1000.,t,g)       
        outf.write(fstr)

    return
