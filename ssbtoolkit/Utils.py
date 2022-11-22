
__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"

import os, warnings
from glob import glob
import numpy as np
import pandas as pd
warnings.simplefilter(action='ignore')
 

def FastaSequence(uniprotID):
    from bioservices import UniProt
    u = UniProt(verbose=False)

    fasta = u.search(uniprotID, frmt='fasta', limit=None)
    return "".join(fasta.split("\n")[1:])

def PrintProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
            
    :parameter iteration:Required: current iteration (Int)
    :parameter total:    Required: total iterations (Int)
    :parameter prefix:   Optional: prefix string (Str)
    :parameter suffix:   Optional: suffix string (Str)
    :parameter decimals: Optional: positive number of decimals in percent complete (Int)
    :parameter length:   Optional: character length of bar (Int)
    :parameter fill:     Optional: bar fill character (Str)
    :parameter printEnd: Optional: end character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

def CalcOccupancy(receptor_conc, agonist_conc, antagonist_conc, pkd_agonist, pkd_antagonist):
    """
    This function calculates the fraction of occupited receptors at equilibrium.
    
    :parameter receptor_conc: Required (int): concentration of the receptor
    :parameter agonists_conc: Required (int): concentration of the agonist
    :parameter antagonists_conc: Required (int): concentration of the antagonists (0 if antagonist shoul not be considered)
    :parameter pkd_agonist: Required (int): pKd of agonist
    :parameter pkd_antagonist: Required (int): pKd of antagonists (if antagonist shoul not be considered)
    :return int: fraction of occupied receptors in the equilibrium
    
    """
    import math
    #from pKd to Kd
    kd_anatagonist = (10**(-(float(pkd_antagonist)))) * 10**6

    if pkd_antagonist == 0:
        kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)
    else:
        kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)*(1+(antagonist_conc/kd_anatagonist))

    #LR determination
    a = 1
    b = float(agonist_conc)+float(receptor_conc)+kd_agonist
    c = float(receptor_conc)*float(agonist_conc)
    delta = (b**2) - (4*a*c)
    LR = (b-math.sqrt(delta))/(2*a)
    return LR

def MicrogramsToNanomolar(uniprotID, concentration):
    """
    This function converts micrograms of protein in nanomolar. 
    
    :parameter uniprotID: Required (str)
    :parameter concentration: Required (int): concentration od protein in micrograms
    :return: (flt) concentration of protein in nM
    
    .. note:: This function will obtain the sequence of the protein from UNIPROT and calculate automatically its molecular mass
    """

    from scipy.constants import Avogadro, micro, nano 
    #from Bio.Seq import Seq
    from Bio.SeqUtils import molecular_weight
    

    def Da2gr(x):
        return x*1.6605300000013E-24
    
    #Get protein sequence from UniProt
    seq = FastaSequence(uniprotID)
    
    #Prot MW
    prot_Da = molecular_weight(seq, "protein") #Da (daltons)
    
    #convert dalton to gram and then to picograms
    prot_gr = Da2gr(prot_Da)

    #convert grams to picograms
    prot_microgr = prot_gr/micro
    prot_microgr

    #calculate protein units (experimental)
    prot_Na = concentration/prot_microgr

    #Molar concentration of prot (experimental)
    prot_M = prot_Na/Avogadro
    prot_nM = prot_M/nano
    
    return round(prot_nM,4)

def KineticTempScale(kon, koff, T1, T2, Tu='K', *kwargs):
    """
    This function rescales the kinetics constants to a specific temperature. 
    
    :parameter kon:  Required (flt): foward kinetic constant
    :parameter koff: Required (flt): reverse kinetic constant
    :parameter T1:   Required (flt): Initial temperature
    :parameter T2:   Required (flt): Final temperature
    :paramter Tu:    Optional (kwarg str): Temperature Units (kelvin='K', celsius='C')
    :return: (flt, flt)
    
    """
    from scipy.constants import R
    from numpy import log

    #convert temperature units
    if Tu == 'K':
        pass 

    elif Tu =='C':
        T1 = T1 + 273.15     
        T2 = T2 + 273.15
    else :raise TypeError("Temperature must me in Kelvin (Tu ='K') or Celsius (Tu='C')") 


    #calculate free energies at different temperatures
    DG1 = -R*T1*log(koff/kon)

    DG2 = -R*T2*log(koff/kon)

    """if the T2 is higer than T1 (upscale)"""
    if T2 > T1:

        sf = DG2/DG1

        skon = kon*sf
        skoff = koff*sf

    """if T2 is lower than T1 (downscale)"""
    if T2 < T1:

        sf = DG2/DG1

        skon = kon*sf
        skoff = koff*sf

    if T2 == T1:
        skon = kon
        skoff = koff

    return round(skon, 3), round(skoff,3)

def GetGProtein(uniprotID):
        """
        This function query the SSBtoolkit internal database to extract the G protein associated to GPCR. 
        
        .. warning:: it just works for Human GPCRS!
        
        :parameter uniprotID:  Required (str)
        :return: (str)
        """
        import sqlite3
        
        #query HuTRdb 
        dbpath = os.path.join(os.path.split(os.path.realpath(__file__))[0],'HuTRdb.sqlite3')
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        cur.execute("SELECT * FROM gpcr WHERE uniprotid=?", (uniprotID,))
        rows = cur.fetchall()
        conn.close()
        gprotein=rows[0][7]
        return gprotein

class tauRAMD:
    """
    Implementation of the tRAMD method by Kokh et al., 2018.
    """

    def __init__(self):
        self._files = None
        self._dt = 2e-6
        self._softwr = 'GROMACS'
        self._prefix = None
        self.citation = '''Estimation of Drug-Target Residence Times by τ-Random Acceleration Molecular Dynamics Simulations
                            Daria B. Kokh, Marta Amaral, Joerg Bomke, Ulrich Grädler, Djordje Musil, Hans-Peter Buchstaller, Matthias K. Dreyer, Matthias Frech, Maryse Lowinski, Francois Vallee, Marc Bianciotto, Alexey Rak, and Rebecca C. Wade
                            Journal of Chemical Theory and Computation 2018 14 (7), 3859-3869
                            DOI: 10.1021/acs.jctc.8b00230 '''
        
    def Run(self, **kwargs):
        """
        Calulates the residence time of a ligand from RAMD simualtions.

        :parameter prefix: Required (kwarg str): directory path of .dat files
        :parameter dt:     Optional (kwarg flt): MD simulations time step in ns (defaul is 2E-6)
        :parameter softwr: Optional (kwarg str): software used to perform RAMD simulations: NAMD, GROMACS (default)
        :return (str): residence time
        """

        from scipy.stats import norm

        def bootstrapp(t, rounds=50000):
            """
            This function makes part of implementation of the tRAMD method by Kokh et al., 2018.
            
            :parameter t: Required (int): time
            :parameter rounds: Optional (int): default 50000
            """
            import numpy as np
            max_shuffle = rounds
            alpha = 0.8
            sub_set = int(alpha*len(t))
            tau_bootstr = []
            for i in range(1,max_shuffle):
                # generate a sub-set
                np.random.shuffle(t)
                t_b = t[:sub_set]
                # find residence time from a sub-stet
                t_b_sorted_50 =(np.sort(t_b)[int(len(t_b)/2.0-0.5)]+np.sort(t_b)[int(len(t_b)/2)])/2.0
                tau_bootstr.append(t_b_sorted_50)
            return(tau_bootstr)

        if 'prefix' not in kwargs: raise TypeError("ERROR: prefix is missing")
        if 'dt' in kwargs: self._dt = kwargs.pop('dt')
        if 'softwr' in kwargs: self._softwr = kwargs.pop('softwr')
        self._prefix = kwargs.pop('prefix')
        self._files = glob(self._prefix+'*.dat')
        self._times_set = []
        
        #Get Data
        for t,d in enumerate(self._files):
            with open(d) as f:
                read_data = f.readlines()
            self._times = []
            for r in read_data:
                if self._softwr == "NAMD":
                    self._times.append(int(r[r.find("EXIT:")+6:r.find(">")-2]))   # if NAMD  was used to generate RAMD trajectories
                elif self._softwr == "GROMACS":
                    self._times.append(int(r[r.find("after")+6:r.find("steps")-1]))   # if Gromacs was used to generate RAMD trajectories
                else: raise TypeError("ERROR: sofware unknown. options: NAMD, GROMACS")
            self._times = np.asarray(self._times)*self._dt
            self._times_set.append(self._times)
        
        #Parse Data
        self._mue_set=[]
        RTrelatives=[]
        for t, times in enumerate(self._times_set):
            if len(times) < 0: raise TypeError('ERROR: empty time values')
            else:
                bt2 = bootstrapp(times, rounds=50000)
                mu, std = norm.fit(bt2)
                
                bins = len(times)
                times = np.asarray(times)
                hist, bin_edges = np.histogram(times,bins=bins)
                hist_center = []
                for i,b in enumerate(bin_edges):
                    if i > 0: hist_center.append((bin_edges[i-1]+bin_edges[i])/2.0)
                CD = np.cumsum(hist)/np.max(np.cumsum(hist))
                KS = np.round(np.max(np.abs(1-np.exp(-(np.asarray(hist_center))/mu) - CD)),2)
                RTrelatives.append([t+1,mu,std,KS])
                
                self._mue_set.append(np.round(mu,1))
                self._RTmean = np.round(np.mean(self._mue_set),2)
                self._RTstd = np.round(np.std(self._mue_set),2)
                self.RT = self._RTmean
        
        self.RTdataframe = pd.DataFrame(RTrelatives, columns=['Replica no.', 'Relative res. time', 'SD', 'KS test'])
        print("Residence time:", str(self._RTmean),'±',str(self._RTstd),'ns')
        return
    
    def PlotRTDistribuitons(self,save=False, filename=None):
        """
        Plots the residence time distributions

        :parameter save:     Optional (kwarg boolean): default False
        :parameter filename: Optional (kwarg str)
        """

        import pylab as plt
        from IPython.core.display import display, HTML
        display(HTML("<style>.container { width:90% !important; }</style>"))

        fig  = plt.figure(figsize = (12,8))
        meanpointprops = dict(linestyle='--', linewidth=1.5, color='firebrick')
        medianpointprops = dict(linestyle='-', linewidth=2.0, color='orange')
        plt.boxplot(self._times_set,showmeans=True, meanline=True,meanprops=meanpointprops,medianprops = medianpointprops, bootstrap=5000)
        ymin, ymax = plt.ylim()
        plt.grid(linestyle = '--',linewidth=0.5)
        plt.yticks(np.linspace(0,int(ymax),min(int(ymax)+1,11)), fontsize=9)
        plt.ylabel('residence time [ns]', fontsize=10)
        plt.title("Residence times for "+str(len(self._files))+" replicas, mean: "+str(self._RTmean)+"  std: "+str(self._RTstd),fontsize=12)
        if save==True:
            if filename==None: 
                filename='plot.png'
                plt.savefig(filename, dpi=300)
            else:
                ext = os.path.splitext(filename)[-1]
                if ext == '.png': plt.savefig(filename, dpi=300)
                else: raise TypeError("extension not valid. Use png")
        
        return
        
    def PlotRTStats(self,save=False, filename=None):
        """
        Plots the residence time statistics

        :parameter save:     Optional (kwarg boolean): default False
        :parameter filename: Optional (kwarg str)
        """

        from matplotlib import gridspec
        from scipy.stats import norm

        import pylab as plt
        from IPython.core.display import display, HTML
        display(HTML("<style>.container { width:90% !important; }</style>"))


        def ret_time(t):
            """
            This function makes part of implementation of the tRAMD method by Kokh et al., 2018.
            
            :parameter t: Required (int): time   
            """
            import numpy as np
            t_sorted_50 =(np.sort(t)[int(len(t)/2.0-0.5)]+np.sort(t)[int(len(t)/2)])/2.0
            tau = t_sorted_50
            return(tau)

        def bootstrapp(t, rounds=50000):
            """
            This function makes part of implementation of the tRAMD method by Kokh et al., 2018.
            
            :parameter t: Required (int): time
            :parameter rounds: Optional (int): default 50000
            """
            import numpy as np
            max_shuffle = rounds
            alpha = 0.8
            sub_set = int(alpha*len(t))
            tau_bootstr = []
            for i in range(1,max_shuffle):
                # generate a sub-set
                np.random.shuffle(t)
                t_b = t[:sub_set]
                # find residence time from a sub-stet
                t_b_sorted_50 =(np.sort(t_b)[int(len(t_b)/2.0-0.5)]+np.sort(t_b)[int(len(t_b)/2)])/2.0
                tau_bootstr.append(t_b_sorted_50)
            return(tau_bootstr)

        fig  = plt.figure(figsize = (12,10))
        gs = gridspec.GridSpec(nrows=3, ncols=len(self._files), wspace=0.3,hspace=0.6)
        
        mue_set=[]
        for t, times in enumerate(self._times_set):

            if len(times) > 0:
                
            #First Row
                ax0 = fig.add_subplot(gs[0, t])

                #histogram
                bins = int(len(times)/2)
                s = ax0.hist(times,bins=bins,cumulative=True,histtype="step",color='k',lw=1)

                #plot redline at 50% of data
                ax0.plot([min(times), max(times)],[len(times)/2,len(times)/2], color='red', alpha = 0.5, linestyle='dashed',)
                plt.title("raw CDF",fontsize=12)
                ax0.set_xlabel('dissociation time [ns]', fontsize=10)
                tau = ret_time(times)
                ax0.plot([tau,tau],[0,len(times)/2.0], color='red', alpha = 0.5)
            
            #Second Row
                bt2 = bootstrapp(times, rounds=50000)
                bins = 6
                ax1 = fig.add_subplot(gs[1, t])
                ax1.hist(x=bt2,bins=bins, alpha=0.8,density=True,histtype="step")
                mu, std = norm.fit(bt2)
                mue_set.append(np.round(mu,1))
                xmin, xmax = plt.xlim()
                ymin, ymax = plt.ylim()
                x = np.linspace(0.8*xmin, xmax, 100)
                p = norm.pdf(x, mu, std)
                ax1.plot(x, p, 'k', linewidth=2)
                ax1.plot([mu,mu],[0, max(p)], color='red', alpha = 0.5)
                ax1.plot([xmin, xmax],[max(p)/2.0,max(p)/2.0], color='red', alpha = 0.5)
                ax1.plot([0.8*xmin, mu],[max(p),max(p)], color='red', linestyle='dashed',alpha = 0.5)
                ax1.set_xlabel('res. time [ns]', fontsize=10)
                plt.title("tau distribution",fontsize=12)
                ax1.set_yticks([])
                
            #Third Row
            ax2 = fig.add_subplot(gs[2, t])
            xmin = min(times)
            xmax = np.round(max(times))
            if(xmax==0): xmax = 0.5
            tp = np.linspace(xmin*0.5,xmax*1.5,100)
            poisson = 1-np.exp(-tp/mu) #np.cumsum(1-np.exp(-np.linspace(xmin,xmax,10)/mu))
            points=len(times)
            bins = len(times)
            times = np.asarray(times)
            hist, bin_edges = np.histogram(times,bins=bins)
            hist_center = []
            for i,b in enumerate(bin_edges):
                if i > 0: hist_center.append((bin_edges[i-1]+bin_edges[i])/2.0)
            CD = np.cumsum(hist)/np.max(np.cumsum(hist))
            ax2.scatter(np.log10(np.asarray(hist_center)),CD,marker='o')
            ax2.set_xlabel('log(res. time [ns])', fontsize=10)
            ax2.plot(np.log10(tp),poisson,color = 'k')
            ax2.set_ylim(0,1)
            ax2.set_xlim(-1.5,1.5)
            ax2.set_yticks(np.linspace(0,1,5))
            if (t> 0): ax2.set_yticklabels( [])
            plt.grid(linestyle = '--',linewidth=0.5)
            ax2.plot([np.log10(mu),np.log10(mu)],[0, 1], color='red', alpha = 0.5)
            KS = np.round(np.max(np.abs(1-np.exp(-(np.asarray(hist_center))/mu) - CD)),2)
            plt.title("KS test:"+str(KS),fontsize=12)
        
        
        if save==True:
            if filename==None: 
                filename='plot.png'
                plt.savefig(filename, dpi=300)
            else:
                ext = os.path.splitext(filename)[-1]
                if ext == '.png': plt.savefig(filename, dpi=300)
                else: raise TypeError("extension not valid. Use png")
        
        return

