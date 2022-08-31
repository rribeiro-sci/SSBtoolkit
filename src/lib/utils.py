
__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"

class utils:
    def FastaSequence(uniprotID):
        from bioservices import UniProt
        u = UniProt(verbose=False)

        fasta = u.search(uniprotID, frmt='fasta', limit=None)
        return "".join(fasta.split("\n")[1:])

    def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
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

    def LR_eq_conc(receptor_conc, agonist_conc, antagonist_conc, pkd_agonist, pkd_antagonist):
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

    def maxbend(drug_receptor, lig_conc_range):
        """
        This function calculates the maximum bending point of a sigmoid-shaped curve according to the mathod of Sebaugh et al., 2003.
        
        :parameter drug_receptor: Required (int): concentration of the receptor
        :parameter lig_conc_range: Required (array): array of a range of ligand concentration
        :return: int
        
        The minimization uses the Nelder-Mead method.
        """


        from scipy.optimize import curve_fit, minimize
        import numpy as np

        def equation_binding(X, Bottom, Top, Kd, p):
            return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

        xfit = np.geomspace(1E-3, 1E4, 50000)
        popt, pcov = curve_fit(equation_binding, lig_conc_range, drug_receptor, bounds=([np.min(drug_receptor),-np.inf,-np.inf, 0.5],[np.inf,np.max(drug_receptor),np.inf, 2.5]))


        def equation_binding_deriv_b(x, a,d,c,b):
            return (x/c)**b*(a - d)*np.log(x/c)/((x/c)**b + 1)**2

        min_value = minimize(equation_binding_deriv_b, np.max(xfit), args=(popt[0],popt[1],popt[2],popt[3]), method = 'Nelder-Mead')

        submaximal = round(min_value.x[0],3)
        return submaximal

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

    def ret_time(t):
        """
        This function makes part of implementation of the tRAMD method by Kokh et al., 2018.
        
        :parameter t: Required (int): time   
        """
        import numpy as np
        t_sorted_50 =(np.sort(t)[int(len(t)/2.0-0.5)]+np.sort(t)[int(len(t)/2)])/2.0
        tau = t_sorted_50
        return(tau)