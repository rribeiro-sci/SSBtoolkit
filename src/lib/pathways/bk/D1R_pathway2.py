##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                              D1 receptor                                     ##
##                                                                              ##
##                               PATHWAY 2                                      ##
##                                                                              ##
##                        last modification 10 April 2019                       ##
##                                                                              ##
##   References:                                                                ##
##   Nair, A. G. et al (2015). Sensing Positive versus Negative Reward Signals  ##
##   through Adenylyl Cyclase-Coupled GPCRs in Direct and Indirect Pathway      ##
##   Striatal Medium Spiny Neurons. The Journal of Neuroscience : the Official  ##
##   Journal of the Society for Neuroscience, 35(41), 14017–14030.              ##                                                          ##
##################################################################################
from pysb import *
from pysb.macros import *
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator
import math

def network(pKd,ligand,epsilon):


    Rtotal = 2000
    #from pKd to Kd
    Kd = (10**(-(float(pKd)))) * 10**9
    #LR determination
    a = 1
    b = float(ligand)+float(Rtotal)+Kd
    c = float(Rtotal)*float(ligand)
    delta = (b**2) - (4*a*c)
    LR = (b-math.sqrt(delta))/(2*a)

    #Start a model
    Model()

    #MONOMERS
    Monomer('D1R', ['D1R_b1', 'D1R_p', 'D1R_s'], {'D1R_p':['p0','p1'], 'D1R_s':['i','a']})    # Dopamine receptor 1 with two binding sites: for DA and G-proteins
    Monomer('Golf', ['Golf_b1'])            # Golf protein (alpha and beta/gamma complex)
    Monomer('Ach', ['Ach_b1'])
    Monomer('M4R', ['M4R_b1', 'M4R_b2'])
    Monomer('Gi', ['Gi_b1'])
    Monomer('Gbgolf')
    Monomer('GaolfGTP', ['GaolfGTP_b1'])
    Monomer('GaolfGDP')
    Monomer('GaiGTP', ['GaiGTP_b1'])
    Monomer('Gbgi')
    Monomer('GaiGDP')
    Monomer('AC5', ['AC5_b1', 'AC5_b2', 'AC5_b3', 'AC5_s'],  {'AC5_s': ['i', 'a']})
    Monomer('Ca', ['Ca_b1', 'Ca_l', 'Ca_s'], {'Ca_l':['cytos', 'ext'], 'Ca_s':['free' , 'buff']}) # l (location), s (state: free or buffed)
    Monomer('ATP')
    Monomer('cAMP')
    Monomer('PDE4', ['PDE4_s'], {'PDE4_s':['i', 'a']})
    Monomer('AMP')
    Monomer('PDE10', ['PDE10_c', 'PDE10_s'], {'PDE10_c':['N', 'Y'], 'PDE10_s':['i', 'a']})
    Monomer('PKA', ['PKA_s'], {'PKA_s':['cAMP0', 'cAMP2', 'cAMP4']})
    Monomer('PKAc', ['PKAc_b1'])
    Monomer('PKAreg')
    Monomer('CaM', ['CaM_b1','CaM_s'], {'CaM_s':['Ca0', 'Ca2', 'Ca4']})
    Monomer('PP2B', ['PP2B_b1', 'PP2B_b2'])
    Monomer('B72PP2A', ['B72PP2A_b1', 'B72PP2A_b2'])
    Monomer('DARPP32', ['DARPP32_b1'])
    Monomer('D32p34', ['D32p34_b1'])
    Monomer('B56PP2A', ['B56PP2A_b1', 'B56PP2A_p'], {'B56PP2A_p':['N', 'Y']}) # B56PP2A_p(phospholation status)
    Monomer('PP1', ['PP1_b1'])
    Monomer('CDK5', ['CDK5_b1'])
    Monomer('D32p75', ['D32p75_b1'])
    Monomer('AKAR3', ['AKAR3_b1', 'AKAR3_p'], {'AKAR3_p':['N', 'Y']})

    #INITIAL CONDITIONS
    Initial(D1R(D1R_b1=None, D1R_p='p0', D1R_s='a'), Parameter('D1RDA_0', LR))   # nM
    Initial(D1R(D1R_b1=None, D1R_p='p0', D1R_s='i'), Parameter('D1R_0', 2000.0-LR))   # nM
    Initial(Golf(Golf_b1=None), Parameter('Golf_0', 2000.0))             # nM
    Initial(Ach(Ach_b1=None), Parameter('Ach_0', 100.0))                 # nM
    Initial(M4R(M4R_b1=None, M4R_b2=None), Parameter('M4R_0', 2000.0))   # nM
    Initial(Gi(Gi_b1=None), Parameter('Gi_0', 2000.0))                   # nM
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', 700))  # nM
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', 60)) # nM
    Initial(ATP(), Parameter('ATP_0', 5000000))
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', 2000))
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', 700))
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', 1200))
    Initial(CaM(CaM_b1=None, CaM_s='Ca0'), Parameter('CaM_0', 10000))
    Initial(PP2B(PP2B_b1=None, PP2B_b2=None), Parameter('PP2B_0', 2000))
    Initial(B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None), Parameter('B72PP2A_0', 2000))
    Initial(DARPP32(DARPP32_b1=None), Parameter('DARPP32_0', 50000))
    Initial(B56PP2A(B56PP2A_b1=None, B56PP2A_p='N'), Parameter('B56PP2A_0', 2000))
    Initial(PP1(PP1_b1=None), Parameter('PP1_0', 2000))
    Initial(CDK5(CDK5_b1=None), Parameter('CDK5_0', 1800))
    Initial(AKAR3(AKAR3_b1=None, AKAR3_p='N'), Parameter('AKA3_0', 3000))

    #PARAMETERS
    #Parameter('kD1R_DA_1', 0.005)       # 1/(nM*s) |Association constant of the complex D1R_DA
    #Parameter('kD1R_DA_2', 5.0)         # 1/s      |Dissociation constant of the complex D1R_DA
    Parameter('kD1R_Golf_1', 0.003)     # 1/(nM*s) |Association constant of the complex D1R_Golf
    Parameter('kD1R_Golf_2', 5.0)       # 1/s      |Dissociation constant of the complex D1R_Golf
    Parameter('kD1R_DA_Golf_1', 0.003)  # 1/(nM*s) |Association constant of the complex D1R_DA_Golf
    Parameter('kD1R_DA_Golf_2', 5.0)    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
    Parameter('kD1R_Golf_DA_1', 0.005)  # 1/(nM*s) |Association constant of the complex D1R_Golf_DA
    Parameter('kD1R_Golf_DA_2', 5.0)    # 1/s      |Dissociation constant of the complex D1R_Golf_DA
    Parameter('kM4R_Ach_1',0.01)        # 1/(nM*s) |Association constant of the complex M4R_Ach
    Parameter('kM4R_Ach_2',90)          # 1/s      |Dissociation constant of the complex M4R_Ach
    Parameter('kM4R_Gi_1',0.012)        # 1/(nM*s) |Association constant of the complex M4R_Gi
    Parameter('kM4R_Gi_2',90)           # 1/s      |Dissociation constant of the complex M4R_Gi
    Parameter('kM4R_Ach_Gi_1',1.2)      # 1/(nM*s) |Association constant of the complex M4R_Ach_Gi
    Parameter('kM4R_Ach_Gi_2',90)       # 1/s      |Dissociation constant of the complex M4R_Ach_Gi
    Parameter('kM4R_Gi_Ach_1',1.0)      # 1/(nM*s) |Association constant of the complex M4R_Gi_Ach
    Parameter('kM4R_Gi_Ach_2',90)       # 1/s      |Dissociation constant of the complex M4R_Gi_Ach
    Parameter('kD1R_DA_Golf_decay', 15) # 1/s      |Rate of formation of M4R_DA + Golf
    Parameter('kGaolfGTP_decay', 30)    # 1/s      |Rate of convertion of GaolfGTP into GaolfGDP
    Parameter('kGolf_formation', 100)   # 1/s      |Rate of formation of Golf
    Parameter('kM4R_DA_Gi_decay', 60)   # 1/s      |Rate of formation of M4R_DA + Gi
    Parameter('kGaiGTP_decay', 30)      # 1/s      |Rate of formation of GaiGDP
    Parameter('kGi_formation', 100)     # 1/s      |Rate of formation of Gi
    Parameter('kAC5_GaolfGTP_1', 0.2)   # 1/(nM*s) |Association constant of the complex AC5_GaolfGTP
    Parameter('kAC5_GaolfGTP_2', 0.1)   # 1/s      |Dissociation constant of the complex AC5_GaolfGTP
    Parameter('kAC5_Ca_1', 0.001)       # 1/(nM*s) |Association constant of the complex AC5_Ca
    Parameter('kAC5_Ca_2', 0.9)         # 1/s      |Dissociation constant of the complex AC5_Ca
    Parameter('kAC5_Ca_GaolfGTP_1', 0.2)# 1/(nM*s) |Association constant of the complex AC5_Ca_GaolfGTP
    Parameter('kAC5_Ca_GaolfGTP_2', 0.1)# 1/s      |Dissociation constant of the complex AC5_Ca_GaolfGTP
    Parameter('kAC5_GaolfGTP_ATP_1', 0.00105) # 1/(nM*s) |Association constant of the complex AC5_GaolfGTP-ATP
    Parameter('kAC5_GaolfGTP_ATP_2', 1)       # 1/s      |Dissociation constant of the complex AC5_GaolfGTP-ATP
    Parameter('kcAMP_formation', 20)    # 1/s      |Rate constant of convertion of  ATP into cAMP by AC5
    Parameter('kcAMP_reverse', 0.084)   # 1/s      |Rate constant of reverse cAMP to ATP
    Parameter('kAC5_ATP_1', 0.0001)     # 1/(nM*s) |Association constant of the complex AC5_ATP
    Parameter('kAC5_ATP_2', 1)          # 1/s      |Dissociation constant of the complex AC5_ATP
    Parameter('kAC5_basal', 1)          # 1/s      |Rate constant of basal formation of cAMP
    Parameter('kAC5_reverse_basal', 0.0004) # 1/s  |Reverse Rate constant of basal formation of cAMP
    Parameter('kAC5_Ca_ATP_1', 7.50E-5) # 1/(nM*s) |Association constant of the complex AC5_Ca_ATP
    Parameter('kAC5_Ca_ATP_2', 1)       # 1/s      |Dissociation constant of the complex AC5_Ca_ATP
    Parameter('kAC5_Ca_ATP_to_cAMP', 0.5 )  # 1/s      |Convertion rate of ATP into cAMP by AC5_Ca
    Parameter('kAC5_Ca_ATP_to_cAMP_reverse', 0.00015) # 1/s      |Reverse Convertion rate of ATP into cAMP by AC5_Ca
    Parameter('kAC5_GaolfGTP_ATP_3', 0.2)   # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_ATP_4', 0.1)   # 1/s      |
    Parameter('kAC5_GaolfGTP_decay', 0.2)   # 1/s      |
    Parameter('kAC5_GaolfGTP_ATP_decay', 0.2)# 1/s     |
    Parameter('kAC5_ATP_Ca_1', 0.001)       # 1/(nM*s) |
    Parameter('kAC5_ATP_Ca_2', 0.9)         # 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_ATP_1', 0.00055)   # 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_ATP_2', 1)   # 1/s      |
    Parameter('kAC5_CA_ATP_GaolfGTP_1', 0.2) # 1/(nM*s) |
    Parameter('kAC5_CA_ATP_GaolfGTP_2', 0.1) # 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP', 10)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse', 0.022) # 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_decay', 0.2) # 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_ATP_decay', 0.2) # 1/s      |
    Parameter('kAC5_GaiGTP_1', 50)          # 1/(nM*s) |
    Parameter('kAC5_GaiGTP_2', 5)           # 1/s      |
    Parameter('kAC5_GaiGTP_ATP_1', 6.25E-5) # 1/(nM*s) |
    Parameter('kAC5_GaiGTP_ATP_2', 1)       # 1/s      |
    Parameter('kAC5_ATP_GaiGTP_1', 50)      # 1/(nM*s) |
    Parameter('kAC5_ATP_GaiGTP_2', 5)       # 1/s      |
    Parameter('kAC5_GaiGTP_ATP_to_cAMP', 0.25) # 1/s      |
    Parameter('kAC5_GaiGTP_ATP_to_cAMP_reverse', 0.00105) # 1/s      |
    Parameter('kAC5_GaiGTP_decay', 30)       # 1/s      |
    Parameter('kAC5_GaiGTP_decay_2', 30)     # 1/s      |
    Parameter('kAC5_GaolfGTP_GaiGTP_1', 0.01) # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_2', 0.01) # 1/s      |
    Parameter('kAC5_GaolfGTP_GaiGTP_3', 0.002) # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_4', 0.01) # 1/s      |
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_1', 0.0003)# 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_2', 1)     # 1/s      |
    Parameter('kAC5_GaiGTP_GaolfGTP_ATP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_GaiGTP_GaolfGTP_ATP_2', 0.01)# 1/s      |
    Parameter('kAC5_GaolfGTP_ATP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_ATP_GaiGTP_2', 0.01)# 1/s      |
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP', 5)# 1/s      |
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse', 0.006) # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_1', 0.2)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_2', 30)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_3', 0.2)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_4', 30)    # 1/s      |
    Parameter('kAC5_Ca_GaiGTP_1', 50)   # 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_2', 5)    # 1/s      |
    Parameter('kAC5_Ca_GaiGTP_ATP_1', 5.63E-5)  # 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_ATP_2', 1)    # 1/s      |
    Parameter('kAC5_Ca_ATP_GaiGTP_1', 50)   # 1/(nM*s) |
    Parameter('kAC5_Ca_ATP_GaiGTP_2', 5)    # 1/s      |
    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP', 0.125) # 1/s      |
    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse', 2.81E-5)# 1/s      |
    Parameter('kAC5_Ca_GaiGTP_decay', 30)# 1/s      |
    Parameter('kAC5_Ca_GaiGTP_ATP_decay', 30)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_2', 0.01)# 1/s      |
    Parameter('kAC5_Ca_GaiGTP_GaolfGTP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_GaolfGTP_2',0.01)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_1', 0.000175)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_2', 1)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_ATP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_ATP_GaiGTP_2', 0.01)# 1/s      |
    Parameter('kAC5_Ca_GaiGTP_ATP_GaolfGTP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_ATP_GaolfGTP_2', 0.01)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP', 2.5)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse', 0.00175)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_decay_1', 30)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_decay_2', 0.2)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_1', 30)# 1/s      |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_2', 0.2)# 1/s      |
    Parameter('kPDE4_cAMP_1', 0.01)# 1/(nM*s) |
    Parameter('kPDE4_cAMP_2', 1)# 1/s      |
    Parameter('kPDE4_cAMP_to_AMP', 2)
    Parameter('kPDE10_2cAMP_1', 1.0E-6)# 1/(nM*s) |
    Parameter('kPDE10_2cAMP_2', 9)# 1/s      |
    Parameter('kPDE10_cAMP_1', 0.1)# 1/(nM*s) |
    Parameter('kPDE10_cAMP_2', 2)# 1/s      |
    Parameter('kPDE10_2cAMP_cAMP_1', 0.13)# 1/(nM*s) |
    Parameter('kPDE10_2cAMP_cAMP_2', 2)# 1/s      |
    Parameter('kPDE10_cAMP_decay', 3)# 1/s      |
    Parameter('kPDE10_2cAMP_cAMP_decay', 10)# 1/s      |
    Parameter('kPKA_cAMP2_1', 0.00026)# 1/(nM*s) |
    Parameter('kPKA_cAMP2_2', 1)# 1/s |
    Parameter('kPKA_cAMP4_1', 0.000346)# 1/(nM*s) |
    Parameter('kPKA_cAMP4_2', 1)# 1/s |
    Parameter('kPKA_activation', 10)# 1/(nM*s) |
    Parameter('kPKA_activation_reverse', 0.01)# 1/s |
    Parameter('kCaM_Ca2_1', 0.006)# 1/(nM*s) |
    Parameter('kCaM_Ca2_2', 9.1)# 1/s |
    Parameter('kCaM_Ca4_1', 0.1)# 1/(nM*s) |
    Parameter('kCaM_Ca4_2', 1000)# 1/s |
    Parameter('kCaM_PP2B_1', 1)# 1/(nM*s) |
    Parameter('kCaM_PP2B_2', 30)# 1/s |
    Parameter('kCaMCa4_PP2B_1', 1)# 1/(nM*s) |
    Parameter('kCaMCa4_PP2B_2', 3)# 1/s |
    Parameter('kCaMCa2_PP2B_1', 0.006)# 1/(nM*s) |
    Parameter('kCaMCa2_PP2B_2', 0.91)# 1/s |
    Parameter('kCaMCa4_PP2B_3', 0.1)# 1/(nM*s) |
    Parameter('kCaMCa4_PP2B_4', 10)# 1/s |
    Parameter('kCaMCa2_PP2B_3', 1)# 1/(nM*s) |
    Parameter('kCaMCa2_PP2B_4', 3)# 1/s |
    Parameter('kB72PP2A_Ca_1', 0.0001)# 1/(nM*s) |
    Parameter('kB72PP2A_Ca_2', 0.1)# 1/s |
    Parameter('kPKAc_DARPP32_1', 0.001)# 1/(nM*s) |
    Parameter('kPKAc_DARPP32_2', 8)# 1/s |
    Parameter('kPKAc_DARPP32_decay', 10)# 1/s |
    Parameter('kPKAc_B56PP2A_1', 0.001)# 1/(nM*s) |
    Parameter('kPKAc_B56PP2A_2', 0.3)# 1/s |
    Parameter('kD23p34_PP1_1', 0.04)# 1/(nM*s) |
    Parameter('kD23p34_PP1_2', 0.08)# 1/s |
    Parameter('kB56PP2A_phosphorylation', 0.2) # 1/s |
    Parameter('kB56PP2Ap_decay', 0.008) # 1/s |
    Parameter('kCDK5_DARPP32_1', 0.0009)# 1/(nM*s) |
    Parameter('kCDK5_DARPP32_2', 2)# 1/s |
    Parameter('kCDK5_DARPP32_decay', 3)# 1/s |
    Parameter('kPKAc_D32p75_1', 0.00037)# 1/(nM*s) |
    Parameter('kPKAc_D32p75_2', 1)# 1/s |
    Parameter('kB56PP2Ap_D32p75_1', 0.0015) # 1/(nM*s) |
    Parameter('kB56PP2Ap_D32p75_2', 10)# 1/s |
    Parameter('kB72PP2A_D32p75_1', 0.0008) # 1/(nM*s) |
    Parameter('kB72PP2A_D32p75_2', 6.4)# 1/s |
    Parameter('kB56PP2Ap_D32p75_decay', 8)# 1/s |
    Parameter('kB72PP2A_D32p75_decay', 1)# 1/s |
    Parameter('kB72PP2ACa_D32p75_1', 0.0015)# 1/(nM*s) |
    Parameter('kB72PP2ACa_D32p75_2', 10 )# 1/s |
    Parameter('kB72PP2ACa_D32p75_decay', 5)# 1/s |
    Parameter('kB56PP2A_D32p75_decay', 2.3)# 1/s |
    Parameter('kPP2BCaM_D23p34_1', 0.002)# 1/(nM*s) |
    Parameter('kPP2BCaM_D23p34_2', 1)# 1/s |
    Parameter('kPP2BCaM_D23p34_decay', 10)# 1/s |
    Parameter('kD32p75_B56PP2A_1', 0.0008)# 1/(nM*s) |
    Parameter('kD32p75_B56PP2A_2', 6.4)# 1/s |
    Parameter('kB72PP2ACa_D32p34_1', 0)# 1/(nM*s) |
    Parameter('kB72PP2ACa_D32p34_2', 8) # 1/s |
    Parameter('kB72PP2A_D32p34_1', 0.5)# 1/(nM*s) |
    Parameter('kB72PP2A_D32p34_2', 1) # 1/s |
    Parameter('kB72PP2ACa_D32p34_decay', 0.3)# 1/s |
    Parameter('kB72PP2A_D32p34_decay', 0.3)# 1/s
    Parameter('kAKAR3_PKAc_1', 0.003)# 1/(nM*s) |
    Parameter('kAKAR3_PKAc_2', 1)# 1/s
    Parameter('kAKAR3_PKAc_decay', 3) # 1/s
    Parameter('kAKAR3p_PP1_1', 1.0E-4)# 1/(nM*s) |
    Parameter('kAKAR3p_PP1_2', 10)# 1/s
    Parameter('kAKAR3p_PP1_decay', 1.2) # 1/s

    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The D1R can bind either the inactive G protein first, and then
    dopamine, or dopamine first and then and then the inactivate G protein.'''

    Rule('reaction3', D1R(D1R_b1=None, D1R_p='p0', D1R_s='a') + Golf(Golf_b1=None) | D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50), kD1R_DA_Golf_1, kD1R_DA_Golf_2)
    Rule('reaction9', D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50) >> D1R(D1R_b1=None, D1R_p='p0', D1R_s='a') + Gbgolf() + GaolfGTP(GaolfGTP_b1=None) , kD1R_DA_Golf_decay)
    Rule('reaction10', GaolfGTP(GaolfGTP_b1=None) >> GaolfGDP(), kGaolfGTP_decay)
    Rule('reaction11', GaolfGDP() + Gbgolf() >> Golf(Golf_b1=None), kGolf_formation)

    '''The M4R can bind either the inactive G protein first, and then
    acetylcholine, or acetylcholine first and then and then the inactivate G protein.'''

    Rule('reaction5', M4R(M4R_b1=None, M4R_b2=None) + Ach(Ach_b1=None) | M4R(M4R_b1=50, M4R_b2=None) % Ach(Ach_b1=50), kM4R_Ach_1, kM4R_Ach_2)
    Rule('reaction6', M4R(M4R_b1=None, M4R_b2=None) + Gi(Gi_b1=None) | M4R(M4R_b1=None, M4R_b2=20) % Gi(Gi_b1=20), kM4R_Gi_1, kM4R_Gi_2)
    Rule('reaction7', M4R(M4R_b1=50, M4R_b2=None) % Ach(Ach_b1=50) + Gi(Gi_b1=None) | M4R(M4R_b1=50, M4R_b2=20) % Ach(Ach_b1=50) % Gi(Gi_b1=20), kM4R_Ach_Gi_1, kM4R_Ach_Gi_2)
    Rule('reaction8', M4R(M4R_b1=None, M4R_b2=20) % Gi(Gi_b1=20) + Ach(Ach_b1=None) | M4R(M4R_b1=50, M4R_b2=20) % Ach(Ach_b1=50) % Gi(Gi_b1=20), kM4R_Gi_Ach_1, kM4R_Gi_Ach_2)
    Rule('reaction12', M4R(M4R_b1=50, M4R_b2=20) % Ach(Ach_b1=50) % Gi(Gi_b1=20) >> M4R(M4R_b1=50, M4R_b2=None) % Ach(Ach_b1=50) + GaiGTP(GaiGTP_b1=None) + Gbgi() , kM4R_DA_Gi_decay)
    Rule('reaction13', GaiGTP(GaiGTP_b1=None) >> GaiGDP(), kGaiGTP_decay)
    Rule('reaction14', Gbgi() + GaiGDP() >> Gi(Gi_b1=None), kGi_formation)

    '''CYCLIC AMP FORMATION AND PKA ACTIVATION:
    - active GaolfGTP binds to and activate cyclase type V (AC5), which produces cyclic AMP (cAMP)
    - Ca2+ causes a 50% reduction in the activity of AC5 before or after binding to Gaolf. In the model
        calcium is allowed to bind to the inactive AC5, which then can be activated by GaolfGTP.
    - A step of regeneration of ATP is included that allows for a steady state concentration of ATP in
    the absence of stimulation. It is not assumed to be a rate-limiting substrate.
    '''

    # AC5 basal activity
    Rule('reaction21', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') , kAC5_ATP_1, kAC5_ATP_2)
    Rule('reaction22', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP(), kAC5_basal)
    Rule('reaction23', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'), kAC5_reverse_basal )

    # Interaction between AC5 and Ca2+
    Rule('reaction16', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_1, kAC5_Ca_2)
    Rule('reaction24', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_1, kAC5_Ca_ATP_2)
    Rule('reaction25', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') >> cAMP() + AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP)
    Rule('reaction26', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP_reverse)
    Rule('reaction30', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') , kAC5_ATP_Ca_1, kAC5_ATP_Ca_2)

    # Interaction between AC5 and Gaolf
    Rule('reaction15', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_1, kAC5_GaolfGTP_2)
    Rule('reaction18', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_ATP_1, kAC5_GaolfGTP_ATP_2)
    Rule('reaction19', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , kcAMP_formation)
    Rule('reaction20', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kcAMP_reverse)
    Rule('reaction27', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_ATP_3, kAC5_GaolfGTP_ATP_4)
    Rule('reaction28', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGDP(), kAC5_GaolfGTP_decay)
    Rule('reaction29', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGDP(), kAC5_GaolfGTP_ATP_decay)

    # Interaction between AC5 and Gai
    Rule('reaction37', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) , kAC5_GaiGTP_1, kAC5_GaiGTP_2)
    Rule('reaction38', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_1, kAC5_GaiGTP_ATP_2)
    Rule('reaction39', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_ATP_GaiGTP_1, kAC5_ATP_GaiGTP_2)
    Rule('reaction40', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_GaiGTP_ATP_to_cAMP)
    Rule('reaction41', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction42', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGDP(), kAC5_GaiGTP_decay)
    Rule('reaction43', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGDP(), kAC5_GaiGTP_decay_2)

    # Interaction between AC5, Ca2+, and Gaolf
    Rule('reaction17', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_1 , kAC5_Ca_GaolfGTP_2)
    Rule('reaction31', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_ATP_1, kAC5_Ca_GaolfGTP_ATP_2)
    Rule('reaction32', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_CA_ATP_GaolfGTP_1, kAC5_CA_ATP_GaolfGTP_2)
    Rule('reaction33', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , kAC5_Ca_GaolfGTP_ATP_to_cAMP)
    Rule('reaction34', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse )
    Rule('reaction35', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None), kAC5_Ca_GaolfGTP_decay)
    Rule('reaction36', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGDP(), kAC5_Ca_GaolfGTP_ATP_decay)

    # Interaction between AC5, Gaolf, and Gai
    Rule('reaction44', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaolfGTP_GaiGTP_1, kAC5_GaolfGTP_GaiGTP_2)
    Rule('reaction45', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaolfGTP_GaiGTP_3, kAC5_GaolfGTP_GaiGTP_4)
    Rule('reaction46', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaolfGTP_GaiGTP_ATP_1, kAC5_GaolfGTP_GaiGTP_ATP_2 )
    Rule('reaction47', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_GaolfGTP_ATP_1, kAC5_GaiGTP_GaolfGTP_ATP_2)
    Rule('reaction48', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaolfGTP_ATP_GaiGTP_1, kAC5_GaolfGTP_ATP_GaiGTP_2)
    Rule('reaction49', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP)
    Rule('reaction50', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction51', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + GaolfGDP(), kAC5_GaolfGTP_GaolfGTP_decay_1)
    Rule('reaction52', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + GaiGDP(), kAC5_GaolfGTP_GaolfGTP_decay_2)
    Rule('reaction53', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) + GaolfGDP(), kAC5_GaolfGTP_GaolfGTP_decay_3)
    Rule('reaction54', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) + GaiGDP(), kAC5_GaolfGTP_GaolfGTP_decay_4)

    # Interaction between AC5, Ca2+, and Gai
    Rule('reaction55', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_1 , kAC5_Ca_GaiGTP_2)
    Rule('reaction56', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_1 , kAC5_Ca_GaiGTP_ATP_2 )
    Rule('reaction57', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_ATP_GaiGTP_1 , kAC5_Ca_ATP_GaiGTP_2)
    Rule('reaction58', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_Ca_GaiGTP_ATP_to_cAMP)
    Rule('reaction59', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction60', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_decay)
    Rule('reaction61', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_ATP_decay)

    # Interaction between AC5, Gaolf, Gai, and Ca2+
    Rule('reaction62', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaolfGTP_GaiGTP_1, kAC5_Ca_GaolfGTP_GaiGTP_2)
    Rule('reaction63', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=50) + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_GaolfGTP_1, kAC5_Ca_GaiGTP_GaolfGTP_2)
    Rule('reaction64', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaolfGTP_GaiGTP_ATP_1, kAC5_Ca_GaolfGTP_GaiGTP_ATP_2)
    Rule('reaction65', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaolfGTP_ATP_GaiGTP_1, kAC5_Ca_GaolfGTP_ATP_GaiGTP_2)
    Rule('reaction66', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=50) + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_GaolfGTP_1, kAC5_Ca_GaiGTP_ATP_GaolfGTP_2)
    Rule('reaction67', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP)
    Rule('reaction68', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction69', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + GaiGDP(), kAC5_Ca_GaolfGTP_GaiGTP_decay_1)
    Rule('reaction70', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + GaolfGDP(), kAC5_Ca_GaolfGTP_GaiGTP_decay_2)
    Rule('reaction71', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + GaiGDP(), kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_1)
    Rule('reaction72', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + GaolfGDP(), kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_2)

    '''Several Subtypes of Phosphodiesterases (PDE) degrade cAMP'''

    Rule('reaction73', PDE4(PDE4_s='i') + cAMP() | PDE4(PDE4_s='a'), kPDE4_cAMP_1, kPDE4_cAMP_2)
    Rule('reaction74', PDE4(PDE4_s='a') >> PDE4(PDE4_s='i') + AMP(), kPDE4_cAMP_to_AMP)
    Rule('reaction75', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() + cAMP() | PDE10(PDE10_c='Y', PDE10_s='i'), kPDE10_2cAMP_1, kPDE10_2cAMP_2)
    Rule('reaction76', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() | PDE10(PDE10_c='N', PDE10_s='a'),  kPDE10_cAMP_1, kPDE10_cAMP_2)
    Rule('reaction77', PDE10(PDE10_c='Y', PDE10_s='i') + cAMP() | PDE10(PDE10_c='Y', PDE10_s='a'), kPDE10_2cAMP_cAMP_1, kPDE10_2cAMP_cAMP_2)
    Rule('reaction78', PDE10(PDE10_c='N', PDE10_s='a') >> PDE10(PDE10_c='N', PDE10_s='i') + AMP(), kPDE10_cAMP_decay)
    Rule('reaction79', PDE10(PDE10_c='Y', PDE10_s='a') >> PDE10(PDE10_c='Y', PDE10_s='i') + AMP(), kPDE10_2cAMP_cAMP_decay)

    '''cAMP Activates PKA'''
    Rule('reaction80', PKA(PKA_s='cAMP0') + cAMP() + cAMP() | PKA(PKA_s='cAMP2'), kPKA_cAMP2_1, kPKA_cAMP2_2)
    Rule('reaction81', PKA(PKA_s='cAMP2') + cAMP() + cAMP() | PKA(PKA_s='cAMP4'), kPKA_cAMP4_1, kPKA_cAMP4_2)
    Rule('reaction82', PKA(PKA_s='cAMP4') | PKAc(PKAc_b1=None) + PKAreg(), kPKA_activation, kPKA_activation_reverse)

    '''GLUTAMATE AND CALCIUM-ACTIVATED ENZYMES
    - Calcium binds to CaM, which has four calcium-Binding sites_ 2 N sites (fast, high affinity), and 2 C sites (slow, low affinity).
    - Calcium binds to CaM in pairs, producing the active Ca4CaM in two consecutive steps. Ca4CaM binds to and activates several molecules in the striatum.
    - PP2B has a very high affinity for Ca4CaM, and a lower affinity for CaM than for CA4CaM, and the calcium dissossiation rate of PP2B-CaM and PP2B-Ca2CaM is slower than for PP2B-Ca4CaM.
    - CaMKII is a buffer for CaM:
        * reversibe bingind of CaM to CaMKII
        * slow autophosphorylation into CaMKIIp, in wich CaM is trapped
        * dephosphorylation of CaMKIIp by PP1 '''

    Rule('reaction83', CaM(CaM_b1=None, CaM_s='Ca0') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | CaM(CaM_b1=None, CaM_s='Ca2'), kCaM_Ca2_1, kCaM_Ca2_2)
    Rule('reaction84', CaM(CaM_b1=None, CaM_s='Ca2') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | CaM(CaM_b1=None, CaM_s='Ca4'), kCaM_Ca4_1, kCaM_Ca4_2)
    Rule('reaction85', CaM(CaM_b1= None, CaM_s='Ca0') + PP2B(PP2B_b1=None, PP2B_b2=None) |  CaM(CaM_b1=50, CaM_s='Ca0') % PP2B(PP2B_b1=50, PP2B_b2=None), kCaM_PP2B_1, kCaM_PP2B_2)
    Rule('reaction86', CaM(CaM_b1=None, CaM_s='Ca4') + PP2B(PP2B_b1=None, PP2B_b2=None) | CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=None), kCaMCa4_PP2B_1, kCaMCa4_PP2B_2)
    Rule('reaction87', CaM(CaM_b1=50, CaM_s='Ca0') % PP2B(PP2B_b1=50, PP2B_b2=None) + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | CaM(CaM_b1=50, CaM_s='Ca2') % PP2B(PP2B_b1=50, PP2B_b2=None), kCaMCa2_PP2B_1, kCaMCa2_PP2B_2)
    Rule('reaction88', CaM(CaM_b1=50, CaM_s='Ca2') % PP2B(PP2B_b1=50, PP2B_b2=None) + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=None), kCaMCa4_PP2B_3, kCaMCa4_PP2B_4)
    Rule('reaction89', CaM(CaM_b1=None, CaM_s='Ca2') + PP2B(PP2B_b1=None, PP2B_b2=None) | CaM(CaM_b1=50, CaM_s='Ca2') % PP2B(PP2B_b1=50, PP2B_b2=None), kCaMCa2_PP2B_3, kCaMCa2_PP2B_4 )
    Rule('reaction90', B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None) + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff'), kB72PP2A_Ca_1, kB72PP2A_Ca_2)

    '''PHOSPHORYLATON OF DARPP-32
    DARPP-32 activity is regulated by phosphorylation on two sites: Thr34 and Thr75.
    The catalytic subunit of PKA catalyzes phosphorilation of DARPP-32 at Thr34.
    Cdk5 catalyses phosphorylation of DARPP-32 at Thr75. experiments show that there is a basal level of phosphorylation at Thr75 due to Cdk5 activity.
    '''

    Rule('reaction91', PKAc(PKAc_b1=None) + DARPP32(DARPP32_b1=None) | PKAc(PKAc_b1=50) % DARPP32(DARPP32_b1=50), kPKAc_DARPP32_1, kPKAc_DARPP32_2)
    Rule('reaction92', PKAc(PKAc_b1=50) % DARPP32(DARPP32_b1=50) >> D32p34(D32p34_b1=None) + PKAc(PKAc_b1=None), kPKAc_DARPP32_decay)
    Rule('reaction93', PKAc(PKAc_b1=None) + B56PP2A(B56PP2A_b1=None, B56PP2A_p='N') | PKAc(PKAc_b1=50) % B56PP2A(B56PP2A_b1=50, B56PP2A_p='N'), kPKAc_B56PP2A_1, kPKAc_B56PP2A_2)
    Rule('reaction94', D32p34(D32p34_b1=None) + PP1(PP1_b1=None) | D32p34(D32p34_b1=50) % PP1(PP1_b1=50), kD23p34_PP1_1, kD23p34_PP1_2)
    Rule('reaction95', PKAc(PKAc_b1=50) % B56PP2A(B56PP2A_b1=50, B56PP2A_p='N') >> PKAc(PKAc_b1=None) + B56PP2A(B56PP2A_b1=None, B56PP2A_p='Y'), kB56PP2A_phosphorylation)
    Rule('reaction96', B56PP2A(B56PP2A_b1=None, B56PP2A_p='Y') >> B56PP2A(B56PP2A_b1=None, B56PP2A_p='N'), kB56PP2Ap_decay)
    Rule('reaction97', CDK5(CDK5_b1=None) + DARPP32(DARPP32_b1=None) | CDK5(CDK5_b1=50) % DARPP32(DARPP32_b1=50), kCDK5_DARPP32_1, kCDK5_DARPP32_2 )
    Rule('reaction98', CDK5(CDK5_b1=50) % DARPP32(DARPP32_b1=50) >> CDK5(CDK5_b1=None) + D32p75(D32p75_b1=None), kCDK5_DARPP32_decay)
    Rule('reaction99', PKAc(PKAc_b1=None) + D32p75(D32p75_b1=None) | PKAc(PKAc_b1=50) % D32p75(D32p75_b1=50), kPKAc_D32p75_1, kPKAc_D32p75_2)
    Rule('reaction100', B56PP2A(B56PP2A_b1=None, B56PP2A_p='Y') + D32p75(D32p75_b1=None) | B56PP2A(B56PP2A_b1=50, B56PP2A_p='Y') % D32p75(D32p75_b1=50), kB56PP2Ap_D32p75_1, kB56PP2Ap_D32p75_2)
    Rule('reaction101', B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None) + D32p75(D32p75_b1=None) | B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % D32p75(D32p75_b1=50), kB72PP2A_D32p75_1, kB72PP2A_D32p75_2)
    Rule('reaction102', B56PP2A(B56PP2A_b1=50, B56PP2A_p='Y') % D32p75(D32p75_b1=50) >> B56PP2A(B56PP2A_b1=None, B56PP2A_p='Y') + DARPP32(DARPP32_b1=None), kB56PP2Ap_D32p75_decay)
    Rule('reaction103', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % D32p75(D32p75_b1=50) >> B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None) + DARPP32(DARPP32_b1=None), kB72PP2A_D32p75_decay)
    Rule('reaction104', D32p75(D32p75_b1=None) + B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') | B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p75(D32p75_b1=20), kB72PP2ACa_D32p75_1, kB72PP2ACa_D32p75_2)
    Rule('reaction105', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p75(D32p75_b1=20) >> B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') + DARPP32(DARPP32_b1=None), kB72PP2ACa_D32p75_decay)
    Rule('reaction106', B56PP2A(B56PP2A_b1=50, B56PP2A_p='N') % D32p75(D32p75_b1=50) >> B56PP2A(B56PP2A_b1=None, B56PP2A_p='N') + DARPP32(DARPP32_b1=None), kB56PP2A_D32p75_decay)
    Rule('reaction107', D32p34(D32p34_b1=None) + CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=None) | D32p34(D32p34_b1=20) % CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=20), kPP2BCaM_D23p34_1, kPP2BCaM_D23p34_2)
    Rule('reaction108', D32p34(D32p34_b1=20) % CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=20) >> CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=None) + DARPP32(DARPP32_b1=None), kPP2BCaM_D23p34_decay)
    Rule('reaction109', D32p75(D32p75_b1=None) + B56PP2A(B56PP2A_b1=None, B56PP2A_p='N') | D32p75(D32p75_b1=50) % B56PP2A(B56PP2A_b1=50, B56PP2A_p='N'), kD32p75_B56PP2A_1, kD32p75_B56PP2A_2)
    Rule('reaction110', D32p34(D32p34_b1=None) + B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') | B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p34(D32p34_b1=20), kB72PP2ACa_D32p34_1, kB72PP2ACa_D32p34_2)
    Rule('reaction111', D32p34(D32p34_b1=None) + B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None) | B72PP2A(B72PP2A_b1=None, B72PP2A_b2=20) % D32p34(D32p34_b1=20), kB72PP2A_D32p34_1, kB72PP2A_D32p34_2  )
    Rule('reaction112', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p34(D32p34_b1=20) >> B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') + DARPP32(DARPP32_b1=None), kB72PP2ACa_D32p34_decay)
    Rule('reaction113', B72PP2A(B72PP2A_b1=None, B72PP2A_b2=20) % D32p34(D32p34_b1=20) >> B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None) + DARPP32(DARPP32_b1=None), kB72PP2A_D32p34_decay)

    ''''''
    Rule('reaction114', AKAR3(AKAR3_b1=None, AKAR3_p='N') + PKAc(PKAc_b1=None) | AKAR3(AKAR3_b1=50, AKAR3_p='N') % PKAc(PKAc_b1=50), kAKAR3_PKAc_1, kAKAR3_PKAc_2)
    Rule('reaction115', AKAR3(AKAR3_b1=50, AKAR3_p='N') % PKAc(PKAc_b1=50) >> AKAR3(AKAR3_b1=None, AKAR3_p='Y') + PKAc(PKAc_b1=None), kAKAR3_PKAc_decay)
    Rule('reaction116', AKAR3(AKAR3_b1=None, AKAR3_p='Y') + PP1(PP1_b1=None) | AKAR3(AKAR3_b1=50, AKAR3_p='Y') % PP1(PP1_b1=50), kAKAR3p_PP1_1, kAKAR3p_PP1_2)
    Rule('reaction117', AKAR3(AKAR3_b1=50, AKAR3_p='Y') % PP1(PP1_b1=50) >> AKAR3(AKAR3_b1=None, AKAR3_p='N') + PP1(PP1_b1=None), kAKAR3p_PP1_decay)

    # OBSERVABLES
    Observable('obs_D1R',D1R(D1R_b1=None, D1R_p='p0', D1R_s='i'))
    Observable('obs_D1RDA',D1R(D1R_b1=None, D1R_p='p0', D1R_s='a'))
    Observable('obs_Golf',Golf(Golf_b1=None))
    Observable('obs_Ach',Ach(Ach_b1=None))
    Observable('obs_M4R',M4R(M4R_b1=None, M4R_b2=None))
    Observable('obs_Gi',Gi(Gi_b1=None))
    Observable('obs_Gbgolf',Gbgolf())
    #Observable('obs_D1R_Golf',D1R(D1R_b1=None, D1R_b2=20) % Golf(Golf_b1=20))
    Observable('obs_D1R_DA_Golf',D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50))
    Observable('obs_M4R_Ach',M4R(M4R_b1=50, M4R_b2=None) % Ach(Ach_b1=50))
    Observable('obs_M4R_Gi',M4R(M4R_b1=None, M4R_b2=20) % Gi(Gi_b1=20))
    Observable('obs_M4R_Ach_Gi',M4R(M4R_b1=50, M4R_b2=20) % Ach(Ach_b1=50) % Gi(Gi_b1=20))
    Observable('obs_GaolfGTP',GaolfGTP(GaolfGTP_b1=None))
    Observable('obs_GaolfGDP',GaolfGDP())
    Observable('obs_GaiGTP',GaiGTP(GaiGTP_b1=None))
    Observable('obs_Gbgi',Gbgi())
    Observable('obs_GaiGDP', GaiGDP())
    Observable('obs_AC5',AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'))
    Observable('obs_AC5_GaolfGTP',AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50))
    Observable('obs_AC5_Ca',AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_Ca',Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'))
    Observable('obs_AC5_Ca_GaolfGTP', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50))
    Observable('obs_AC5_Ca_GaolfGTP_ATP', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50))
    Observable('obs_AC5_GaolfGTP_ATP',AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50))
    Observable('obs_ATP', ATP())
    Observable('obs_cAMP', cAMP())
    Observable('obs_AC5_ATP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'))
    Observable('obs_AC5_Ca_ATP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_AC5_GaiGTP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_GaiGTP_ATP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_GaolfGTP_GaiGTP', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_GaolfGTP_GaiGTP_ATP', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaiGTP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaiGTP_ATP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaolfGTP_GaiGTP', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaolfGTP_GaiGTP_ATP', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) % GaiGTP(GaiGTP_b1=10))
    Observable('obs_PDE4', PDE4(PDE4_s='i'))
    Observable('obs_PDE4_cAMP', PDE4(PDE4_s='a'))
    Observable('obs_AMP', AMP())
    Observable('obs_PDE10', PDE10(PDE10_c='N', PDE10_s='i'))
    Observable('obs_PDE10_2cAMP', PDE10(PDE10_c='Y', PDE10_s='i'))
    Observable('obs_PDE10_cAMP', PDE10(PDE10_c='N', PDE10_s='a'))
    Observable('obs_PDE10_2cAMP_cAMP', PDE10(PDE10_c='Y', PDE10_s='a'))
    Observable('obs_PKA', PKA(PKA_s='cAMP0'))
    Observable('obs_PKAcAMP2', PKA(PKA_s='cAMP2'))
    Observable('obs_PKAcAMP4', PKA(PKA_s='cAMP4'))
    Observable('obs_PKAc', PKAc(PKAc_b1=None))
    Observable('obs_PKAreg', PKAreg())
    Observable('obs_CaM', CaM(CaM_b1=None, CaM_s='Ca0'))
    Observable('obs_CaMCa2', CaM(CaM_b1=None, CaM_s='Ca2'))
    Observable('obs_CaMCa4', CaM(CaM_b1=None, CaM_s='Ca4'))
    Observable('obs_PP2B', PP2B(PP2B_b1=None, PP2B_b2=None))
    Observable('obs_CaM_PP2B', CaM(CaM_b1=50, CaM_s='Ca0') % PP2B(PP2B_b1=50, PP2B_b2=None))
    Observable('obs_CaMCa4_PP2B', CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=None))
    Observable('obs_CaMCa2_PP2B', CaM(CaM_b1=50, CaM_s='Ca2') % PP2B(PP2B_b1=50, PP2B_b2=None))
    Observable('obs_B72PP2A', B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None))
    Observable('obs_B72PP2A_Ca', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=None) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_DARPP32', DARPP32(DARPP32_b1=None))
    Observable('obs_PKAc_DARPP32', PKAc(PKAc_b1=50) % DARPP32(DARPP32_b1=50))
    Observable('obs_D23p34', D32p34(D32p34_b1=None))
    Observable('obs_B56PP2A', B56PP2A(B56PP2A_b1=None, B56PP2A_p='N'))
    Observable('obs_PKAc_B56PP2A', PKAc(PKAc_b1=None) % B56PP2A(B56PP2A_b1=None, B56PP2A_p='N'))
    Observable('obs_PP1', PP1(PP1_b1=None))
    Observable('obs_D32p34_PP1', D32p34(D32p34_b1=50) % PP1(PP1_b1=50))
    Observable('obs_B56PP2Ap', B56PP2A(B56PP2A_b1=None, B56PP2A_p='Y'))
    Observable('obs_CDK5', CDK5(CDK5_b1=None))
    Observable('obs_CDK5_DARPP32', CDK5(CDK5_b1=50) % DARPP32(DARPP32_b1=50))
    Observable('obs_D32p75', D32p75(D32p75_b1=None))
    Observable('obs_PKAc_D32p75', PKAc(PKAc_b1=50) % D32p75(D32p75_b1=50))
    Observable('obs_B56PP2Ap_D32p75', B56PP2A(B56PP2A_b1=50, B56PP2A_p='Y') % D32p75(D32p75_b1=50))
    Observable('obs_B72PP2A_D32p75', B56PP2A(B56PP2A_b1=50, B56PP2A_p='N') % D32p75(D32p75_b1=50))
    Observable('obs_B72PP2ACa_D32p75', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p75(D32p75_b1=20))
    Observable('obs_PP2BCaM_D23p34', D32p34(D32p34_b1=20) % CaM(CaM_b1=50, CaM_s='Ca4') % PP2B(PP2B_b1=50, PP2B_b2=20))
    Observable('obs_D32p75_B56PP2A', D32p75(D32p75_b1=50) % B56PP2A(B56PP2A_b1=50, B56PP2A_p='N'))
    Observable('obs_B72PP2ACa_D32p34', B72PP2A(B72PP2A_b1=50, B72PP2A_b2=20) % Ca(Ca_b1=50, Ca_l='cytos', Ca_s='buff') % D32p34(D32p34_b1=20))
    Observable('obs_B72PP2A_D32p34', B72PP2A(B72PP2A_b1=None, B72PP2A_b2=20) % D32p34(D32p34_b1=20))
    Observable('obs_AKAR3', AKAR3(AKAR3_b1=None, AKAR3_p='N'))
    Observable('obs_AKAR3_PKAc', AKAR3(AKAR3_b1=50, AKAR3_p='N') % PKAc(PKAc_b1=50))
    Observable('obs_AKAR3p', AKAR3(AKAR3_b1=None, AKAR3_p='Y'))
    Observable('obs_AKAR3p_PP1', AKAR3(AKAR3_b1=50, AKAR3_p='Y') % PP1(PP1_b1=50))

    return model



list_of_observables=['obs_D1R','obs_D1RDA','obs_Golf','obs_Ach','obs_M4R','obs_Gi','obs_Gbgolf','obs_D1R_DA_Golf','obs_M4R_Ach','obs_M4R_Gi','obs_M4R_Ach_Gi','obs_GaolfGTP','obs_GaolfGDP','obs_GaiGTP','obs_Gbgi','obs_GaiGDP','obs_AC5','obs_AC5_GaolfGTP','obs_AC5_Ca','obs_Ca','obs_AC5_Ca_GaolfGTP','obs_AC5_Ca_GaolfGTP_ATP','obs_AC5_GaolfGTP_ATP','obs_ATP','obs_cAMP','obs_AC5_ATP','obs_AC5_Ca_ATP','obs_AC5_GaiGTP','obs_AC5_GaiGTP_ATP','obs_AC5_GaolfGTP_GaiGTP','obs_AC5_GaolfGTP_GaiGTP_ATP','obs_AC5_Ca_GaiGTP','obs_AC5_Ca_GaiGTP_ATP','obs_AC5_Ca_GaolfGTP_GaiGTP','obs_AC5_Ca_GaolfGTP_GaiGTP_ATP','obs_PDE4','obs_PDE4_cAMP','obs_AMP','obs_PDE10','obs_PDE10_2cAMP','obs_PDE10_cAMP','obs_PDE10_2cAMP_cAMP','obs_PKA','obs_PKAcAMP2','obs_PKAcAMP4','obs_PKAc','obs_PKAreg','obs_CaM','obs_CaMCa2','obs_CaMCa4','obs_PP2B','obs_CaM_PP2B','obs_CaMCa4_PP2B','obs_CaMCa2_PP2B','obs_B72PP2A','obs_B72PP2A_Ca','obs_DARPP32','obs_PKAc_DARPP32','obs_D23p34','obs_B56PP2A','obs_PKAc_B56PP2A','obs_PP1','obs_D32p34_PP1','obs_B56PP2Ap','obs_CDK5','obs_CDK5_DARPP32','obs_D32p75','obs_PKAc_D32p75','obs_B56PP2Ap_D32p75','obs_B72PP2A_D32p75','obs_B72PP2ACa_D32p75','obs_PP2BCaM_D23p34','obs_D32p75_B56PP2A','obs_B72PP2ACa_D32p34','obs_B72PP2A_D32p34','obs_AKAR3','obs_AKAR3_PKAc','obs_AKAR3p','obs_AKAR3p_PP1']
