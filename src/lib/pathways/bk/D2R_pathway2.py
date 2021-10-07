##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                              D2 receptor                                     ##
##                                                                              ##
##                               PATHWAY 1                                      ##
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
    Monomer('D2R', ['D2R_b1', 'D2R_p', 'D2R_s'], {'D2R_p':['p0','p1'], 'D2R_s':['i','a']})    # Dopamine receptor 1 with two binding sites: for DA and G-proteins
    Monomer('Golf', ['Golf_b1'])            # Golf protein (alpha and beta/gamma complex)
    Monomer('Gi', ['Gi_b1'])
    Monomer('Gbgolf')
    Monomer('GaolfGTP', ['GaolfGTP_b1'])
    Monomer('GaolfGDP')
    Monomer('GaiGTP', ['GaiGTP_b1'])
    Monomer('Gbgi')
    Monomer('Gbi')
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

    #INITIAL CONDITIONS
    Initial(D2R(D2R_b1=None, D2R_p='p0', D2R_s='a'), Parameter('D2RDA_0', LR))   # nM
    Initial(D2R(D2R_b1=None, D2R_p='p0', D2R_s='i'), Parameter('D2R_0', 2000.0-LR))   # nM    Initial(Golf(Golf_b1=None), Parameter('Golf_0', 1200.0))             # nM
    Initial(Gi(Gi_b1=None), Parameter('Gi_0', 3000.0))                   # nM
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', 700))  # nM
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', 60)) # nM
    Initial(ATP(), Parameter('ATP_0', 5000000))
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', 2000))
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', 1000))
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', 1200))
    Initial(CaM(CaM_b1=None, CaM_s='Ca0'), Parameter('CaM_0', 10000))
    Initial(PP2B(PP2B_b1=None, PP2B_b2=None), Parameter('PP2B_0', 4000))
    Initial(B72PP2A(B72PP2A_b1=None, B72PP2A_b2=None), Parameter('B72PP2A_0', 700))
    Initial(DARPP32(DARPP32_b1=None), Parameter('DARPP32_0', 50000))
    Initial(B56PP2A(B56PP2A_b1=None, B56PP2A_p='N'), Parameter('B56PP2A_0', 2000))
    Initial(PP1(PP1_b1=None), Parameter('PP1_0', 1000))
    Initial(CDK5(CDK5_b1=None), Parameter('CDK5_0', 1800))

    Initial(cAMP(), Parameter('cAMP_0', 65000))

    #PARAMETERS

    #reaction2
    Parameter('kD2R_Gi_1', 0.055)     # 1/(nM*s) |Association constant of the complex D2R_Golf
    Parameter('kD2R_Gi_2', 200)       # 1/s      |Dissociation constant of the complex D2R_Golf

    #reaction9
    Parameter('kD2R_DA_Gi_decay', 60) # 1/s      |Rate of formation of M4R_DA + Golf
    #reaction10
    Parameter('kGaiGTP_decay', 30)    # 1/s      |Rate of convertion of GaolfGTP into GaolfGDP
    #reaction11
    Parameter('kGi_formation', 100)   # 1/s      |Rate of formation of Golf
    #reaction16
    Parameter('kAC5_Ca_1', 0.001)       # 1/(nM*s) |Association constant of the complex AC5_Ca
    Parameter('kAC5_Ca_2', 0.9)         # 1/s      |Dissociation constant of the complex AC5_Ca
    #reaction21
    Parameter('kAC5_ATP_1', 0.0001)     # 1/(nM*s) |Association constant of the complex AC5_ATP
    Parameter('kAC5_ATP_2', 1)          # 1/s      |Dissociation constant of the complex AC5_ATP
    #reaction22
    Parameter('kAC5_basal', 1)          # 1/s      |Rate constant of basal formation of cAMP
    #reaction23
    Parameter('kAC5_reverse_basal', 0.0004) # 1/s  |Reverse Rate constant of basal formation of cAMP
    #reaction24
    Parameter('kAC5_Ca_ATP_1', 7.50E-5) # 1/(nM*s) |Association constant of the complex AC5_Ca_ATP
    Parameter('kAC5_Ca_ATP_2', 1)       # 1/s      |Dissociation constant of the complex AC5_Ca_ATP
    #reaction25
    Parameter('kAC5_Ca_ATP_to_cAMP', 0.5 )  # 1/s      |Convertion rate of ATP into cAMP by AC5_Ca
    #reaction26
    Parameter('kAC5_Ca_ATP_to_cAMP_reverse', 0.00015) # 1/s      |Reverse Convertion rate of ATP into cAMP by AC5_Ca
    #reaction27
    Parameter('kAC5_GaolfGTP_ATP_3', 0.2)   # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_ATP_4', 0.1)   # 1/s      |
    #reaction28
    Parameter('kAC5_GaolfGTP_decay', 0.25)   # 1/s      |
    #reaction29
    Parameter('kAC5_GaolfGTP_ATP_decay', 0.25)# 1/s     |
    #reaction30
    Parameter('kAC5_ATP_Ca_1', 0.001)       # 1/(nM*s) |
    Parameter('kAC5_ATP_Ca_2', 0.9)         # 1/s      |
    #reaction31
    Parameter('kAC5_Ca_GaolfGTP_ATP_1', 0.00055)   # 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_ATP_2', 1)   # 1/s      |
    #reaction32
    Parameter('kAC5_CA_ATP_GaolfGTP_1', 0.2) # 1/(nM*s) |
    Parameter('kAC5_CA_ATP_GaolfGTP_2', 0.1) # 1/s      |
    #reaction33
    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP', 10)# 1/s      |
    #reaction34
    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse', 0.022) # 1/s      |
    #reaction35
    Parameter('kAC5_Ca_GaolfGTP_decay', 0.25) # 1/s      |
    #reaction36
    Parameter('kAC5_Ca_GaolfGTP_ATP_decay', 0.25) # 1/s      |
    #reaction37
    Parameter('kAC5_GaiGTP_1', 50)          # 1/(nM*s) |
    Parameter('kAC5_GaiGTP_2', 5)           # 1/s      |
    #reaction38
    Parameter('kAC5_GaiGTP_ATP_1', 6.25E-5) # 1/(nM*s) |
    Parameter('kAC5_GaiGTP_ATP_2', 1)       # 1/s      |
    #reaction39
    Parameter('kAC5_ATP_GaiGTP_1', 50)      # 1/(nM*s) |
    Parameter('kAC5_ATP_GaiGTP_2', 5)       # 1/s      |
    #reaction40
    Parameter('kAC5_GaiGTP_ATP_to_cAMP', 0.25) # 1/s      |
    #reaction41
    Parameter('kAC5_GaiGTP_ATP_to_cAMP_reverse', 0.00105) # 1/s      |
    #reaction42
    Parameter('kAC5_GaiGTP_decay', 30)       # 1/s      |
    #reaction43
    Parameter('kAC5_GaiGTP_decay_2', 30)     # 1/s      |
    #reaction44
    Parameter('kAC5_GaolfGTP_GaiGTP_1', 0.01) # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_2', 0.01) # 1/s      |
    #reaction45
    Parameter('kAC5_GaolfGTP_GaiGTP_3', 0.002) # 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_4', 0.01) # 1/s      |
    #reaction46
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_1', 0.0003)# 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_2', 1)     # 1/s      |
    #reaction47
    Parameter('kAC5_GaiGTP_GaolfGTP_ATP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_GaiGTP_GaolfGTP_ATP_2', 0.01)# 1/s      |
    #reaction48
    Parameter('kAC5_GaolfGTP_ATP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_GaolfGTP_ATP_GaiGTP_2', 0.01)# 1/s      |
    #reaction49
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP', 5)# 1/s      |
    #reaction50
    Parameter('kAC5_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse', 0.006) # 1/s      |
    #reaction51
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_1', 0.25)    # 1/s      |
    #reaction52
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_2', 30)    # 1/s      |
    #reaction53
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_3', 30)    # 1/s      |
    #reaction54
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_4', 30)    # 1/s      |
    #reaction55
    Parameter('kAC5_Ca_GaiGTP_1', 50)   # 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_2', 5)    # 1/s      |
    #reaction56
    Parameter('kAC5_Ca_GaiGTP_ATP_1', 5.63E-5)  # 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_ATP_2', 1)    # 1/s      |
    #reaction57
    Parameter('kAC5_Ca_ATP_GaiGTP_1', 50)   # 1/(nM*s) |
    Parameter('kAC5_Ca_ATP_GaiGTP_2', 5)    # 1/s      |
    #reaction58
    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP', 0.125) # 1/s      |
    #reaction59
    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse', 2.81E-5)# 1/s      |
    #reaction60
    Parameter('kAC5_Ca_GaiGTP_decay', 30)# 1/s      |
    #reaction61
    Parameter('kAC5_Ca_GaiGTP_ATP_decay', 30)# 1/s      |
    #reaction62
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_2', 0.01)# 1/s      |
    #reaction63
    Parameter('kAC5_Ca_GaiGTP_GaolfGTP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_GaolfGTP_2',0.01)# 1/s      |
    #reaction64
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_1', 0.000175)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_2', 1)# 1/s      |
    #reaction65
    Parameter('kAC5_Ca_GaolfGTP_ATP_GaiGTP_1', 0.01)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaolfGTP_ATP_GaiGTP_2', 0.01)# 1/s      |
    #reaction66
    Parameter('kAC5_Ca_GaiGTP_ATP_GaolfGTP_1', 0.002)# 1/(nM*s) |
    Parameter('kAC5_Ca_GaiGTP_ATP_GaolfGTP_2', 0.01)# 1/s      |
    #reaction67
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP', 2.5)# 1/s      |
    #reaction68
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_to_cAMP_reverse', 0.00175)# 1/s      |
    #reaction69
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_decay_1', 30)# 1/s      |
    #reaction70
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_decay_2', 0.25)# 1/s      |
    #reaction71
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_1', 30)# 1/s      |
    #reaction72
    Parameter('kAC5_Ca_GaolfGTP_GaiGTP_ATP_decay_2', 0.25)# 1/s      |
    #reaction73
    Parameter('kPDE4_cAMP_1', 0.01)# 1/(nM*s) |
    Parameter('kPDE4_cAMP_2', 1)# 1/s      |
    #reaction74
    Parameter('kPDE10_2cAMP_1', 1.0E-6)# 1/(nM*s) |
    Parameter('kPDE10_2cAMP_2', 9)# 1/s      |
    #reaction75
    Parameter('kPDE4_cAMP_to_AMP', 2)
    #reaction76
    Parameter('kPDE10_cAMP_1', 0.1)# 1/(nM*s) |
    Parameter('kPDE10_cAMP_2', 2)# 1/s      |
    #reaction77
    Parameter('kPDE10_2cAMP_cAMP_1', 0.13)# 1/(nM*s) |
    Parameter('kPDE10_2cAMP_cAMP_2', 2)# 1/s      |
    #reaction78
    Parameter('kPDE10_cAMP_decay', 3)# 1/s      |
    #reaction79
    Parameter('kPDE10_2cAMP_cAMP_decay', 10)# 1/s      |
    #reaction80
    Parameter('kPKA_cAMP2_1', 0.00026)# 1/(nM*s) |
    Parameter('kPKA_cAMP2_2', 1)# 1/s |
    #reaction81
    Parameter('kPKA_cAMP4_1', 0.000346)# 1/(nM*s) |
    Parameter('kPKA_cAMP4_2', 1)# 1/s |
    #reaction82
    Parameter('kPKA_activation', 10)# 1/(nM*s) |
    Parameter('kPKA_activation_reverse', 0.025)# 1/s |
    #reaction83
    Parameter('kCaM_Ca2_1', 0.006)# 1/(nM*s) |
    Parameter('kCaM_Ca2_2', 9.1)# 1/s |
    #reaction84
    Parameter('kCaM_Ca4_1', 0.1)# 1/(nM*s) |
    Parameter('kCaM_Ca4_2', 1000)# 1/s |
    #reaction85
    Parameter('kCaM_PP2B_1', 1)# 1/(nM*s) |
    Parameter('kCaM_PP2B_2', 30)# 1/s |
    #reaction86
    Parameter('kCaMCa4_PP2B_1', 1)# 1/(nM*s) |
    Parameter('kCaMCa4_PP2B_2', 3)# 1/s |
    #reaction87
    Parameter('kCaMCa2_PP2B_1', 0.006)# 1/(nM*s) |
    Parameter('kCaMCa2_PP2B_2', 0.91)# 1/s |
    #reaction88
    Parameter('kCaMCa4_PP2B_3', 0.1)# 1/(nM*s) |
    Parameter('kCaMCa4_PP2B_4', 10)# 1/s |
    #reaction89
    Parameter('kCaMCa2_PP2B_3', 1)# 1/(nM*s) |
    Parameter('kCaMCa2_PP2B_4', 3)# 1/s |
    #reaction90
    Parameter('kB72PP2A_Ca_1', 0.0001)# 1/(nM*s) |
    Parameter('kB72PP2A_Ca_2', 0.1)# 1/s |
    #reaction91
    Parameter('kPKAc_DARPP32_1', 0.001)# 1/(nM*s) |
    Parameter('kPKAc_DARPP32_2', 8)# 1/s |
    #reaction92
    Parameter('kPKAc_DARPP32_decay', 10)# 1/s |
    #reaction93
    Parameter('kPKAc_B56PP2A_1', 0.005)# 1/(nM*s) |
    Parameter('kPKAc_B56PP2A_2', 0.3)# 1/s |
    #reaction94
    Parameter('kD23p34_PP1_1', 0.02)# 1/(nM*s) |
    Parameter('kD23p34_PP1_2', 0.04)# 1/s |
    #reaction95
    Parameter('kB56PP2A_phosphorylation', 0.1) # 1/s |
    #reaction96
    Parameter('kB56PP2Ap_decay', 0.006) # 1/s |
    #reaction97
    Parameter('kCDK5_DARPP32_1', 0.0009)# 1/(nM*s) |
    Parameter('kCDK5_DARPP32_2', 2)# 1/s |
    #reaction98
    Parameter('kCDK5_DARPP32_decay', 3)# 1/s |
    #reaction99
    Parameter('kPKAc_D32p75_1', 0.00037)# 1/(nM*s) |
    Parameter('kPKAc_D32p75_2', 1)# 1/s |
    #reaction100
    Parameter('kB56PP2Ap_D32p75_1', 0.001)# 1/(nM*s) |
    Parameter('kB56PP2Ap_D32p75_2', 10)# 1/s |
    #reaction101
    Parameter('kB72PP2A_D32p75_1', 0.0008) # 1/(nM*s) |
    Parameter('kB72PP2A_D32p75_2', 6.4)# 1/s |
    #reaction102
    Parameter('kB56PP2Ap_D32p75_decay', 8)# 1/s |
    #reaction103
    Parameter('kB72PP2A_D32p75_decay', 1)# 1/s |
    #reaction104
    Parameter('kB72PP2ACa_D32p75_1', 0.001)# 1/(nM*s) |
    Parameter('kB72PP2ACa_D32p75_2', 10 )# 1/s |
    #reaction105
    Parameter('kB72PP2ACa_D32p75_decay', 5)# 1/s |
    #reaction106
    Parameter('kB56PP2A_D32p75_decay', 2.3)# 1/s |
    #reaction107
    Parameter('kPP2BCaM_D23p34_1', 0.005)# 1/(nM*s) |
    Parameter('kPP2BCaM_D23p34_2', 1)# 1/s |
    #reaction108
    Parameter('kPP2BCaM_D23p34_decay', 7)# 1/s |
    #reaction109
    Parameter('kD32p75_B56PP2A_1', 0.0008)# 1/(nM*s) |
    Parameter('kD32p75_B56PP2A_2', 6.4)# 1/s |
    #reaction110
    Parameter('kB72PP2ACa_D32p34_1', 0.1)# 1/(nM*s) |
    Parameter('kB72PP2ACa_D32p34_2', 1) # 1/s |
    #reaction111
    Parameter('kB72PP2A_D32p34_1', 0.1)# 1/(nM*s) |
    Parameter('kB72PP2A_D32p34_2', 1) # 1/s |
    #reaction112
    Parameter('kB72PP2ACa_D32p34_decay', 0.2)# 1/s |
    #reaction113
    Parameter('kB72PP2A_D32p34_decay', 0.2)# 1/s


    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The D2R can bind either the inactive G protein first, and then
    dopamine, or dopamine first and then and then the inactivate G protein.'''

    Rule('reaction3', D2R(D2R_b1=None, D2R_p='p0', D2R_s='a') + Gi(Gi_b1=None) | D2R(D2R_b1=20, D2R_p='p0', D2R_s='a') % Gi(Gi_b1=20), kD2R_Gi_1, kD2R_Gi_2)
    Rule('reaction9', D2R(D2R_b1=20, D2R_p='p0', D2R_s='a') % Gi(Gi_b1=20) >> D2R(D2R_b1=None, D2R_p='p0', D2R_s='a') + Gbi() + GaiGTP(GaiGTP_b1=None) , kD2R_DA_Gi_decay)
    Rule('reaction10', GaiGTP(GaiGTP_b1=None) >> GaiGDP(), kGaiGTP_decay)
    Rule('reaction11', GaiGDP() + Gbi() >> Gi(Gi_b1=None), kGi_formation)



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


    # Interaction between AC5 and Gai
    Rule('reaction37', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) , kAC5_GaiGTP_1, kAC5_GaiGTP_2)
    Rule('reaction38', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_1, kAC5_GaiGTP_ATP_2)
    Rule('reaction39', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_ATP_GaiGTP_1, kAC5_ATP_GaiGTP_2)
    Rule('reaction40', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_GaiGTP_ATP_to_cAMP)
    Rule('reaction41', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction42', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGDP(), kAC5_GaiGTP_decay)
    Rule('reaction43', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGDP(), kAC5_GaiGTP_decay_2)


    # Interaction between AC5, Ca2+, and Gai
    Rule('reaction55', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_1 , kAC5_Ca_GaiGTP_2)
    Rule('reaction56', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_1 , kAC5_Ca_GaiGTP_ATP_2 )
    Rule('reaction57', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_ATP_GaiGTP_1 , kAC5_Ca_ATP_GaiGTP_2)
    Rule('reaction58', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_Ca_GaiGTP_ATP_to_cAMP)
    Rule('reaction59', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse)
    Rule('reaction60', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_decay)
    Rule('reaction61', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_ATP_decay)


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


    # OBSERVABLES
    Observable('obs_AC5_active', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'))

    Observable('obs_D2R',D2R(D2R_b1=None, D2R_p='p0', D2R_s='a'))
    Observable('obs_Gi',Gi(Gi_b1=None))
    Observable('obs_D2R_Gi',D2R(D2R_b1=20, D2R_p='p0', D2R_s='a') % Gi(Gi_b1=20))
    Observable('obs_GaiGTP',GaiGTP(GaiGTP_b1=None))
    Observable('obs_Gbgi',Gbgi())
    Observable('obs_GaiGDP', GaiGDP())
    Observable('obs_AC5',AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'))
    Observable('obs_AC5_Ca',AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_Ca',Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'))
    Observable('obs_ATP', ATP())
    Observable('obs_cAMP', cAMP()+PDE4(PDE4_s='a')+PDE10(PDE10_c='Y', PDE10_s='i')+PDE10(PDE10_c='N', PDE10_s='a')+PDE10(PDE10_c='Y', PDE10_s='a')+PKA(PKA_s='cAMP0')+PKA(PKA_s='cAMP2')+PKA(PKA_s='cAMP4'))
    Observable('obs_AC5_ATP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'))
    Observable('obs_AC5_Ca_ATP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_AC5_GaiGTP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_GaiGTP_ATP', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaiGTP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10))
    Observable('obs_AC5_Ca_GaiGTP_ATP', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10))
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
    return model



list_of_observables=['obs_AC5_active','obs_D2R','obs_Gi','obs_D2R_Gi','obs_GaiGTP','obs_Gbgi',
'obs_GaiGDP','obs_AC5','obs_AC5_Ca','obs_Ca','obs_ATP','obs_cAMP','obs_AC5_ATP','obs_AC5_Ca_ATP',
'obs_AC5_GaiGTP', 'obs_AC5_GaiGTP_ATP','obs_AC5_Ca_GaiGTP','obs_AC5_Ca_GaiGTP_ATP','obs_PDE4',
'obs_PDE4_cAMP','obs_AMP','obs_PDE10','obs_PDE10_2cAMP','obs_PDE10_cAMP','obs_PDE10_2cAMP_cAMP',
'obs_PKA','obs_PKAcAMP2','obs_PKAcAMP4','obs_PKAc','obs_PKAreg']


labels=[{'label': 'Receptor*Ligand', 'value':'D2R'},
        {'label': 'ATP', 'value':'ATP'},
        {'label': 'cAMP', 'value':'cAMP'},
        {'label': 'AC5', 'value':'AC5_active'},
        {'label': 'PDE4 (activated)', 'value':'PDE4_cAMP'},
        {'label': 'PKA (activated)', 'value':'PKAc'},
        ]
