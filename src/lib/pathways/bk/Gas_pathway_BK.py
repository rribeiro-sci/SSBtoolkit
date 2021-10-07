##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                        G alpha s coupled receptor                            ##
##                                                                              ##
##                      PATHWAY based on D1 receptor                            ##
##                                                                              ##
##                    last modification 22 June 2019                            ##
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
    Monomer('Gbgolf')
    Monomer('GaolfGTP', ['GaolfGTP_b1'])
    Monomer('GaolfGDP')
    Monomer('AC5', ['AC5_b1', 'AC5_b2', 'AC5_b3', 'AC5_s'],  {'AC5_s': ['i', 'a']})
    Monomer('Ca', ['Ca_b1', 'Ca_l', 'Ca_s'], {'Ca_l':['cytos', 'ext'], 'Ca_s':['free' , 'buff']}) # l (location), s (state: free or buffed)
    Monomer('ATP')
    Monomer('cAMP')
    Monomer('AMP')
    Monomer('PDE4', ['PDE4_s'], {'PDE4_s':['i', 'a']})
    Monomer('PDE10', ['PDE10_c', 'PDE10_s'], {'PDE10_c':['N', 'Y'], 'PDE10_s':['i', 'a']})
    Monomer('PKA', ['PKA_s'], {'PKA_s':['cAMP0', 'cAMP2', 'cAMP4']})
    Monomer('PKAc', ['PKAc_b1'])
    Monomer('PKAreg')


    #INITIAL CONDITIONS
    Initial(D1R(D1R_b1=None, D1R_p='p0', D1R_s='a'), Parameter('D1RDA_0', LR))   # nM
    Initial(D1R(D1R_b1=None, D1R_p='p0', D1R_s='i'), Parameter('D1R_0', 2000.0-LR))   # nM
    Initial(Golf(Golf_b1=None), Parameter('Golf_0', 2000.0))             # nM
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', 700))  # nM
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', 60)) # nM
    Initial(ATP(), Parameter('ATP_0', 5000000))
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', 2000))
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', 700))
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', 1200))

    #PARAMETERS
    Parameter('kD1R_Golf_1', 0.003)     # 1/(nM*s) |Association constant of the complex D1R_Golf
    Parameter('kD1R_Golf_2', 5.0)       # 1/s      |Dissociation constant of the complex D1R_Golf
    Parameter('kD1R_DA_Golf_1', 0.003)  # 1/(nM*s) |Association constant of the complex D1R_DA_Golf
    Parameter('kD1R_DA_Golf_2', 5.0)    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
    Parameter('kD1R_Golf_DA_1', 0.005)  # 1/(nM*s) |Association constant of the complex D1R_Golf_DA
    Parameter('kD1R_Golf_DA_2', 5.0)    # 1/s      |Dissociation constant of the complex D1R_Golf_DA
    Parameter('kM4R_Ach_1',0.01)        # 1/(nM*s) |Association constant of the complex M4R_Ach
    Parameter('kM4R_Ach_2',90)          # 1/s      |Dissociation constant of the complex M4R_Ach

    Parameter('kD1R_DA_Golf_decay', 15) # 1/s      |Rate of formation of M4R_DA + Golf
    Parameter('kGaolfGTP_decay', 30)    # 1/s      |Rate of convertion of GaolfGTP into GaolfGDP
    Parameter('kGolf_formation', 100)   # 1/s      |Rate of formation of Golf

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

    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_1', 0.2)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_2', 30)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_3', 0.2)    # 1/s      |
    Parameter('kAC5_GaolfGTP_GaolfGTP_decay_4', 30)    # 1/s      |

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

    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The D1R can bind either the inactive G protein first, and then
    dopamine, or dopamine first and then and then the inactivate G protein.'''

    Rule('reaction3', D1R(D1R_b1=None, D1R_p='p0', D1R_s='a') + Golf(Golf_b1=None) | D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50), kD1R_DA_Golf_1, kD1R_DA_Golf_2)
    Rule('reaction9', D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50) >> D1R(D1R_b1=None, D1R_p='p0', D1R_s='a') + Gbgolf() + GaolfGTP(GaolfGTP_b1=None) , kD1R_DA_Golf_decay)
    Rule('reaction10', GaolfGTP(GaolfGTP_b1=None) >> GaolfGDP(), kGaolfGTP_decay)
    Rule('reaction11', GaolfGDP() + Gbgolf() >> Golf(Golf_b1=None), kGolf_formation)

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

    # Interaction between AC5, Ca2+, and Gaolf
    Rule('reaction17', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_1 , kAC5_Ca_GaolfGTP_2)
    Rule('reaction31', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_ATP_1, kAC5_Ca_GaolfGTP_ATP_2)
    Rule('reaction32', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_CA_ATP_GaolfGTP_1, kAC5_CA_ATP_GaolfGTP_2)
    Rule('reaction33', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , kAC5_Ca_GaolfGTP_ATP_to_cAMP)
    Rule('reaction34', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse )
    Rule('reaction35', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None), kAC5_Ca_GaolfGTP_decay)
    Rule('reaction36', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGDP(), kAC5_Ca_GaolfGTP_ATP_decay)

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
    Observable('obs_D1R',D1R(D1R_b1=None, D1R_p='p0', D1R_s='i'))
    Observable('obs_D1RDA',D1R(D1R_b1=None, D1R_p='p0', D1R_s='a'))
    #Observable('obs_D1R_Golf',D1R(D1R_b1=None, D1R_b2=20) % Golf(Golf_b1=20))
    Observable('obs_D1R_DA_Golf',D1R(D1R_b1=50, D1R_p='p0', D1R_s='a') % Golf(Golf_b1=50))
    Observable('obs_cAMP', cAMP())
    Observable('obs_AC5', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'))
    return model

list_of_observables=['obs_D1R','obs_D1RDA','obs_cAMP', 'obs_AC5']
