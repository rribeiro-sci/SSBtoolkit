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


def network(Rtotal, LR):

    #Start a model
    Model()

    #MONOMERS
    Monomer('R', ['R_b1', 'R_p', 'R_s'], {'R_p':['p0','p1'], 'R_s':['i','a']})    # Dopamine receptor 1 with two binding sites: for ligand and G-proteins
    Monomer('Gi', ['Gi_b1'])
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


    #INITIAL CONDITIONS μM
    Initial(R(R_b1=None, R_p='p0', R_s='a'), Parameter('RL_0', LR))   # μM
    Initial(R(R_b1=None, R_p='p0', R_s='i'), Parameter('R_0', Rtotal-LR))   # μM
    Initial(Gi(Gi_b1=None), Parameter('Gi_0', 2))                   # μM
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', 0.7))  # μM
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', 0.06)) # μM
    Initial(ATP(), Parameter('ATP_0', 5000))# μM
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', 2))# μM
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', 0.7))# μM
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', 1.2))# μM

    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The R can bind either the inactive G protein first, and then
    dopamine, or dopamine first and then and then the inactivate G protein.'''
    Parameter('kRL_Gi_1',1.2*1E3)      # 1/(μM*s) |Association constant of the complex M4R_Ach_Gi
    Parameter('kRL_Gi_2',90)       # 1/s      |Dissociation constant of the complex M4R_Ach_Gi
    Rule('reaction7', R(R_b1=None, R_p='p0', R_s='a') + Gi(Gi_b1=None) | R(R_b1=50, R_p='p0', R_s='a') % Gi(Gi_b1=50), kRL_Gi_1, kRL_Gi_2)

    Parameter('kRL_Gi_decay', 60)   # 1/s      |Rate of formation of M4R_DA + Gi
    Rule('reaction8', R(R_b1=50, R_p='p0', R_s='a') % Gi(Gi_b1=50) >> R(R_b1=None, R_p='p0', R_s='a') + Gbgi() + GaiGTP(GaiGTP_b1=None) , kRL_Gi_decay)

    Parameter('kGaiGTP_decay', 30)      # 1/s      |Rate of formation of GaiGDP
    Rule('reaction13', GaiGTP(GaiGTP_b1=None) >> GaiGDP(), kGaiGTP_decay)

    Parameter('kGi_formation', 100)     # 1/s      |Rate of formation of Gi
    Rule('reaction14', Gbgi() + GaiGDP() >> Gi(Gi_b1=None), kGi_formation)

    '''CYCLIC AMP FORMATION AND PKA ACTIVATION:
    - active GaolfGTP binds to and activate cyclase type V (AC5), which produces cyclic AMP (cAMP)
    - Ca2+ causes a 50% reduction in the activity of AC5 before or after binding to Gaolf. In the model
        calcium is allowed to bind to the inactive AC5, which then can be activated by GaolfGTP.
    - A step of regeneration of ATP is included that allows for a steady state concentration of ATP in
    the absence of stimulation. It is not assumed to be a rate-limiting substrate.
    '''

    # AC5 basal activity
    Parameter('kAC5_ATP_1', 0.0001*1E3)     # 1/(μM*s) |Association constant of the complex AC5_ATP
    Parameter('kAC5_ATP_2', 1)          # 1/s      |Dissociation constant of the complex AC5_ATP
    Rule('reaction21', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') , kAC5_ATP_1, kAC5_ATP_2)

    Parameter('kAC5_basal', 1)          # 1/s      |Rate constant of basal formation of cAMP
    Rule('reaction22', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP(), kAC5_basal)

    Parameter('kAC5_reverse_basal', 0.0004) # 1/s  |Reverse Rate constant of basal formation of cAMP
    Rule('reaction23', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'), kAC5_reverse_basal )

    # Interaction between AC5 and Ca2+
    Parameter('kAC5_Ca_1', 0.001*1E3)       # 1/(μM*s) |Association constant of the complex AC5_Ca
    Parameter('kAC5_Ca_2', 0.9)         # 1/s      |Dissociation constant of the complex AC5_Ca
    Rule('reaction16', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_1, kAC5_Ca_2)

    Parameter('kAC5_Ca_ATP_1', 7.50E-5*1E3) # 1/(μM*s) |Association constant of the complex AC5_Ca_ATP
    Parameter('kAC5_Ca_ATP_2', 1)       # 1/s      |Dissociation constant of the complex AC5_Ca_ATP
    Rule('reaction24', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_1, kAC5_Ca_ATP_2)

    Parameter('kAC5_Ca_ATP_to_cAMP', 0.5 )  # 1/s      |Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction25', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') >> cAMP() + AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP)

    Parameter('kAC5_Ca_ATP_to_cAMP_reverse', 0.00015) # 1/s      |Reverse Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction26', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP_reverse)

    Parameter('kAC5_ATP_Ca_1', 0.001*1E3)       # 1/(μM*s) |
    Parameter('kAC5_ATP_Ca_2', 0.9)         # 1/s      |
    Rule('reaction30', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') , kAC5_ATP_Ca_1, kAC5_ATP_Ca_2)

    # Interaction between AC5 and Gai

    Parameter('kAC5_GaiGTP_1', 50*1E3)          # 1/(μM*s) |
    Parameter('kAC5_GaiGTP_2', 5)           # 1/s
    Rule('reaction37', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) , kAC5_GaiGTP_1, kAC5_GaiGTP_2)

    Parameter('kAC5_GaiGTP_ATP_1', 6.25E-5*1E3) # 1/(μM*s) |
    Parameter('kAC5_GaiGTP_ATP_2', 1)       # 1/s      |
    Rule('reaction38', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_1, kAC5_GaiGTP_ATP_2)

    Parameter('kAC5_ATP_GaiGTP_1', 50*1E3)      # 1/(μM*s) |
    Parameter('kAC5_ATP_GaiGTP_2', 5)       # 1/s
    Rule('reaction39', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_ATP_GaiGTP_1, kAC5_ATP_GaiGTP_2)

    Parameter('kAC5_GaiGTP_ATP_to_cAMP', 0.25) # 1/s      |
    Rule('reaction40', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_GaiGTP_ATP_to_cAMP)

    Parameter('kAC5_GaiGTP_ATP_to_cAMP_reverse', 0.00105) # 1/s      |
    Rule('reaction41', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), kAC5_GaiGTP_ATP_to_cAMP_reverse)

    Parameter('kAC5_GaiGTP_decay', 30)       # 1/s      |
    Rule('reaction42', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGDP(), kAC5_GaiGTP_decay)

    Parameter('kAC5_GaiGTP_decay_2', 30)     # 1/s      |
    Rule('reaction43', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGDP(), kAC5_GaiGTP_decay_2)

    # Interaction between AC5, Ca2+, and Gai

    Parameter('kAC5_Ca_GaiGTP_1', 50*1E3)   # 1/(μM*s) |
    Parameter('kAC5_Ca_GaiGTP_2', 5)    # 1/s      |
    Rule('reaction55', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_1 , kAC5_Ca_GaiGTP_2)

    Parameter('kAC5_Ca_GaiGTP_ATP_1', 5.63E-5*1E3)  # 1/(μM*s) |
    Parameter('kAC5_Ca_GaiGTP_ATP_2', 1)    # 1/s      |
    Rule('reaction56', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_1 , kAC5_Ca_GaiGTP_ATP_2 )

    Parameter('kAC5_Ca_ATP_GaiGTP_1', 50*1E3)   # 1/(μM*s)
    Parameter('kAC5_Ca_ATP_GaiGTP_2', 5)    # 1/s      |
    Rule('reaction57', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_ATP_GaiGTP_1 , kAC5_Ca_ATP_GaiGTP_2)

    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP', 0.125) # 1/s      |
    Rule('reaction58', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP(), kAC5_Ca_GaiGTP_ATP_to_cAMP)

    Parameter('kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse', 2.81E-5)# 1/s      |
    Rule('reaction59', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), kAC5_Ca_GaiGTP_ATP_to_cAMP_reverse)

    Parameter('kAC5_Ca_GaiGTP_decay', 30)# 1/s      |
    Rule('reaction60', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_decay)

    Parameter('kAC5_Ca_GaiGTP_ATP_decay', 30)# 1/s      |
    Rule('reaction61', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), kAC5_Ca_GaiGTP_ATP_decay)

    '''Several Subtypes of Phosphodiesterases (PDE) degrade cAMP'''
    Parameter('kPDE4_cAMP_1', 0.01*1E3)# 1/(μM*s) |
    Parameter('kPDE4_cAMP_2', 1)# 1/s      |
    Rule('reaction73', PDE4(PDE4_s='i') + cAMP() | PDE4(PDE4_s='a'), kPDE4_cAMP_1, kPDE4_cAMP_2)

    Parameter('kPDE4_cAMP_to_AMP', 2)
    Rule('reaction74', PDE4(PDE4_s='a') >> PDE4(PDE4_s='i') + AMP(), kPDE4_cAMP_to_AMP)

    Parameter('kPDE10_2cAMP_1', 1.0E-6*1E3)# 1/(μM*s) |
    Parameter('kPDE10_2cAMP_2', 9)# 1/s
    Rule('reaction75', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() + cAMP() | PDE10(PDE10_c='Y', PDE10_s='i'), kPDE10_2cAMP_1, kPDE10_2cAMP_2)

    Parameter('kPDE10_cAMP_1', 0.1*1E3)# 1/(μM*s) |
    Parameter('kPDE10_cAMP_2', 2)# 1/s      |
    Rule('reaction76', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() | PDE10(PDE10_c='N', PDE10_s='a'),  kPDE10_cAMP_1, kPDE10_cAMP_2)

    Parameter('kPDE10_2cAMP_cAMP_1', 0.13*1E3)# 1/(μM*s) |
    Parameter('kPDE10_2cAMP_cAMP_2', 2)# 1/s      |
    Rule('reaction77', PDE10(PDE10_c='Y', PDE10_s='i') + cAMP() | PDE10(PDE10_c='Y', PDE10_s='a'), kPDE10_2cAMP_cAMP_1, kPDE10_2cAMP_cAMP_2)

    Parameter('kPDE10_cAMP_decay', 3)# 1/s      |
    Rule('reaction78', PDE10(PDE10_c='N', PDE10_s='a') >> PDE10(PDE10_c='N', PDE10_s='i') + AMP(), kPDE10_cAMP_decay)

    Parameter('kPDE10_2cAMP_cAMP_decay', 10)# 1/s      |
    Rule('reaction79', PDE10(PDE10_c='Y', PDE10_s='a') >> PDE10(PDE10_c='Y', PDE10_s='i') + AMP(), kPDE10_2cAMP_cAMP_decay)

    '''cAMP Activates PKA'''
    Parameter('kPKA_cAMP2_1', 0.00026*1E3)# 1/(μM*s) |
    Parameter('kPKA_cAMP2_2', 1)# 1/s |
    Rule('reaction80', PKA(PKA_s='cAMP0') + cAMP() + cAMP() | PKA(PKA_s='cAMP2'), kPKA_cAMP2_1, kPKA_cAMP2_2)

    Parameter('kPKA_cAMP4_1', 0.000346*1E3)# 1/(μM*s) |
    Parameter('kPKA_cAMP4_2', 1)# 1/s |
    Rule('reaction81', PKA(PKA_s='cAMP2') + cAMP() + cAMP() | PKA(PKA_s='cAMP4'), kPKA_cAMP4_1, kPKA_cAMP4_2)

    Parameter('kPKA_activation', 10*1E3)# 1/(μM*s) |
    Parameter('kPKA_activation_reverse', 0.01)# 1/s |
    Rule('reaction82', PKA(PKA_s='cAMP4') | PKAc(PKAc_b1=None) + PKAreg(), kPKA_activation, kPKA_activation_reverse)


    # OBSERVABLES
    Observable('obs_R',R(R_b1=None, R_p='p0', R_s='i'))
    Observable('obs_RL',R(R_b1=None, R_p='p0', R_s='a'))
    Observable('obs_Gi',Gi(Gi_b1=None))
    Observable('obs_GaiGTP',GaiGTP(GaiGTP_b1=None))
    Observable('obs_Gbgi',Gbgi())
    Observable('obs_GaiGDP', GaiGDP())
    Observable('obs_AC5',AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'))
    Observable('obs_AC5_Ca',AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'))
    Observable('obs_Ca',Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'))
    Observable('obs_ATP', ATP())
    Observable('obs_cAMP', cAMP())
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

list_of_observables=['obs_R','obs_RL','obs_cAMP', 'obs_AC5']
