##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                        G alpha s coupled receptor                            ##
##                                                                              ##
##                      PATHWAY based on D1 receptor                            ##
##                                                                              ##
##                    last modification 22 June 2021                            ##
##                                                                              ##
##   References:                                                                ##
##   Nair, A. G. et al (2015). Sensing Positive versus Negative Reward Signals  ##
##   through Adenylyl Cyclase-Coupled GPCRs in Direct and Indirect Pathway      ##
##   Striatal Medium Spiny Neurons. The Journal of Neuroscience : the Official  ##
##   Journal of the Society for Neuroscience, 35(41), 14017–14030.              ##
##################################################################################

from pysb import *
from pysb.macros import *


def network(LR=None, kinetics=True, **kwargs):

    defaultKwargs = {
        'L_init':0.01,
        'R_init': 2,
        'Golf_init': 2,
        'AC5_init': 0.7,
        'Ca_init': 0.06,
        'ATP_init': 5000,
        'PDE4_init': 2,
        'PDE10_init': 0.7,
        'PKA_init': 1.2,
        'RL_kon':0.35,
        'RL_koff':0.00055,
        'RL_Golf_kon': 0.003*1E3,
        'RL_Golf_koff': 5.0,
        'RL_Golf_decay': 15,
        'GaolfGTP_decay': 30,
        'Golf_formation': 100,
        'AC5_ATP_kon': 0.0001*1E3,
        'AC5_ATP_koff': 1,
        'AC5_basal': 1,
        'AC5_reverse_basal': 0.0004,
        'AC5_Ca_kon': 0.001*1E3,
        'AC5_Ca_koff': 0.9,
        'AC5_Ca_ATP_kon':  7.50E-5*1E3,
        'AC5_Ca_ATP_koff': 1,
        'AC5_Ca_ATP_to_cAMP': 0.5,
        'AC5_Ca_ATP_to_cAMP_reverse': 0.00015,
        'AC5_ATP_Ca_kon': 0.001*1E3,
        'AC5_ATP_Ca_koff': 0.9,
        'AC5_GaolfGTP_kon': 0.2*1E3,
        'AC5_GaolfGTP_koff': 0.1,
        'AC5_GaolfGTP_ATP_i_kon': 0.00105*1E3,
        'AC5_GaolfGTP_ATP_i_koff': 1,
        'cAMP_formation': 20,
        'cAMP_reverse': 0.084,
        'AC5_GaolfGTP_ATP_a_kon': 0.2*1E3,
        'AC5_GaolfGTP_ATP_a_koff': 0.1,
        'AC5_GaolfGTP_decay': 0.2,
        'AC5_GaolfGTP_ATP_decay': 0.2,
        'AC5_Ca_GaolfGTP_kon': 0.2*1E3,
        'AC5_Ca_GaolfGTP_koff': 0.1,
        'AC5_Ca_GaolfGTP_ATP_kon': 0.00055*1E3,
        'AC5_Ca_GaolfGTP_ATP_koff': 1,
        'AC5_CA_ATP_GaolfGTP_kon': 0.2*1E3,
        'AC5_CA_ATP_GaolfGTP_koff': 0.1,
        'AC5_Ca_GaolfGTP_ATP_to_cAMP': 10,
        'AC5_Ca_GaolfGTP_ATP_to_cAMP_reverse': 0.022,
        'AC5_Ca_GaolfGTP_decay': 0.2,
        'AC5_Ca_GaolfGTP_ATP_decay': 0.2,
        'PDE4_cAMP_kon': 0.01*1E3,
        'PDE4_cAMP_koff': 1,
        'PDE4_cAMP_to_AMP': 2,
        'PDE10_2cAMP_kon': 1.0E-6*1E3,
        'PDE10_2cAMP_koff': 9,
        'PDE10_cAMP_kon': 0.1*1E3,
        'PDE10_cAMP_koff': 2,
        'PDE10_2cAMP_cAMP_kon': 0.13*1E3,
        'PDE10_2cAMP_cAMP_koff': 2,
        'PDE10_cAMP_decay': 3,
        'PDE10_2cAMP_cAMP_decay': 10,
        'PKA_cAMP2_kon': 0.00026*1E3,
        'PKA_cAMP2_koff': 1,
        'PKA_cAMP4_kon': 0.000346*1E3,
        'PKA_cAMP4_koff': 1,
        'PKA_activation': 10*1E3,
        'PKA_activation_reverse': 0.01
    }
    parameters={**defaultKwargs, **kwargs}

    
    #Start a model
    Model()

    #MONOMERS
    Monomer('L', ['L_b1'])
    Monomer('R', ['R_b1', 'R_p', 'R_s'], {'R_p':['p0','p1'], 'R_s':['i','a']})    # Dopamine receptor 1 with two binding sites: for ligand and G-proteins
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
    if kinetics==True:
        Initial(R(R_b1=None, R_p='p0', R_s='i'), Parameter('R_0', parameters['R_init'])) 
        Initial(L(L_b1=None), Parameter('L_0', parameters['L_init']))   
    else:
        Initial(R(R_b1=None, R_p='p0', R_s='a'), Parameter('RL_0', LR))   
        Initial(R(R_b1=None, R_p='p0', R_s='i'), Parameter('R_0', parameters['R_init']-LR))   
    
    Initial(Golf(Golf_b1=None), Parameter('Golf_0', parameters['Golf_init']))             
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', parameters['AC5_init']))  
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', parameters['Ca_init'])) 
    Initial(ATP(), Parameter('ATP_0', parameters['ATP_init']))
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', parameters['PDE4_init']))
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', parameters['PDE10_init']))
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', parameters['PKA_init']))


    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The receptor can bind either the inactive G protein first, and then
    dopamine, or dopamine first and then and then the inactivate G protein.'''
    if kinetics == True:
        Parameter('kRL_1', parameters['RL_kon'])  # 1/(μM*s) |Association constant of the complex D1R_DA_Golf
        Parameter('kRL_2', parameters['RL_koff'])    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
        Rule('reaction1', R(R_b1=None, R_p='p0', R_s='i') + L(L_b1=None) | R(R_b1=None, R_p='p0', R_s='a'), kRL_1, kRL_2)
    else: pass

    Parameter('kRL_Golf_1', parameters['RL_Golf_kon'])  # 1/(μM*s) |Association constant of the complex D1R_DA_Golf
    Parameter('kRL_Golf_2', parameters['RL_Golf_koff'])    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
    Rule('reaction3', R(R_b1=None, R_p='p0', R_s='a') + Golf(Golf_b1=None) | R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50), kRL_Golf_1, kRL_Golf_2)

    Parameter('kRL_Golf_decay', parameters['RL_Golf_decay']) # 1/s      |Rate of formation of M4R_DA + Golf
    Rule('reaction9', R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50) >> R(R_b1=None, R_p='p0', R_s='a') + Gbgolf() + GaolfGTP(GaolfGTP_b1=None) , kRL_Golf_decay)

    Parameter('kGaolfGTP_decay', parameters['GaolfGTP_decay'])    # 1/s      |Rate of convertion of GaolfGTP into GaolfGDP
    Rule('reaction10', GaolfGTP(GaolfGTP_b1=None) >> GaolfGDP(), kGaolfGTP_decay)

    Parameter('kGolf_formation', parameters['Golf_formation'])   # 1/s      |Rate of formation of Golf
    Rule('reaction11', GaolfGDP() + Gbgolf() >> Golf(Golf_b1=None), kGolf_formation)

    '''CYCLIC AMP FORMATION AND PKA ACTIVATION:
    - active GaolfGTP binds to and activate cyclase type V (AC5), which produces cyclic AMP (cAMP)
    - Ca2+ causes a 50% reduction in the activity of AC5 before or after binding to Gaolf. In the model
        calcium is allowed to bind to the inactive AC5, which then can be activated by GaolfGTP.
    - A step of regeneration of ATP is included that allows for a steady state concentration of ATP in
    the absence of stimulation. It is not assumed to be a rate-limiting substrate.
    '''

    # AC5 basal activity
    Parameter('kAC5_ATP_1', parameters['AC5_ATP_kon'])
    Parameter('kAC5_ATP_2', parameters['AC5_ATP_koff'])        
    Rule('reaction21', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') , kAC5_ATP_1, kAC5_ATP_2)

    Parameter('kAC5_basal', parameters['AC5_basal'])          # 1/s      |Rate constant of basal formation of cAMP
    Rule('reaction22', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP(), kAC5_basal)

    Parameter('kAC5_reverse_basal', parameters['AC5_reverse_basal']) # 1/s  |Reverse Rate constant of basal formation of cAMP
    Rule('reaction23', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'), kAC5_reverse_basal )

    # Interaction between AC5 and Ca2+
    Parameter('kAC5_Ca_1',parameters['AC5_Ca_kon'])       # 1/(μM*s) |Association constant of the complex AC5_Ca
    Parameter('kAC5_Ca_2', parameters['AC5_Ca_koff'])         # 1/s      |Dissociation constant of the complex AC5_Ca
    Rule('reaction16', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_1, kAC5_Ca_2)

    Parameter('kAC5_Ca_ATP_1', parameters['AC5_Ca_ATP_kon']) # 1/(μM*s) |Association constant of the complex AC5_Ca_ATP
    Parameter('kAC5_Ca_ATP_2', parameters['AC5_Ca_ATP_koff'])       # 1/s      |Dissociation constant of the complex AC5_Ca_ATP
    Rule('reaction24', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_1, kAC5_Ca_ATP_2)

    Parameter('kAC5_Ca_ATP_to_cAMP', parameters['AC5_Ca_ATP_to_cAMP'] )  # 1/s      |Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction25', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') >> cAMP() + AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP)

    Parameter('kAC5_Ca_ATP_to_cAMP_reverse', parameters['AC5_Ca_ATP_to_cAMP_reverse']) # 1/s      |Reverse Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction26', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), kAC5_Ca_ATP_to_cAMP_reverse)

    Parameter('kAC5_ATP_Ca_1', parameters['AC5_ATP_Ca_kon'])       # 1/(μM*s) |
    Parameter('kAC5_ATP_Ca_2', parameters['AC5_ATP_Ca_koff'])         # 1/s      |
    Rule('reaction30', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') , kAC5_ATP_Ca_1, kAC5_ATP_Ca_2)

    # Interaction between AC5 and Gaolf
    Parameter('kAC5_GaolfGTP_1', parameters['AC5_GaolfGTP_kon'] )   # 1/(μM*s) |Association constant of the complex AC5_GaolfGTP
    Parameter('kAC5_GaolfGTP_2', parameters['AC5_GaolfGTP_koff'])   # 1/s      |Dissociation constant of the complex AC5_GaolfGTP
    Rule('reaction15', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_1, kAC5_GaolfGTP_2)

    Parameter('kAC5_GaolfGTP_ATP_1', parameters['AC5_GaolfGTP_ATP_i_kon']) # 1/(μM*s) |Association constant of the complex AC5_GaolfGTP-ATP
    Parameter('kAC5_GaolfGTP_ATP_2', parameters['AC5_GaolfGTP_ATP_i_koff'])       # 1/s      |Dissociation constant of the complex AC5_GaolfGTP-ATP
    Rule('reaction18', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_ATP_1, kAC5_GaolfGTP_ATP_2)

    Parameter('kcAMP_formation', parameters['cAMP_formation'])    # 1/s      |Rate constant of convertion of  ATP into cAMP by AC5
    Rule('reaction19', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , kcAMP_formation)

    Parameter('kcAMP_reverse', parameters['cAMP_reverse'])   # 1/s      |Rate constant of reverse cAMP to ATP
    Rule('reaction20', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kcAMP_reverse)

    Parameter('kAC5_GaolfGTP_ATP_3', parameters['AC5_GaolfGTP_ATP_a_kon'])   # 1/(μM*s) |
    Parameter('kAC5_GaolfGTP_ATP_4', parameters['AC5_GaolfGTP_ATP_a_koff'])   # 1/s      |
    Rule('reaction27', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), kAC5_GaolfGTP_ATP_3, kAC5_GaolfGTP_ATP_4)

    Parameter('kAC5_GaolfGTP_decay', parameters['AC5_GaolfGTP_decay'])   # 1/s      |
    Rule('reaction28', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGDP(), kAC5_GaolfGTP_decay)

    Parameter('kAC5_GaolfGTP_ATP_decay', parameters['AC5_GaolfGTP_ATP_decay'])# 1/s     |
    Rule('reaction29', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGDP(), kAC5_GaolfGTP_ATP_decay)

    # Interaction between AC5, Ca2+, and Gaolf
    Parameter('kAC5_Ca_GaolfGTP_1', parameters['AC5_Ca_GaolfGTP_kon'])# 1/(μM*s) |Association constant of the complex AC5_Ca_GaolfGTP
    Parameter('kAC5_Ca_GaolfGTP_2', parameters['AC5_Ca_GaolfGTP_koff'])# 1/s      |Dissociation constant of the complex AC5_Ca_GaolfGTP
    Rule('reaction17', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_1 , kAC5_Ca_GaolfGTP_2)

    Parameter('kAC5_Ca_GaolfGTP_ATP_1', parameters['AC5_Ca_GaolfGTP_ATP_kon'])   
    Parameter('kAC5_Ca_GaolfGTP_ATP_2', parameters['AC5_Ca_GaolfGTP_ATP_koff'])   
    Rule('reaction31', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), kAC5_Ca_GaolfGTP_ATP_1, kAC5_Ca_GaolfGTP_ATP_2)

    Parameter('kAC5_CA_ATP_GaolfGTP_1', parameters['AC5_CA_ATP_GaolfGTP_kon']) 
    Parameter('kAC5_CA_ATP_GaolfGTP_2', parameters['AC5_CA_ATP_GaolfGTP_koff'])
    Rule('reaction32', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_CA_ATP_GaolfGTP_1, kAC5_CA_ATP_GaolfGTP_2)

    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP', parameters['AC5_Ca_GaolfGTP_ATP_to_cAMP'])
    Rule('reaction33', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , kAC5_Ca_GaolfGTP_ATP_to_cAMP)

    Parameter('kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse', parameters['AC5_Ca_GaolfGTP_ATP_to_cAMP_reverse']) 
    Rule('reaction34', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , kAC5_Ca_GaolfGTP_ATP_to_cAMP_reverse )

    Parameter('kAC5_Ca_GaolfGTP_decay', parameters['AC5_Ca_GaolfGTP_decay'])
    Rule('reaction35', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None), kAC5_Ca_GaolfGTP_decay)

    Parameter('kAC5_Ca_GaolfGTP_ATP_decay', parameters['AC5_Ca_GaolfGTP_ATP_decay']) 
    Rule('reaction36', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGDP(), kAC5_Ca_GaolfGTP_ATP_decay)

    '''Several Subtypes of Phosphodiesterases (PDE) degrade cAMP'''
    Parameter('kPDE4_cAMP_1', parameters['PDE4_cAMP_kon'])
    Parameter('kPDE4_cAMP_2', parameters['PDE4_cAMP_koff'])
    Rule('reaction73', PDE4(PDE4_s='i') + cAMP() | PDE4(PDE4_s='a'), kPDE4_cAMP_1, kPDE4_cAMP_2)

    Parameter('kPDE4_cAMP_to_AMP', parameters['PDE4_cAMP_to_AMP'])
    Rule('reaction74', PDE4(PDE4_s='a') >> PDE4(PDE4_s='i') + AMP(), kPDE4_cAMP_to_AMP)

    Parameter('kPDE10_2cAMP_1', parameters['PDE10_2cAMP_kon'])
    Parameter('kPDE10_2cAMP_2', parameters['PDE10_2cAMP_koff'])
    Rule('reaction75', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() + cAMP() | PDE10(PDE10_c='Y', PDE10_s='i'), kPDE10_2cAMP_1, kPDE10_2cAMP_2)

    Parameter('kPDE10_cAMP_1', parameters['PDE10_cAMP_kon'])
    Parameter('kPDE10_cAMP_2',  parameters['PDE10_cAMP_koff'])
    Rule('reaction76', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() | PDE10(PDE10_c='N', PDE10_s='a'),  kPDE10_cAMP_1, kPDE10_cAMP_2)

    Parameter('kPDE10_2cAMP_cAMP_1', parameters['PDE10_2cAMP_cAMP_kon'])
    Parameter('kPDE10_2cAMP_cAMP_2', parameters['PDE10_2cAMP_cAMP_koff'])
    Rule('reaction77', PDE10(PDE10_c='Y', PDE10_s='i') + cAMP() | PDE10(PDE10_c='Y', PDE10_s='a'), kPDE10_2cAMP_cAMP_1, kPDE10_2cAMP_cAMP_2)

    Parameter('kPDE10_cAMP_decay', parameters['PDE10_cAMP_decay'])
    Rule('reaction78', PDE10(PDE10_c='N', PDE10_s='a') >> PDE10(PDE10_c='N', PDE10_s='i') + AMP(), kPDE10_cAMP_decay)

    Parameter('kPDE10_2cAMP_cAMP_decay', parameters['PDE10_2cAMP_cAMP_decay'])
    Rule('reaction79', PDE10(PDE10_c='Y', PDE10_s='a') >> PDE10(PDE10_c='Y', PDE10_s='i') + AMP(), kPDE10_2cAMP_cAMP_decay)

    '''cAMP Activates PKA'''
    Parameter('kPKA_cAMP2_1', parameters['PKA_cAMP2_kon'])
    Parameter('kPKA_cAMP2_2', parameters['PKA_cAMP2_koff'])
    Rule('reaction80', PKA(PKA_s='cAMP0') + cAMP() + cAMP() | PKA(PKA_s='cAMP2'), kPKA_cAMP2_1, kPKA_cAMP2_2)

    Parameter('kPKA_cAMP4_1', parameters['PKA_cAMP4_kon'])
    Parameter('kPKA_cAMP4_2', parameters['PKA_cAMP4_koff'])
    Rule('reaction81', PKA(PKA_s='cAMP2') + cAMP() + cAMP() | PKA(PKA_s='cAMP4'), kPKA_cAMP4_1, kPKA_cAMP4_2)

    Parameter('kPKA_activation', parameters['PKA_activation'])
    Parameter('kPKA_activation_reverse', parameters['PKA_activation_reverse'])
    Rule('reaction82', PKA(PKA_s='cAMP4') | PKAc(PKAc_b1=None) + PKAreg(), kPKA_activation, kPKA_activation_reverse)

    # OBSERVABLES
    
    Observable('obs_R',R(R_b1=None, R_p='p0', R_s='i'))
    Observable('obs_RL',R(R_b1=None, R_p='p0', R_s='a'))
    #Observable('obs_R_Golf',R(R_b1=None, R_b2=20) % Golf(Golf_b1=20))
    Observable('obs_RL_Golf',R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50))
    Observable('obs_cAMP', cAMP())
    Observable('obs_AC5', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'))
    if kinetics==True:
        Observable('obs_L', L(L_b1=None))
    return model

list_of_observables=['obs_R','obs_RL','obs_cAMP', 'obs_AC5']
