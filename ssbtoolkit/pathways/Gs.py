##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                        G alpha s coupled receptor                            ##
##                                                                              ##
##                      PATHWAY based on D1 receptor                            ##
##                                                                              ##
##                    last modification 26 October 2021                         ##
##                                                                              ##
##   References:                                                                ##
##   Nair, A. G. et al (2015). Sensing Positive versus Negative Reward Signals  ##
##   through Adenylyl Cyclase-Coupled GPCRs in Direct and Indirect Pathway      ##
##   Striatal Medium Spiny Neurons. The Journal of Neuroscience : the Official  ##
##   Journal of the Society for Neuroscience, 35(41), 14017–14030.              ##
##################################################################################

from pysb import *
from pysb.macros import *
from sympy import Piecewise
from pysb.macros import create_t_obs

defaultParameters = {
        'time_in':0,
        'time_out':0,
        'L_init':0.01,
        'R_init': 2,
        'Golf_init': 2,
        'AC5_init': 0.7,
        'Ca_init': 0.06,
        'ATP_init': 5000,
        'PDE4_init': 2,
        'PDE10_init': 0.7,
        'PKA_init': 1.2,
        'RL_kon':0.005*1E3,
        'RL_koff':5,
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

def network(LR=None, kinetics=True, **kwargs):

    parameters={**defaultParameters, **kwargs}
 
    def myeval(x):
        try:
            y = eval(x)
        except:
            y=x
        return y

    parameters = dict(zip(parameters.keys(), map(myeval, parameters.values())))


    #Start a model
    Model()

    ##TIME OBSERVABLE
    components_time_obs = create_t_obs()
    time_obs = components_time_obs.t

    """IMPORTANT INFO:
        
        The receptor MONOMER is represented by:
            * 1 binding site;
            * 1 phosporilation site;
            * 1 state site 

        LR is represented by the state site: R(R_b1=None, R_p='p0', R_s='a'), REACTION1
        Golf binds to the biding site: R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50), REACTION3
        Arrestin acts on the phosporilate site (not described in this short network).
        
    """
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
        if parameters['time_in'] !=0 and parameters['time_out']!=0:
            Expression('RL_kon', Piecewise((0, time_obs > parameters['time_out']),(parameters['RL_kon'], time_obs > parameters['time_in']),(0, True)))
            Parameter('RL_koff', parameters['RL_koff'])    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
            Rule('reaction1', R(R_b1=None, R_p='p0', R_s='i') + L(L_b1=None) | R(R_b1=None, R_p='p0', R_s='a'), RL_kon, RL_koff)
        else:
            Parameter('RL_kon', parameters['RL_kon'])  # 1/(μM*s) |Association constant of the complex D1R_DA_Golf
            Parameter('RL_koff', parameters['RL_koff'])    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
            Rule('reaction1', R(R_b1=None, R_p='p0', R_s='i') + L(L_b1=None) | R(R_b1=None, R_p='p0', R_s='a'), RL_kon, RL_koff)
    else: pass

    Parameter('RL_Golf_kon', parameters['RL_Golf_kon'])  # 1/(μM*s) |Association constant of the complex D1R_DA_Golf
    Parameter('RL_Golf_koff', parameters['RL_Golf_koff'])    # 1/s      |Dissociation constant of the complex D1R_DA_Golf
    Rule('reaction3', R(R_b1=None, R_p='p0', R_s='a') + Golf(Golf_b1=None) | R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50), RL_Golf_kon, RL_Golf_koff)

    Parameter('RL_Golf_decay', parameters['RL_Golf_decay']) # 1/s      |Rate of formation of M4R_DA + Golf
    Rule('reaction9', R(R_b1=50, R_p='p0', R_s='a') % Golf(Golf_b1=50) >> R(R_b1=None, R_p='p0', R_s='a') + Gbgolf() + GaolfGTP(GaolfGTP_b1=None) , RL_Golf_decay)

    Parameter('GaolfGTP_decay', parameters['GaolfGTP_decay'])    # 1/s      |Rate of convertion of GaolfGTP into GaolfGDP
    Rule('reaction10', GaolfGTP(GaolfGTP_b1=None) >> GaolfGDP(), GaolfGTP_decay)

    Parameter('Golf_formation', parameters['Golf_formation'])   # 1/s      |Rate of formation of Golf
    Rule('reaction11', GaolfGDP() + Gbgolf() >> Golf(Golf_b1=None), Golf_formation)

    '''CYCLIC AMP FORMATION AND PKA ACTIVATION:
    - active GaolfGTP binds to and activate cyclase type V (AC5), which produces cyclic AMP (cAMP)
    - Ca2+ causes a 50% reduction in the activity of AC5 before or after binding to Gaolf. In the model
        calcium is allowed to bind to the inactive AC5, which then can be activated by GaolfGTP.
    - A step of regeneration of ATP is included that allows for a steady state concentration of ATP in
    the absence of stimulation. It is not assumed to be a rate-limiting substrate.
    '''

    # AC5 basal activity
    Parameter('AC5_ATP_kon', parameters['AC5_ATP_kon'])
    Parameter('AC5_ATP_koff', parameters['AC5_ATP_koff'])        
    Rule('reaction21', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') , AC5_ATP_kon, AC5_ATP_koff)

    Parameter('AC5_basal', parameters['AC5_basal'])          # 1/s      |Rate constant of basal formation of cAMP
    Rule('reaction22', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP(), AC5_basal)

    Parameter('AC5_reverse_basal', parameters['AC5_reverse_basal']) # 1/s  |Reverse Rate constant of basal formation of cAMP
    Rule('reaction23', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'), AC5_reverse_basal )

    # Interaction between AC5 and Ca2+
    Parameter('AC5_Ca_kon',parameters['AC5_Ca_kon'])       # 1/(μM*s) |Association constant of the complex AC5_Ca
    Parameter('AC5_Ca_koff', parameters['AC5_Ca_koff'])         # 1/s      |Dissociation constant of the complex AC5_Ca
    Rule('reaction16', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_kon, AC5_Ca_koff)

    Parameter('AC5_Ca_ATP_kon', parameters['AC5_Ca_ATP_kon']) # 1/(μM*s) |Association constant of the complex AC5_Ca_ATP
    Parameter('AC5_Ca_ATP_koff', parameters['AC5_Ca_ATP_koff'])       # 1/s      |Dissociation constant of the complex AC5_Ca_ATP
    Rule('reaction24', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_kon, AC5_Ca_ATP_koff)

    Parameter('AC5_Ca_ATP_to_cAMP', parameters['AC5_Ca_ATP_to_cAMP'] )  # 1/s      |Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction25', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') >> cAMP() + AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_to_cAMP)

    Parameter('AC5_Ca_ATP_to_cAMP_reverse', parameters['AC5_Ca_ATP_to_cAMP_reverse']) # 1/s      |Reverse Convertion rate of ATP into cAMP by AC5_Ca
    Rule('reaction26', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_to_cAMP_reverse)

    Parameter('AC5_ATP_Ca_kon', parameters['AC5_ATP_Ca_kon'])       # 1/(μM*s) |
    Parameter('AC5_ATP_Ca_koff', parameters['AC5_ATP_Ca_koff'])         # 1/s      |
    Rule('reaction30', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') , AC5_ATP_Ca_kon, AC5_ATP_Ca_koff)

    # Interaction between AC5 and Gaolf
    Parameter('AC5_GaolfGTP_kon', parameters['AC5_GaolfGTP_kon'] )   # 1/(μM*s) |Association constant of the complex AC5_GaolfGTP
    Parameter('AC5_GaolfGTP_koff', parameters['AC5_GaolfGTP_koff'])   # 1/s      |Dissociation constant of the complex AC5_GaolfGTP
    Rule('reaction15', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50), AC5_GaolfGTP_kon, AC5_GaolfGTP_koff)

    Parameter('AC5_GaolfGTP_ATP_i_kon', parameters['AC5_GaolfGTP_ATP_i_kon']) # 1/(μM*s) |Association constant of the complex AC5_GaolfGTP-ATP
    Parameter('AC5_GaolfGTP_ATP_i_koff', parameters['AC5_GaolfGTP_ATP_i_koff'])       # 1/s      |Dissociation constant of the complex AC5_GaolfGTP-ATP
    Rule('reaction18', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), AC5_GaolfGTP_ATP_i_kon, AC5_GaolfGTP_ATP_i_koff)

    Parameter('cAMP_formation', parameters['cAMP_formation'])    # 1/s      |Rate constant of convertion of  ATP into cAMP by AC5
    Rule('reaction19', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , cAMP_formation)

    Parameter('cAMP_reverse', parameters['cAMP_reverse'])   # 1/s      |Rate constant of reverse cAMP to ATP
    Rule('reaction20', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), cAMP_reverse)

    Parameter('AC5_GaolfGTP_ATP_a_kon', parameters['AC5_GaolfGTP_ATP_a_kon'])   # 1/(μM*s) |
    Parameter('AC5_GaolfGTP_ATP_a_koff', parameters['AC5_GaolfGTP_ATP_a_koff'])   # 1/s      |
    Rule('reaction27', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50), AC5_GaolfGTP_ATP_a_kon, AC5_GaolfGTP_ATP_a_koff)

    Parameter('AC5_GaolfGTP_decay', parameters['AC5_GaolfGTP_decay'])   # 1/s      |
    Rule('reaction28', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='i') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaolfGDP(), AC5_GaolfGTP_decay)

    Parameter('AC5_GaolfGTP_ATP_decay', parameters['AC5_GaolfGTP_ATP_decay'])# 1/s     |
    Rule('reaction29', AC5(AC5_b1=50, AC5_b2=None, AC5_b3=None, AC5_s='a') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaolfGDP(), AC5_GaolfGTP_ATP_decay)

    # Interaction between AC5, Ca2+, and Gaolf
    Parameter('AC5_Ca_GaolfGTP_kon', parameters['AC5_Ca_GaolfGTP_kon'])# 1/(μM*s) |Association constant of the complex AC5_Ca_GaolfGTP
    Parameter('AC5_Ca_GaolfGTP_koff', parameters['AC5_Ca_GaolfGTP_koff'])# 1/s      |Dissociation constant of the complex AC5_Ca_GaolfGTP
    Rule('reaction17', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), AC5_Ca_GaolfGTP_kon , AC5_Ca_GaolfGTP_koff)

    Parameter('AC5_Ca_GaolfGTP_ATP_kon', parameters['AC5_Ca_GaolfGTP_ATP_kon'])   
    Parameter('AC5_Ca_GaolfGTP_ATP_koff', parameters['AC5_Ca_GaolfGTP_ATP_koff'])   
    Rule('reaction31', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + ATP() | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50), AC5_Ca_GaolfGTP_ATP_kon, AC5_Ca_GaolfGTP_ATP_koff)

    Parameter('AC5_CA_ATP_GaolfGTP_kon', parameters['AC5_CA_ATP_GaolfGTP_kon']) 
    Parameter('AC5_CA_ATP_GaolfGTP_koff', parameters['AC5_CA_ATP_GaolfGTP_koff'])
    Rule('reaction32', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None) | AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , AC5_CA_ATP_GaolfGTP_kon, AC5_CA_ATP_GaolfGTP_koff)

    Parameter('AC5_Ca_GaolfGTP_ATP_to_cAMP', parameters['AC5_Ca_GaolfGTP_ATP_to_cAMP'])
    Rule('reaction33', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() , AC5_Ca_GaolfGTP_ATP_to_cAMP)

    Parameter('AC5_Ca_GaolfGTP_ATP_to_cAMP_reverse', parameters['AC5_Ca_GaolfGTP_ATP_to_cAMP_reverse']) 
    Rule('reaction34', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) + cAMP() >> AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) , AC5_Ca_GaolfGTP_ATP_to_cAMP_reverse )

    Parameter('AC5_Ca_GaolfGTP_decay', parameters['AC5_Ca_GaolfGTP_decay'])
    Rule('reaction35', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGTP(GaolfGTP_b1=None), AC5_Ca_GaolfGTP_decay)

    Parameter('AC5_Ca_GaolfGTP_ATP_decay', parameters['AC5_Ca_GaolfGTP_ATP_decay']) 
    Rule('reaction36', AC5(AC5_b1=50, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaolfGTP(GaolfGTP_b1=50) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaolfGDP(), AC5_Ca_GaolfGTP_ATP_decay)

    '''Several Subtypes of Phosphodiesterases (PDE) degrade cAMP'''
    Parameter('PDE4_cAMP_kon', parameters['PDE4_cAMP_kon'])
    Parameter('PDE4_cAMP_koff', parameters['PDE4_cAMP_koff'])
    Rule('reaction73', PDE4(PDE4_s='i') + cAMP() | PDE4(PDE4_s='a'), PDE4_cAMP_kon, PDE4_cAMP_koff)

    Parameter('PDE4_cAMP_to_AMP', parameters['PDE4_cAMP_to_AMP'])
    Rule('reaction74', PDE4(PDE4_s='a') >> PDE4(PDE4_s='i') + AMP(), PDE4_cAMP_to_AMP)

    Parameter('PDE10_2cAMP_kon', parameters['PDE10_2cAMP_kon'])
    Parameter('PDE10_2cAMP_koff', parameters['PDE10_2cAMP_koff'])
    Rule('reaction75', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() + cAMP() | PDE10(PDE10_c='Y', PDE10_s='i'), PDE10_2cAMP_kon, PDE10_2cAMP_koff)

    Parameter('PDE10_cAMP_kon', parameters['PDE10_cAMP_kon'])
    Parameter('PDE10_cAMP_koff',  parameters['PDE10_cAMP_koff'])
    Rule('reaction76', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() | PDE10(PDE10_c='N', PDE10_s='a'),  PDE10_cAMP_kon, PDE10_cAMP_koff)

    Parameter('PDE10_2cAMP_cAMP_kon', parameters['PDE10_2cAMP_cAMP_kon'])
    Parameter('PDE10_2cAMP_cAMP_koff', parameters['PDE10_2cAMP_cAMP_koff'])
    Rule('reaction77', PDE10(PDE10_c='Y', PDE10_s='i') + cAMP() | PDE10(PDE10_c='Y', PDE10_s='a'), PDE10_2cAMP_cAMP_kon, PDE10_2cAMP_cAMP_koff)

    Parameter('PDE10_cAMP_decay', parameters['PDE10_cAMP_decay'])
    Rule('reaction78', PDE10(PDE10_c='N', PDE10_s='a') >> PDE10(PDE10_c='N', PDE10_s='i') + AMP(), PDE10_cAMP_decay)

    Parameter('PDE10_2cAMP_cAMP_decay', parameters['PDE10_2cAMP_cAMP_decay'])
    Rule('reaction79', PDE10(PDE10_c='Y', PDE10_s='a') >> PDE10(PDE10_c='Y', PDE10_s='i') + AMP(), PDE10_2cAMP_cAMP_decay)

    '''cAMP Activates PKA'''
    Parameter('PKA_cAMP2_kon', parameters['PKA_cAMP2_kon'])
    Parameter('PKA_cAMP2_koff', parameters['PKA_cAMP2_koff'])
    Rule('reaction80', PKA(PKA_s='cAMP0') + cAMP() + cAMP() | PKA(PKA_s='cAMP2'), PKA_cAMP2_kon, PKA_cAMP2_koff)

    Parameter('PKA_cAMP4_kon', parameters['PKA_cAMP4_kon'])
    Parameter('PKA_cAMP4_koff', parameters['PKA_cAMP4_koff'])
    Rule('reaction81', PKA(PKA_s='cAMP2') + cAMP() + cAMP() | PKA(PKA_s='cAMP4'), PKA_cAMP4_kon, PKA_cAMP4_koff)

    Parameter('PKA_activation', parameters['PKA_activation'])
    Parameter('PKA_activation_reverse', parameters['PKA_activation_reverse'])
    Rule('reaction82', PKA(PKA_s='cAMP4') | PKAc(PKAc_b1=None) + PKAreg(), PKA_activation, PKA_activation_reverse)

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

list_of_observables=['obs_R','obs_RL','obs_cAMP', 'obs_AC5'] #do we need this?

