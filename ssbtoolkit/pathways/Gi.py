##################################################################################
##                 THE DOPAMINE RECEPTORS CASCADE VERSION 1.0                   ##
##                                                                              ##
##                        G alpha i coupled receptor                            ##
##                                                                              ##
##                        PATHWAY based on D2 receptor                          ##
##                                                                              ##
##                        last modification 21 Oct 2021                         ##
##                                                                              ##
##   References:                                                                ##
##   Nair, A. G. et al (2015). Sensing Positive versus Negative Reward Signals  ##
##   through Adenylyl Cyclase-Coupled GPCRs in Direct and Indirect Pathway      ##
##   Striatal Medium Spiny Neurons. The Journal of Neuroscience : the Official  ##
##   Journal of the Society for Neuroscience, 35(41), 14017–14030.              ##                                                          ##
##################################################################################
from pysb import *
from pysb.macros import *
from sympy import Piecewise

defaultParameters = {
    'time_in':0,
    'time_out':0,
    'L_init':0.01,
    'R_init': 2,
    'Gi_init': 3,
    'AC5_init': 0.7,
    'Ca_cytos_free': 0.06,
    'ATP_init': 5000,
    'PDE4_init': 2,
    'PDE10_init': 1,
    'PKA_init': 1.2,        
    'RL_kon':0.1*1E3,
    'RL_koff':200,
    'RL_Gi_kon':6.6*1E3,
    'RL_Gi_koff':200,
    'RL_Gi_decay': 60,
    'GaiGTP_decay': 30,
    'Gi_formation': 100,
    'AC5_ATP_kon': 0.0001*1E3,
    'AC5_ATP_koff': 1,
    'AC5_basal': 1,
    'AC5_reverse_basal': 0.0004,
    'AC5_Ca_kon': 0.001*1E3,
    'AC5_Ca_koff': 0.9,
    'AC5_Ca_ATP_kon': 7.50E-5*1E3,
    'AC5_Ca_ATP_koff': 1,
    'AC5_Ca_ATP_to_cAMP': 0.5,
    'AC5_Ca_ATP_to_cAMP_reverse': 0.00015,
    'AC5_ATP_Ca_kon': 0.001*1E3,
    'AC5_ATP_Ca_koff': 0.9,
    'AC5_GaiGTP_kon': 50*1E3,
    'AC5_GaiGTP_koff': 5,
    'AC5_GaiGTP_ATP_kon': 6.25E-5*1E3,
    'AC5_GaiGTP_ATP_koff': 1,
    'AC5_ATP_GaiGTP_kon': 50*1E3,
    'AC5_ATP_GaiGTP_koff': 5,
    'AC5_GaiGTP_ATP_to_cAMP': 0.25,
    'AC5_GaiGTP_ATP_to_cAMP_reverse': 0.00105,
    'AC5_GaiGTP_decay': 30,
    'AC5_GaiGTP_decay_koff': 30,
    'AC5_Ca_GaiGTP_kon': 50*1E3,
    'AC5_Ca_GaiGTP_koff': 5,
    'AC5_Ca_GaiGTP_ATP_kon': 5.63E-5*1E3,
    'AC5_Ca_GaiGTP_ATP_koff': 1,
    'AC5_Ca_ATP_GaiGTP_kon': 50*1E3,
    'AC5_Ca_ATP_GaiGTP_koff': 5,
    'AC5_Ca_GaiGTP_ATP_to_cAMP': 0.125,
    'AC5_Ca_GaiGTP_ATP_to_cAMP_reverse': 2.81E-5,
    'AC5_Ca_GaiGTP_decay': 30,
    'AC5_Ca_GaiGTP_ATP_decay': 30,
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
    'PKA_activation_reverse': 0.025
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


    #INITIAL CONDITIONS
    if kinetics==True:
        Initial(R(R_b1=None, R_p='p0', R_s='i'), Parameter('R_0', parameters['R_init'])) 
        Initial(L(L_b1=None), Parameter('L_0', parameters['L_init']))   
    else:
        Initial(R(R_b1=None, R_p='p0', R_s='a'), Parameter('RL_0', LR))   
        Initial(R(R_b1=None, R_p='p0', R_s='i'), Parameter('R_0', parameters['R_init']-LR))  

    
    Initial(Gi(Gi_b1=None), )                  
    Initial(AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i'), Parameter('AC5_0', parameters['AC5_init']))  
    Initial(Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free'), Parameter('Ca_cytos_free', parameters['Ca_cytos_free'])) 
    Initial(ATP(), Parameter('ATP_0', parameters['ATP_init']))
    Initial(PDE4(PDE4_s='i'), Parameter('PDE4_0', parameters['PDE4_init']))
    Initial(PDE10(PDE10_c='N', PDE10_s='i'), Parameter('PDE10_0', parameters['PDE10_init']))
    Initial(PKA(PKA_s='cAMP0'), Parameter('PKA_0', parameters['PKA_init']))

    #REACTIONS
    ##Dopamine and the G-protein Coupled RECEPTORS
    '''The R can bind either the inactive G protein first, and then
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

    Parameter('RL_Gi_kon', parameters['RL_Gi_kon'])
    Parameter('RL_Gi_koff', parameters['RL_Gi_koff'])
    Rule('reaction3', R(R_b1=None, R_p='p0', R_s='a') + Gi(Gi_b1=None) | R(R_b1=50, R_p='p0', R_s='a') % Gi(Gi_b1=50), RL_Gi_kon, RL_Gi_koff)

    Parameter('RL_Gi_decay', parameters['RL_Gi_decay'])
    Rule('reaction9', R(R_b1=50, R_p='p0', R_s='a') % Gi(Gi_b1=50) >> R(R_b1=None, R_p='p0', R_s='a') + Gbgi() + GaiGTP(GaiGTP_b1=None) , RL_Gi_decay)

    Parameter('GaiGTP_decay', parameters['GaiGTP_decay'])
    Rule('reaction13', GaiGTP(GaiGTP_b1=None) >> GaiGDP(), GaiGTP_decay)

    Parameter('Gi_formation', parameters['kGi_formation'])
    Rule('reaction14', Gbgi() + GaiGDP() >> Gi(Gi_b1=None), Gi_formation)

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

    Parameter('AC5_basal', parameters['AC5_basal'])
    Rule('reaction22', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP(), AC5_basal)

    Parameter('AC5_reverse_basal', parameters['AC5_reverse_basal'])
    Rule('reaction23', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a'), AC5_reverse_basal)

    # Interaction between AC5 and Ca2+
    Parameter('AC5_Ca_kon', parameters['AC5_Ca_kon'])
    Parameter('AC5_Ca_koff', parameters['AC5_Ca_koff'])
    Rule('reaction16', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_kon, AC5_Ca_koff)

    Parameter('AC5_Ca_ATP_kon', parameters['AC5_Ca_ATP_kon'])
    Parameter('AC5_Ca_ATP_koff', parameters['AC5_Ca_ATP_koff'])
    Rule('reaction24', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_kon, AC5_Ca_ATP_koff)

    Parameter('AC5_Ca_ATP_to_cAMP', parameters['AC5_Ca_ATP_to_cAMP'])
    Rule('reaction25', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') >> cAMP() + AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_to_cAMP)

    Parameter('AC5_Ca_ATP_to_cAMP_reverse', parameters['AC5_Ca_ATP_to_cAMP_reverse'])
    Rule('reaction26', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff'), AC5_Ca_ATP_to_cAMP_reverse)

    Parameter('AC5_ATP_Ca_kon', parameters['AC5_ATP_Ca_kon'])
    Parameter('AC5_ATP_Ca_koff', parameters['AC5_ATP_Ca_koff'])
    Rule('reaction30', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + Ca(Ca_b1=None, Ca_l='cytos', Ca_s='free') | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') , AC5_ATP_Ca_kon, AC5_ATP_Ca_koff)

    # Interaction between AC5 and Gai
    Parameter('AC5_GaiGTP_kon', parameters['AC5_GaiGTP_kon'])
    Parameter('AC5_GaiGTP_koff', parameters['AC5_GaiGTP_koff'])
    Rule('reaction37', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) , AC5_GaiGTP_kon, AC5_GaiGTP_koff)

    Parameter('AC5_GaiGTP_ATP_kon', parameters['AC5_GaiGTP_ATP_kon'])
    Parameter('AC5_GaiGTP_ATP_koff', parameters['AC5_GaiGTP_ATP_koff'])
    Rule('reaction38', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), AC5_GaiGTP_ATP_kon, AC5_GaiGTP_ATP_koff)

    Parameter('AC5_ATP_GaiGTP_kon', parameters['AC5_ATP_GaiGTP_kon'])
    Parameter('AC5_ATP_GaiGTP_koff', parameters['AC5_ATP_GaiGTP_koff'])
    Rule('reaction39', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), AC5_ATP_GaiGTP_kon, AC5_ATP_GaiGTP_koff)

    Parameter('AC5_GaiGTP_ATP_to_cAMP', parameters['AC5_GaiGTP_ATP_to_cAMP'])
    Rule('reaction40', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP(), AC5_GaiGTP_ATP_to_cAMP)

    Parameter('AC5_GaiGTP_ATP_to_cAMP_reverse', parameters['AC5_GaiGTP_ATP_to_cAMP_reverse'])
    Rule('reaction41', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10), AC5_GaiGTP_ATP_to_cAMP_reverse)

    Parameter('AC5_GaiGTP_decay', parameters['AC5_GaiGTP_decay'])
    Rule('reaction42', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='i') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='i') + GaiGDP(), AC5_GaiGTP_decay)

    Parameter('AC5_GaiGTP_decay_koff', parameters['AC5_GaiGTP_decay_koff'])
    Rule('reaction43', AC5(AC5_b1=None, AC5_b2=None, AC5_b3=10, AC5_s='a') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=None, AC5_b3=None, AC5_s='a') + GaiGDP(), AC5_GaiGTP_decay_koff)

    # Interaction between AC5, Ca2+, and Gai
    Parameter('AC5_Ca_GaiGTP_kon', parameters['AC5_Ca_GaiGTP_kon'])
    Parameter('AC5_Ca_GaiGTP_koff', parameters['AC5_Ca_GaiGTP_koff'])
    Rule('reaction55', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), AC5_Ca_GaiGTP_kon , AC5_Ca_GaiGTP_koff)

    Parameter('AC5_Ca_GaiGTP_ATP_kon', parameters['AC5_Ca_GaiGTP_ATP_kon'])
    Parameter('AC5_Ca_GaiGTP_ATP_koff', parameters['AC5_Ca_GaiGTP_ATP_koff'])
    Rule('reaction56', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + ATP() | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), AC5_Ca_GaiGTP_ATP_kon , AC5_Ca_GaiGTP_ATP_koff )

    Parameter('AC5_Ca_ATP_GaiGTP_kon', parameters['AC5_Ca_ATP_GaiGTP_kon'])
    Parameter('AC5_Ca_ATP_GaiGTP_koff', parameters['AC5_Ca_ATP_GaiGTP_koff'])
    Rule('reaction57', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGTP(GaiGTP_b1=None) | AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), AC5_Ca_ATP_GaiGTP_kon , AC5_Ca_ATP_GaiGTP_koff)

    Parameter('AC5_Ca_GaiGTP_ATP_to_cAMP', parameters['AC5_Ca_GaiGTP_ATP_to_cAMP'])
    Rule('reaction58', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP(), AC5_Ca_GaiGTP_ATP_to_cAMP)

    Parameter('AC5_Ca_GaiGTP_ATP_to_cAMP_reverse', parameters['AC5_Ca_GaiGTP_ATP_to_cAMP_reverse'])
    Rule('reaction59', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) + cAMP() >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10), AC5_Ca_GaiGTP_ATP_to_cAMP_reverse)

    Parameter('AC5_Ca_GaiGTP_decay', parameters['AC5_Ca_GaiGTP_decay'])
    Rule('reaction60', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='i') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), AC5_Ca_GaiGTP_decay)

    Parameter('AC5_Ca_GaiGTP_ATP_decay', parameters['AC5_Ca_GaiGTP_ATP_decay'])
    Rule('reaction61', AC5(AC5_b1=None, AC5_b2=20, AC5_b3=10, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') % GaiGTP(GaiGTP_b1=10) >> AC5(AC5_b1=None, AC5_b2=20, AC5_b3=None, AC5_s='a') % Ca(Ca_b1=20, Ca_l='cytos', Ca_s='buff') + GaiGDP(), AC5_Ca_GaiGTP_ATP_decay)

    '''Several Subtypes of Phosphodiesterases (PDE) degrade cAMP'''
    Parameter('PDE104_cAMP_kon', parameters['PDE104_cAMP_kon'])
    Parameter('PDE4_cAMP_koff', parameters['PDE4_cAMP_koff'])
    Rule('reaction73', PDE4(PDE4_s='i') + cAMP() | PDE4(PDE4_s='a'), PDE4_cAMP_kon, PDE4_cAMP_koff)

    Parameter('PDE4_cAMP_to_AMP', parameters['PDE4_cAMP_to_AMP'])
    Rule('reaction75', PDE4(PDE4_s='a') >> PDE4(PDE4_s='i') + AMP(), PDE4_cAMP_to_AMP)

    Parameter('PDE10_2cAMP_kon', parameters['PDE10_2cAMP_kon'])
    Parameter('PDE10_2cAMP_koff', parameters['PDE10_2cAMP_koff'])   
    Rule('reaction74', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() + cAMP() | PDE10(PDE10_c='Y', PDE10_s='i'), PDE10_2cAMP_kon, PDE10_2cAMP_koff)

    Parameter('PDE10_cAMP_kon', parameters['PDE10_cAMP_kon'])
    Parameter('PDE10_cAMP_koff', parameters['PDE10_cAMP_koff'])
    Rule('reaction76', PDE10(PDE10_c='N', PDE10_s='i') + cAMP() | PDE10(PDE10_c='N', PDE10_s='a'),  PDE10_cAMP_kon, PDE10_cAMP_koff)

    Parameter('PDE10_2cAMP_cAMP_kon', parameters['PDE10_2cAMP_cAMP_kon'])
    Parameter('PDE10_2cAMP_cAMP_koff', parameters['PDE10_2cAMP_cAMP_koff'])
    Rule('reaction77', PDE10(PDE10_c='Y', PDE10_s='i') + cAMP() | PDE10(PDE10_c='Y', PDE10_s='a'), PDE10_2cAMP_cAMP_kon, PDE10_2cAMP_cAMP_koff)

    Parameter('PDE10_cAMP_decay', parameters[''])
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
    if kinetics==True:
        Observable('obs_L', L(L_b1=None))
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
