##################################################################################
##               THE G ALPHA Q11 COUPLED RECEPTOR CASCADE VERSION 1.0           ##
##                                                                              ##
##                      PATHWAY based on serotonin receptor                     ##
##                                                                              ##
##                    last modification 26 October 2021                         ##
##                                                                              ##
##   References:                                                                ##
##   Chang, Chiung-wen, Ethan Poteet, John A. Schetz, Zeynep H. Gümüş,          ##
##   and Harel Weinstein. 2009. “Towards a Quantitative Representation          ##
##   of the Cell Signaling Mechanisms of Hallucinogens: Measurement and         ##
##   Mathematical Modeling of 5-HT1A and 5-HT2A Receptor-Mediated ERK1/2        ##
##   Activation.” Neuropharmacology 56 (Suppl 1): 213–25.                       ##
##################################################################################


from pysb import *
from pysb.macros import *
from pysb.macros import create_t_obs
from sympy import Piecewise


__author__ = "Rui Ribeiro"
__organizarion__ = "University of Verona"
__copyright__ = "Copyright 2020, Rui Ribeiro"
__credits__ = ["Rui Ribeiro","Pietro Micheli"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"
__status__ = "Production"

defaultParameters = {
        'time_in':0,
        'time_out':0,
        'L_init':0,
        'R_init':1.4107,
        'Gq_a_GDP_init':0.0027739,
        'Gq_a_GDP_init':6.4172E-4,
        'Gq_bg_init':0.0037173,
        'RGS4_init':0.019994,
        'RGS4_Gq_a_GTP_init':6.4168E-6,
        'Ca_init':0.1,
        'PLCb_init':0.090022,
        'PLCb_Gq_a_GTP_init':1.4492E-4,
        'PLCb_Ca_init':0.0093825,
        'PLCb_Ca_Gq_a_GTP_init':1.5038E-4,
        'PIP2_init':2.6578,
        'IP3_init':0.21952,
        'DAG_init':0.055555,
        'R_Gq_trimer_init':0,
        'R_L_Gq_trimer_init':0,
        'Gq_trimer_init':0.61869,
        'RL_kon':1.00,
        'RL_koff':0.0046,
        'Gqa_Gqbg_kon':6.0,
        'Gqa_Gqbg_koff':0.0001,
        'RL_Gq_decay':0.04,
        'RGS4_Gq_a_GTP_kon': 20.83,
        'RGS4_Gq_a_GTP_koff': 33.32,
        'RGS4_Gq_a_GTP_decay': 8.33, 
        'Gq_a_GTP_decay':0.01,
        'Gq_a_GTP_PLCb_kon':2.52,
        'Gq_a_GTP_PLCb_off':1.00,
        'Gq_a_GTP_PLCb_Ca_kon':30.0,
        'Gq_a_GTP_PLCb_Ca_koff':1.00,
        'Gq_a_GTP_PLCb_Ca_decay':0.013,
        'PLCb_Ca_kon':3.00,
        'PLCb_Ca_koff':1.00,
        'PLCb_Ca_Gq_a_GTP_kon':25.2,
        'PLCb_Ca_Gq_a_GTP_koff':1.00,
        'DAG_decay':0.15,
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
            * 1 state site 

        LR is represented by the state site: R(R_b1=None, R_s='a'), REACTION1
        Golf binds to the biding site: R(R_b1=50, R_s='a') % Golf(Golf_b1=50), REACTION3
        
    """

    #MONOMERS
    Monomer('L', ['L_b1'])
    Monomer('R', ['R_b1', 'R_s'], {'R_s':['inact', 'act']})
    Monomer('Gq_a', ['Gq_a_b1', 'Gq_a_b2', 'Gq_a_s'], {'Gq_a_s' : ['GTP', 'GDP']})
    Monomer('Gq_bg', ['Gq_bg_b1', 'Gq_bg_b2'])
    Monomer('RGS4', ['RGS4_b1'])
    Monomer('Ca', ['Ca_b1'])
    Monomer('PLCb', ['PLCb_b1', 'PLCb_b2'])
    Monomer('PIP2', ['PIP2_b1'])
    Monomer('IP3', ['IP3_b1'])
    Monomer('DAG', ['DAG_b1', 'DAG_s'], {'DAG_s' : ['act', 'inact']})

    #INITIAL CONDITIONS
    if kinetics==True:
        Initial(R(R_b1=None, R_s='inact'), Parameter('R_0', parameters['R_init']))
        Initial(L(L_b1=None), Parameter('L_0', parameters['L_init']))
    else:
        Initial(R(R_b1=None, R_s='act'), Parameter('RL_0', LR))   
        Initial(R(R_b1=None, R_s='inact'), Parameter('R_0', parameters['R_init']-LR))
    
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Parameter('Gq_a_GDP_init', parameters['Gq_a_GDP_init']))
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('Gq_a_GDP_init', parameters['Gq_a_GDP_init']))
    Initial(Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), Parameter('Gq_bg_init', parameters['Gq_bg_init']))
    Initial(RGS4(RGS4_b1=None), Parameter('RGS4_init', parameters['RGS4_init']))
    Initial(RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('RGS4_Gq_a_GTP_init', parameters['RGS4_Gq_a_GTP_init']))
    Initial(Ca(Ca_b1=None), Parameter('Ca_init', parameters['Ca_init']))
    Initial(PLCb(PLCb_b1=None, PLCb_b2=None), Parameter('PLCb_init', parameters['0.090022']))
    Initial(PLCb(PLCb_b1=60, PLCb_b2=None)%Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('PLCb_Gq_a_GTP_init', parameters['PLCb_Gq_a_GTP_init']))
    Initial(PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70), Parameter('PLCb_Ca_init', parameters['PLCb_Ca_init']))
    Initial(Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70), Parameter('PLCb_Ca_Gq_a_GTP_init', parameters['PLCb_Ca_Gq_a_GTP_init']))
    Initial(PIP2(PIP2_b1=None), Parameter('PIP2_init', parameters['PIP2_init']))
    Initial(IP3(IP3_b1=None), Parameter('IP3_init', parameters['IP3_init']))
    Initial(DAG(DAG_b1=None, DAG_s='inact'), Parameter('DAG_init', parameters['DAG_init']))
    Initial(R(R_b1=30,  R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None),Parameter('R_Gq_trimer_init', parameters['R_Gq_trimer_init']))
    Initial(R(R_b1=30,  R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Parameter('R_L_Gq_trimer_init', parameters['R_L_Gq_trimer_init']))
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Parameter('Gq_trimer_init', parameters['Gq_trimer_init']))

    
    

    #Ligand
    #L_conc = 0.1
    
    #
    #
    #Observable('obs_L', L(L_b1=None))

    #OBSERVABLES (NEEDED FOR THE REACTIONS)
    
    Observable('obs_PIP2', PIP2(PIP2_b1=None))
    Observable('obs_IP3', IP3(IP3_b1=None))
    Observable('obs_DAG', DAG(DAG_b1=None, DAG_s='inact'))



################################################################################################################################################################################
    ##RULES

    if kinetics == True:
        if parameters['time_in'] !=0 and parameters['time_out']!=0:
            Expression('RL_kon', Piecewise((0, time_obs > parameters['time_out']),(parameters['RL_kon'], time_obs > parameters['time_in']),(0, True)))
            Parameter('RL_koff', parameters['RL_koff'])    # 1/s      |Dissociation constant of the complex 
            Rule('reaction1', R(R_b1=None, R_s='inact') + L(L_b1=None) | R(R_b1=None, R_s='act'), RL_kon, RL_koff)
        else:
            Parameter('RL_kon', parameters['RL_kon'])  # 1/(μM*s) |Association constant of the complex 
            Parameter('RL_koff', parameters['RL_koff'])    # 1/s      |Dissociation constant of the complex 
            Rule('reaction1', R(R_b1=None, R_s='inact') + L(L_b1=None) | R(R_b1=None, R_s='act'), RL_kon, RL_koff)
    else: pass

    #G-PROTEIN ACTIVATION
    #R+Gtrimer
    #Parameter('R_Gq_trimer_kf', 1.00)
    #Parameter('R_Gq_trimer_kr', 1.67)
    #Rule('R_Gq_trimer', R(R_b1=None, R_b2=None, R_s='inact') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_Gq_trimer_kf, R_Gq_trimer_kr)
    #Observable('obs_RG', R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))

    #RL+Gtrimer
    Parameter('R_L_Gq_trimer_kf', parameters['RL_kon'])
    Parameter('R_L_Gq_trimer_kr', parameters['RL_off'])
    Rule('R_L_Gq_trimer', R(R_b1=None, R_s='act') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_L_Gq_trimer_kf, R_L_Gq_trimer_kr)

    #Gq_trimerization
    Parameter('Gq_trimerization_k', parameters['Gqa_Gqbg_kon'])
    Parameter('Gq_trimer_split_k', parameters['Gqa_Gqbg_koff'])
    Rule('Gq_trimerization', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None) | Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimerization_k, Gq_trimer_split_k)

    #RL_Gq split
    Parameter('R_L_Gq_trimer_split_k', parameters['RL_Gq_decay'])
    Rule('R_L_Gq_trimer_split', R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) >> R(R_b1=None, R_s='act') + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), R_L_Gq_trimer_split_k)

    #Deactivation of Gq_a_GTP by RGS4
    #RGS4 + Gq_a_GTP
    Parameter('RGS4_Gq_a_GTP_kf', parameters['RGS4_Gq_a_GTP_kon'])
    Parameter('RGS4_Gq_a_GTP_kr', parameters['RGS4_Gq_a_GTP_koff'])
    Rule('RGS4_Gq_a_GTP', RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_kf, RGS4_Gq_a_GTP_kr)

    #RGS4_Gq_a_GTP dissociation
    Parameter('RGS4_Gq_a_GTP_diss_k', parameters['RGS4_Gq_a_GTP_decay'])
    Rule('RGS4_Gq_a_GTP_diss', RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP') >> RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), RGS4_Gq_a_GTP_diss_k)

    #Gq_a_GTP decay
    Parameter('Gq_a_GTP_decay_k', parameters['Gq_a_GTP_decay'])
    Rule('Gq_a_GTP_decay', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Gq_a_GTP_decay_k)

    #ACTIVATION OF PLCb
    #Gq_a_GTP + PLCb
    Parameter('Gq_a_GTP_PLCb_kf', parameters['Gq_a_GTP_PLCb_kon'])
    Parameter('Gq_a_GTP_PLCb_kr', parameters['Gq_a_GTP_PLCb_off'])
    Rule('Gq_a_GTP_PLCb', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + PLCb(PLCb_b1=None, PLCb_b2=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None), Gq_a_GTP_PLCb_kf, Gq_a_GTP_PLCb_kr)

    #Gq_a_GTP_PLCb + Ca
    Parameter('Gq_a_GTP_PLCb_Ca_kf', parameters['Gq_a_GTP_PLCb_Ca_kon'])
    Parameter('Gq_a_GTP_PLCb_Ca_kr', parameters['Gq_a_GTP_PLCb_Ca_koff'])
    Rule('Gq_a_GTP_PLCb_Ca', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None) + Ca(Ca_b1=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_kf, Gq_a_GTP_PLCb_Ca_kr)

    #Gq_a_GTP_PLCb_Ca DECAY
    Parameter('Gq_a_GTP_PLCb_Ca_decay_k', parameters['Gq_a_GTP_PLCb_Ca_decay'])
    Rule('Gq_a_GTP_PLCb_Ca_diss', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_decay_k)

    #PLCb + Ca
    Parameter('PLCb_Ca_kf', parameters['PLCb_Ca_kon'])
    Parameter('PLCb_Ca_kr', parameters['PLCb_Ca_koff'])
    Rule('PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=None) + Ca(Ca_b1=None) | PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_kf, PLCb_Ca_kr)

    #PLCb_Ca + Gq_a_GTP
    Parameter('PLCb_Ca_Gq_a_GTP_kf', parameters['PLCb_Ca_Gq_a_GTP_kon'])
    Parameter('PLCb_Ca_Gq_a_GTP_kr', parameters['PLCb_Ca_Gq_a_GTP_koff'])
    Rule('PLCb_Ca_Gq_a_GTP', PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_Gq_a_GTP_kf, PLCb_Ca_Gq_a_GTP_kr)

    #IP3 AND DAG PRODUCTION
    #PIP2_PLCb_Ca
    Expression('PIP2_PLCb_Ca_k', (10.0 * obs_PIP2)/(40.13 + obs_PIP2))
    Rule('PIP2_PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1) >> IP3(IP3_b1=None) + DAG(DAG_b1=None, DAG_s='inact') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PIP2_PLCb_Ca_k )

    #PIP2_Gq_a_GTP_PLCb_Ca
    Expression('PIP2_Gq_a_GTP_PLCb_Ca_k', (48.0 * obs_PIP2)/(5.00 + obs_PIP2))
    Rule('PIP2_Gq_a_GTP_PLCb_Ca', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> IP3(IP3_b1=None) + DAG(DAG_b1=None, DAG_s='inact') + Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), PIP2_Gq_a_GTP_PLCb_Ca_k )

    #IP3 degradation
    Expression('IP3_deg_k', (0.14*(obs_IP3 - IP3_init))/obs_IP3)
    Rule('IP3_deg', IP3(IP3_b1=None) >> None, IP3_deg_k)

    #DAG degradation
    Parameter('DAG_deg_k', parameters=['DAG_decay'])
    Rule('DAG_deg', DAG(DAG_b1=None, DAG_s='inact') >> None, DAG_deg_k)

    #OBSERVABLES
    Observable('obs_R', R(R_b1=None, R_s='inact'))
    Observable('obs_RL', R(R_b1=None, R_s='act'))
    Observable('obs_Gq_a_GDP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'))
    Observable('obs_Gq_a_GTP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'))
    Observable('obs_Gq_bg', Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None))
    Observable('obs_Ca', Ca(Ca_b1=None))
    Observable('obs_PLCb', PLCb(PLCb_b1=None, PLCb_b2=None))
    Observable('obs_PLCb_Ca_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70))
    Observable('obs_PLCb_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None))
    Observable('obs_PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70))
    Observable('obs_trimer', R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))
    if kinetics==True:
        Observable('obs_L', L(L_b1=None))

    return model

list_of_observables=['obs_R','obs_RL','obs_IP3', 'obs_PLCb'] #do we need this?
