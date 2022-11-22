
#!/usr/bin/env python
import os
import sys

import math
import numpy as np
import sympy
from sympy import Piecewise
from pysb import *
from pysb.macros import *
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator
from pysb.macros import create_t_obs, drug_binding


defaultParameters = {
        'time_in':50,
        'time_out':51,
        'L_init':0.1,
        'R_init': 1.4107*0.73,
        'PIP2_init':2.6578 ,
        'Gq_a_GDP_init':0.0027739 ,
        'Gq_a_GTP_init': 6.4172E-4,
        'Gq_bg_init': 0.0037173,
        'RGS4_init': 0.019994,
        'RGS4_Gq_a_GTP_init':6.4168E-6 ,
        'Ca_init':0.1 ,
        'CaER_init':1,
        'PLCb_init':0.090022,
        'PLCb_Gq_a_GTP_init':1.4492E-4,
        'PLCb_Ca_init':0.0093825,
        'PLCb_Ca_Gq_a_GTP_init':1.5038E-4,
        'IP3_init':0.21952,
        'DAG_init':0.055555,
        'IP3R_init':0.119,
        'R_Gq_trimer_init':1.4107*0.27, 
        'R_L_Gq_trimer_init':0,
        'Gq_trimer_init': 0.61869,
        'R_L_kf':1.00,
        'R_L_kr':2.50,
        'R_Gq_trimer_kf':1.00,
        'R_Gq_trimer_kr':1.67,
        'R_L_Gq_trimer_kf':1.0,
        'R_L_Gq_trimer_kr': 0.0046,
        'Gq_trimerization_k':6.0,
        'Gq_trimer_split_k':0.0001,
        'R_L_Gq_trimer_split_k':0.04,
        'RGS4_Gq_a_GTP_kf':20.83,
        'RGS4_Gq_a_GTP_kr':33.32,
        'RGS4_Gq_a_GTP_diss_k':8.33,
        'Gq_a_GTP_decay_k':0.01,
        'Gq_a_GTP_PLCb_kf':2.52,
        'Gq_a_GTP_PLCb_kr':1.00,
        'Gq_a_GTP_PLCb_Ca_kf':30.0,
        'Gq_a_GTP_PLCb_Ca_kr':1.00,
        'Gq_a_GTP_PLCb_Ca_decay_k': 0.013,
        'PLCb_Ca_kf':3.00,
        'PLCb_Ca_kr':1.00,
        'PLCb_Ca_Gq_a_GTP_kf':25.2,
        'PLCb_Ca_Gq_a_GTP_kr':1.00,
        'DAG_deg_k':0.15,
        'IP3R_IP3_kf':50,
        'IPR3_IP3_kr':6.5,
        'IP3R_Ca_kf':20.0,
        'IPR3_Ca_kr':0.0806,
        'IP3R_IP3_Ca_kf':1.0,
        'IPR3_IP3_Ca_kr':0.5,
        'IP3R_Ca_IP3_kf':20.0,
        'IPR3_Ca_IP3_kr':14.5,
        'v1':800,
        'v8':0.15,
        'v4':0.5,
        'k4':0.09,
        'c2':0.185
    }



parameters={**defaultParameters}
Model()

##TIME OBSERVABLE
components_time_obs = create_t_obs()
time_obs = components_time_obs.t

##SPECIES
Monomer('R', ['R_b1', 'R_b2', 'R_s'], {'R_s':['inact', 'act']})
Monomer('L', ['L_b1'])
Monomer('Gq_a', ['Gq_a_b1', 'Gq_a_b2', 'Gq_a_s'], {'Gq_a_s' : ['GTP', 'GDP']})
Monomer('Gq_bg', ['Gq_bg_b1', 'Gq_bg_b2'])
Monomer('RGS4', ['RGS4_b1'])
Monomer('Ca', ['Ca_b1'])
Monomer('CaER')
Monomer('PLCb', ['PLCb_b1', 'PLCb_b2'])
Monomer('PIP2', ['PIP2_b1'])
Monomer('IP3', ['IP3_b1'])
Monomer('DAG', ['DAG_b1', 'DAG_s'], {'DAG_s' : ['act', 'inact']})
Monomer('IP3R', ['IP3R_b1', 'IP3R_b2'])

#INITIAL CONDITIONS
Initial(R(R_b1=None, R_b2=None, R_s='inact'), Parameter('R_init', parameters['R_init']))
Initial(L(L_b1=None), Parameter('L_init',  parameters['L_init']))
Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Parameter('Gq_a_GDP_init', parameters['Gq_a_GDP_init']))
Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('Gq_a_GTP_init', parameters['Gq_a_GTP_init']))
Initial(Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), Parameter('Gq_bg_init', parameters['Gq_bg_init']))
Initial(RGS4(RGS4_b1=None), Parameter('RGS4_init', parameters['RGS4_init']))
Initial(RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('RGS4_Gq_a_GTP_init', parameters['RGS4_Gq_a_GTP_init']))
Initial(Ca(Ca_b1=None), Parameter('Ca_init', parameters['Ca_init']))
Initial(CaER(), Parameter('CaER_init',parameters['CaER_init']))
Initial(PLCb(PLCb_b1=None, PLCb_b2=None), Parameter('PLCb_init', parameters['PLCb_init']))
Initial(PLCb(PLCb_b1=60, PLCb_b2=None)%Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP'), Parameter('PLCb_Gq_a_GTP_init', parameters['PLCb_Gq_a_GTP_init']))
Initial(PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70), Parameter('PLCb_Ca_init', parameters['PLCb_Ca_init']))
Initial(Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70), Parameter('PLCb_Ca_Gq_a_GTP_init', parameters['PLCb_Ca_Gq_a_GTP_init']))
Initial(PIP2(PIP2_b1=None), Parameter('PIP2_init', parameters['PIP2_init']))
Initial(IP3(IP3_b1=None), Parameter('IP3_init', parameters['IP3_init']))
Initial(DAG(DAG_b1=None, DAG_s='inact'), Parameter('DAG_init', parameters['DAG_init']))
Initial(IP3R(IP3R_b1=None, IP3R_b2=None), Parameter('IP3R_init', parameters['IP3R_init']))
Initial(R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Parameter('R_Gq_trimer_init', parameters['R_Gq_trimer_init']))
Initial(R(R_b1=30, R_b2=50, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50), Parameter('R_L_Gq_trimer_init', parameters['R_L_Gq_trimer_init']))
Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Parameter('Gq_trimer_init', parameters['Gq_trimer_init']))

#OBSERVABLES
Observable('obs_R', R(R_b1=None, R_b2=None, R_s='inact'))
Observable('obs_L', L(L_b1=None))
Observable('obs_Gq_a_GDP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'))
Observable('obs_Gq_a_GTP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'))
Observable('obs_Gq_bg', Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None)) 
Observable('obs_Ca', Ca(Ca_b1=None))
Observable('obs_CaER', CaER())
Observable('obs_PLCb', PLCb(PLCb_b1=None, PLCb_b2=None))
Observable('obs_PLCb_Ca_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70))
Observable('obs_PLCb_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None))
Observable('obs_PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70))
Observable('obs_PIP2', PIP2(PIP2_b1=None))    
Observable('obs_IP3', IP3(IP3_b1=None))
Observable('obs_DAG', DAG(DAG_b1=None, DAG_s='inact'))
Observable('obs_IP3R', IP3R(IP3R_b1=None, IP3R_b2=None))

################################################################################################################################################################################
##RULES
#Receptor-Ligand
Expression('R_L_kf', Piecewise((0, time_obs > parameters['time_out']),(parameters['R_L_kf'], time_obs > parameters['time_in']),(0, True)))
Parameter('R_L_kr', parameters['R_L_kr'])
Rule('R_L', R(R_b1=None, R_b2=None, R_s='inact') + L(L_b1=None) | R(R_b1=None, R_b2=50, R_s='act')%L(L_b1=50), R_L_kf, R_L_kr)
Observable('obs_RL', R(R_b1=None, R_b2=50, R_s='act')%L(L_b1=50))

#G-PROTEIN ACTIVATION
#R+Gtrimer
Parameter('R_Gq_trimer_kf', parameters['R_Gq_trimer_kf'])
Parameter('R_Gq_trimer_kr', parameters['R_Gq_trimer_kr'])
Rule('R_Gq_trimer', R(R_b1=None, R_b2=None, R_s='inact') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_Gq_trimer_kf, R_Gq_trimer_kr)
Observable('obs_RG', R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))

#RL+Gtrimer
Parameter('R_L_Gq_trimer_kf', parameters['R_L_Gq_trimer_kf'])
Parameter('R_L_Gq_trimer_kr',parameters['R_L_Gq_trimer_kr'])
Rule('R_L_Gq_trimer', R(R_b1=None, R_b2=50, R_s='act')%L(L_b1=50) + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_b2=50, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50), R_L_Gq_trimer_kf, R_L_Gq_trimer_kr)
Observable('obs_trimer', R(R_b1=30, R_b2=50, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50))

#Gq_trimerization
Parameter('Gq_trimerization_k', parameters['Gq_trimerization_k'])
Parameter('Gq_trimer_split_k', parameters['Gq_trimer_split_k'])
Rule('Gq_trimerization', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None) | Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimerization_k, Gq_trimer_split_k)

#RL_Gq split
Parameter('R_L_Gq_trimer_split_k', parameters['R_L_Gq_trimer_split_k'])
Rule('R_L_Gq_trimer_split', R(R_b1=30, R_b2=50, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50) >> R(R_b1=None, R_b2=50, R_s='act')%L(L_b1=50) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), R_L_Gq_trimer_split_k)

#Deactivation of Gq_a_GTP by RGS4
#RGS4 + Gq_a_GTP
Parameter('RGS4_Gq_a_GTP_kf', parameters['RGS4_Gq_a_GTP_kf'])
Parameter('RGS4_Gq_a_GTP_kr', parameters['RGS4_Gq_a_GTP_kr'])
Rule('RGS4_Gq_a_GTP', RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_kf, RGS4_Gq_a_GTP_kr)

#RGS4_Gq_a_GTP dissociation
Parameter('RGS4_Gq_a_GTP_diss_k', parameters['RGS4_Gq_a_GTP_diss_k'])
Rule('RGS4_Gq_a_GTP_diss', RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP') >> RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), RGS4_Gq_a_GTP_diss_k)

#Gq_a_GTP decay
Parameter('Gq_a_GTP_decay_k', parameters['Gq_a_GTP_decay_k'])
Rule('Gq_a_GTP_decay', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Gq_a_GTP_decay_k)

#ACTIVATION OF PLCb
#Gq_a_GTP + PLCb
Parameter('Gq_a_GTP_PLCb_kf', parameters['Gq_a_GTP_PLCb_kf'])
Parameter('Gq_a_GTP_PLCb_kr', parameters['Gq_a_GTP_PLCb_kr'])
Rule('Gq_a_GTP_PLCb', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + PLCb(PLCb_b1=None, PLCb_b2=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None), Gq_a_GTP_PLCb_kf, Gq_a_GTP_PLCb_kr)

#Gq_a_GTP_PLCb + Ca
Parameter('Gq_a_GTP_PLCb_Ca_kf', parameters['Gq_a_GTP_PLCb_Ca_kf'])
Parameter('Gq_a_GTP_PLCb_Ca_kr', parameters['Gq_a_GTP_PLCb_Ca_kr'])
Rule('Gq_a_GTP_PLCb_Ca', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None) + Ca(Ca_b1=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_kf, Gq_a_GTP_PLCb_Ca_kr)

#Gq_a_GTP_PLCb_Ca DECAY
Parameter('Gq_a_GTP_PLCb_Ca_decay_k',parameters['Gq_a_GTP_PLCb_Ca_decay_k'])
Rule('Gq_a_GTP_PLCb_Ca_diss', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_decay_k)

#PLCb + Ca
Parameter('PLCb_Ca_kf', parameters['PLCb_Ca_kf'])
Parameter('PLCb_Ca_kr', parameters['PLCb_Ca_kr'])
Rule('PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=None) + Ca(Ca_b1=None) | PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_kf, PLCb_Ca_kr)

#PLCb_Ca + Gq_a_GTP
Parameter('PLCb_Ca_Gq_a_GTP_kf', parameters['PLCb_Ca_Gq_a_GTP_kf'])
Parameter('PLCb_Ca_Gq_a_GTP_kr', parameters['PLCb_Ca_Gq_a_GTP_kr'])
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
Parameter('DAG_deg_k', parameters['DAG_deg_k'])
Rule('DAG_deg', DAG(DAG_b1=None, DAG_s='inact') >> None, DAG_deg_k)

#Kinetics of IP3R
#IP3R binds IP3
Parameter('IP3R_IP3_kf', parameters['IP3R_IP3_kf'])
Parameter('IPR3_IP3_kr', parameters['IPR3_IP3_kr'])
Rule('IP3R_IP3', IP3R(IP3R_b1=None, IP3R_b2=None) + IP3(IP3_b1=None) | IP3R(IP3R_b1=50, IP3R_b2=None)%IP3(IP3_b1=50), IP3R_IP3_kf, IPR3_IP3_kr)
Observable('obs_IP3R_IP3', IP3R(IP3R_b1=50, IP3R_b2=None)%IP3(IP3_b1=50))

#IP3R binds Ca
Parameter('IP3R_Ca_kf', parameters['IP3R_Ca_kf'])
Parameter('IPR3_Ca_kr', parameters['IPR3_Ca_kr'])
Rule('IP3R_Ca', IP3R(IP3R_b1=None, IP3R_b2=None) + Ca(Ca_b1=None) | IP3R(IP3R_b1=None, IP3R_b2=30)%Ca(Ca_b1=30), IP3R_Ca_kf, IPR3_Ca_kr)
Observable('obs_IP3R_Ca', IP3R(IP3R_b1=None, IP3R_b2=30)%Ca(Ca_b1=30))

#IP3R:IP3 binds Ca
Parameter('IP3R_IP3_Ca_kf', parameters['IP3R_IP3_Ca_kf'])
Parameter('IPR3_IP3_Ca_kr', parameters['IPR3_IP3_Ca_kr'])
Rule('IP3R_IP3_Ca', IP3R(IP3R_b1=50, IP3R_b2=None)%IP3(IP3_b1=50) + Ca(Ca_b1=None) | IP3R(IP3R_b1=50, IP3R_b2=30)%IP3(IP3_b1=50)%Ca(Ca_b1=30), IP3R_IP3_Ca_kf, IPR3_IP3_Ca_kr)

#IP3R:Ca binds IP3
Parameter('IP3R_Ca_IP3_kf', parameters['IP3R_Ca_IP3_kf'])
Parameter('IPR3_Ca_IP3_kr', parameters['IPR3_Ca_IP3_kr'])
Rule('IP3_Ca_IP3', IP3R(IP3R_b1=None, IP3R_b2=30)%Ca(Ca_b1=30) + IP3(IP3_b1=None) | IP3R(IP3R_b1=50, IP3R_b2=30)%IP3(IP3_b1=50)%Ca(Ca_b1=30), IP3R_Ca_IP3_kf, IPR3_Ca_IP3_kr)
Observable('obs_IP3R_IP3_Ca', IP3R(IP3R_b1=50, IP3R_b2=30)%IP3(IP3_b1=50)%Ca(Ca_b1=30))

#ER Ca release
Parameter('v1', parameters['v1'])
Parameter('v8', parameters['v8'])
Parameter('v4', parameters['v4'])
Parameter('k4', parameters['k4'])
Parameter('c2', parameters['c2'])
Expression('CaER_release_k', (c2*(v1*((obs_IP3R_IP3/IP3_init)**4)+v8)*(obs_CaER-obs_Ca)))
Rule('CaER_release', None >> Ca(Ca_b1=None), CaER_release_k)
Expression('Ca_inwards_k', (v4*(obs_Ca**2)/((obs_Ca**2)+(k4**2)))/obs_Ca)
Rule('Ca_inwards', Ca(Ca_b1=None) >> None, Ca_inwards_k)





