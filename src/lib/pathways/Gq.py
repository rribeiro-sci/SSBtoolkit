'''
                 THE OXYTOCIN RECEPTOR METABOLIC PATHWAY
                               VERSION 1.0
                       G alpha q11 coupled receptor
                   last modification 4 October 2020
  References:
   1. Chang, Chiung-wen, Ethan Poteet, John A. Schetz, Zeynep H. Gümüş,
      and Harel Weinstein. 2009. “Towards a Quantitative Representation
      of the Cell Signaling Mechanisms of Hallucinogens: Measurement and
      Mathematical Modeling of 5-HT1A and 5-HT2A Receptor-Mediated ERK1/2
      Activation.” Neuropharmacology 56 (Suppl 1): 213–25.
   2. Keizer, J, and G W De Young. 1992. “Two Roles of Ca2+ in Agonist
       Stimulated Ca2+ Oscillations.” Biophysical Journal 61 (3): 649–60.
'''

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


__author__ = "Rui Ribeiro"
__organizarion__ = "University of Verona"
__copyright__ = "Copyright 2020, Rui Ribeiro"
__credits__ = ["Rui Ribeiro","Pietro Micheli"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"
__status__ = "Production"

USAGE = __doc__.format(__author__, __email__)

def network(Rtotal, LR):

    Model()

    ##SPECIES
    #Receptor
    Monomer('R', ['R_b1', 'R_s'], {'R_s':['inact', 'act']})
    Parameter('R_init', LR)
    Initial(R(R_b1=None, R_s='act'), R_init)
    Observable('obs_RL', R(R_b1=None, R_s='act'))

    #Ligand
    #L_conc = 0.1
    #Monomer('L', ['L_b1'])
    #Parameter('L_init', L_conc)
    #Initial(L(L_b1=None), L_init)
    #Observable('obs_L', L(L_b1=None))

    #G-Protein
    Monomer('Gq_a', ['Gq_a_b1', 'Gq_a_b2', 'Gq_a_s'], {'Gq_a_s' : ['GTP', 'GDP']})
    Parameter('Gq_a_GDP_init', 0.0027739)
    Parameter('Gq_a_GTP_init', 6.4172E-4)
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Gq_a_GDP_init)
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'), Gq_a_GTP_init)
    Observable('obs_Gq_a_GDP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'))
    Observable('obs_Gq_a_GTP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'))

    Monomer('Gq_bg', ['Gq_bg_b1', 'Gq_bg_b2'])
    Parameter('Gq_bg_init', 0.0037173)
    Initial(Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), Gq_bg_init)
    Observable('obs_Gq_bg', Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None))

    #RGS4
    Monomer('RGS4', ['RGS4_b1'])
    Parameter('RGS4_init', 0.019994)
    Parameter('RGS4_Gq_a_GTP_init', 6.4168E-6)
    Initial(RGS4(RGS4_b1=None), RGS4_init)
    Initial(RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_init)

    #Ca2+
    Monomer('Ca', ['Ca_b1'])
    Parameter('Ca_init', 0.1)
    Initial(Ca(Ca_b1=None), Ca_init)
    Observable('obs_Ca', Ca(Ca_b1=None))

    #PLCb
    Monomer('PLCb', ['PLCb_b1', 'PLCb_b2'])
    Parameter('PLCb_init', 0.090022)
    Parameter('PLCb_Gq_a_GTP_init', 1.4492E-4)
    Parameter('PLCb_Ca_init', 0.0093825)
    Parameter('PLCb_Ca_Gq_a_GTP_init', 1.5038E-4)
    Initial(PLCb(PLCb_b1=None, PLCb_b2=None), PLCb_init)
    Initial(PLCb(PLCb_b1=60, PLCb_b2=None)%Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP'), PLCb_Gq_a_GTP_init)
    Initial(PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70), PLCb_Ca_init)
    Initial(Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70), PLCb_Ca_Gq_a_GTP_init)
    Observable('obs_PLCb', PLCb(PLCb_b1=None, PLCb_b2=None))
    Observable('obs_PLCb_Ca_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70))
    Observable('obs_PLCb_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None))
    Observable('obs_PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70))

    #PIP2
    Monomer('PIP2', ['PIP2_b1'])
    Parameter('PIP2_init', 2.6578)
    Initial(PIP2(PIP2_b1=None), PIP2_init)
    Observable('obs_PIP2', PIP2(PIP2_b1=None))

    #IP3
    Monomer('IP3', ['IP3_b1'])
    Parameter('IP3_init', 0.21952)
    Initial(IP3(IP3_b1=None), IP3_init)
    Observable('obs_IP3', IP3(IP3_b1=None))

    #DAG
    Monomer('DAG', ['DAG_b1', 'DAG_s'], {'DAG_s' : ['act', 'inact']})
    Parameter('DAG_init', 0.055555)
    Initial(DAG(DAG_b1=None, DAG_s='inact'), DAG_init)
    Observable('obs_DAG', DAG(DAG_b1=None, DAG_s='inact'))

    #

    ##Complexes
    Parameter('R_Gq_trimer_init', 0)
    Initial(R(R_b1=30,  R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_Gq_trimer_init)

    Parameter('R_L_Gq_trimer_init', 0)
    Initial(R(R_b1=30,  R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_L_Gq_trimer_init)

    Parameter('Gq_trimer_init', 0.61869)
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimer_init)

################################################################################################################################################################################
    ##RULES

    #G-PROTEIN ACTIVATION
    #R+Gtrimer
    #Parameter('R_Gq_trimer_kf', 1.00)
    #Parameter('R_Gq_trimer_kr', 1.67)
    #Rule('R_Gq_trimer', R(R_b1=None, R_b2=None, R_s='inact') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_Gq_trimer_kf, R_Gq_trimer_kr)
    #Observable('obs_RG', R(R_b1=30, R_b2=None, R_s='inact')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))

    #RL+Gtrimer
    Parameter('R_L_Gq_trimer_kf', 1.00)
    Parameter('R_L_Gq_trimer_kr', 0.0046)
    Rule('R_L_Gq_trimer', R(R_b1=None, R_s='act') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R_L_Gq_trimer_kf, R_L_Gq_trimer_kr)
    Observable('obs_trimer', R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))

    #Gq_trimerization
    Parameter('Gq_trimerization_k', 6.0)
    Parameter('Gq_trimer_split_k', 0.0001)
    Rule('Gq_trimerization', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None) | Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimerization_k, Gq_trimer_split_k)

    #RL_Gq split
    Parameter('R_L_Gq_trimer_split_k', 0.04)
    Rule('R_L_Gq_trimer_split', R(R_b1=30, R_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) >> R(R_b1=None, R_s='act') + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), R_L_Gq_trimer_split_k)

    #Deactivation of Gq_a_GTP by RGS4
    #RGS4 + Gq_a_GTP
    Parameter('RGS4_Gq_a_GTP_kf', 20.83)
    Parameter('RGS4_Gq_a_GTP_kr', 33.32)
    Rule('RGS4_Gq_a_GTP', RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_kf, RGS4_Gq_a_GTP_kr)

    #RGS4_Gq_a_GTP dissociation
    Parameter('RGS4_Gq_a_GTP_diss_k', 8.33)
    Rule('RGS4_Gq_a_GTP_diss', RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP') >> RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), RGS4_Gq_a_GTP_diss_k)

    #Gq_a_GTP decay
    Parameter('Gq_a_GTP_decay_k', 0.01)
    Rule('Gq_a_GTP_decay', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Gq_a_GTP_decay_k)

    #ACTIVATION OF PLCb
    #Gq_a_GTP + PLCb
    Parameter('Gq_a_GTP_PLCb_kf', 2.52)
    Parameter('Gq_a_GTP_PLCb_kr', 1.00)
    Rule('Gq_a_GTP_PLCb', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + PLCb(PLCb_b1=None, PLCb_b2=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None), Gq_a_GTP_PLCb_kf, Gq_a_GTP_PLCb_kr)

    #Gq_a_GTP_PLCb + Ca
    Parameter('Gq_a_GTP_PLCb_Ca_kf', 30.0)
    Parameter('Gq_a_GTP_PLCb_Ca_kr', 1.00)
    Rule('Gq_a_GTP_PLCb_Ca', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None) + Ca(Ca_b1=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_kf, Gq_a_GTP_PLCb_Ca_kr)

    #Gq_a_GTP_PLCb_Ca DECAY
    Parameter('Gq_a_GTP_PLCb_Ca_decay_k', 0.013)
    Rule('Gq_a_GTP_PLCb_Ca_diss', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_decay_k)

    #PLCb + Ca
    Parameter('PLCb_Ca_kf', 3.00)
    Parameter('PLCb_Ca_kr', 1.00)
    Rule('PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=None) + Ca(Ca_b1=None) | PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_kf, PLCb_Ca_kr)

    #PLCb_Ca + Gq_a_GTP
    Parameter('PLCb_Ca_Gq_a_GTP_kf', 25.2)
    Parameter('PLCb_Ca_Gq_a_GTP_kr', 1.00)
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
    Parameter('DAG_deg_k', 0.15)
    Rule('DAG_deg', DAG(DAG_b1=None, DAG_s='inact') >> None, DAG_deg_k)

    return model

list_of_observables=['obs_RL','obs_IP3']
