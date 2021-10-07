from pysb import *
from pysb.macros import *
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator
import math

def network(pKd, ligand_conc,epsilon):

    Model()

    #total receptor
    #Rtotal = Rtotal
    Rtotal = 2 #uM
    #from pKd to Kd
    Kd = (10**(-(float(pKd)))) * 10**6
    #LR determination
    a = 1
    b = float(Rtotal)+float(ligand_conc)+Kd
    c = float(Rtotal)*float(ligand_conc)
    delta = (b**2) - (4*a*c)
    RL = (b-math.sqrt(delta))/(2*a)


# MONOMERS

    Monomer('R5HT2A', ['R5HT2A_b1', 'R5HT2A_b2', 'R5HT2A_s', 'R5HT2A_l'], {'R5HT2A_s':['inact', 'act'], 'R5HT2A_l' : ['intra', 'membrane']})
    Monomer('L', ['L_b1'])
    Monomer('Gq_a', ['Gq_a_b1', 'Gq_a_b2', 'Gq_a_s'], {'Gq_a_s' : ['GTP', 'GDP']})
    Monomer('Gq_bg', ['Gq_bg_b1', 'Gq_bg_b2'])
    Monomer('RGS4', ['RGS4_b1'])
    Monomer('PLCb', ['PLCb_b1', 'PLCb_b2'])
    Monomer('Ca', ['Ca_b1'])
    Monomer('Pi')
    Monomer('PI')
    Monomer('PIP2', ['PIP2_b1'])
    Monomer('IP3')
    Monomer('DAG', ['DAG_b1', 'DAG_s'], {'DAG_s' : ['act', 'inact']})
    Monomer('PA')
    Monomer('PKC', ['PKC_b1', 'PKC_b2'])
    Monomer('PP2A', ['PP2A_b1'])
    Monomer('PLA2', ['PLA2_b1', 'PLA2_b2', 'PLA2_s', 'PIP2_b'], {'PLA2_s' : ['act', 'inact'], 'PIP2_b' : ['True', 'False']})
    Monomer('APC')
    Monomer('AA', ['AA_b1', 'AA_s'], {'AA_s' : ['act', 'inact']})
    Monomer('Raf', ['Raf_b1', 'Raf_s'], {'Raf_s' : ['act', 'act2', 'inact']})
    Monomer('MEK', ['MEK_b1', 'MEK_p'], {'MEK_p' : ['p0', 'p1', 'p2']})
    Monomer('ERK', ['ERK_b1', 'ERK_p'], {'ERK_p' : ['p0', 'p1', 'p2']})
    Monomer('MKP', ['MKP_b1', 'MKP_p', 'MKP_s'], {'MKP_p' : ['p0', 'p2'], 'MKP_s' : ['act', 'inact']})
    Monomer('Ras', ['Ras_b1', 'Ras_s'], {'Ras_s' : ['GTP', 'GDP']})
    Monomer('RasGEF', ['RasGEF_s'], {'RasGEF_s' : ['act', 'inact']})
    Monomer('RasGAP', ['RasGAP_s'], {'RasGAP_s' : ['act', 'inact']})

# PARAMETERS

    # Initial Concentrations (mol)

    # Original R, L model parameters
    # Parameter('R5HT2A_init', 1.0192)
    # Parameter('R5HT2A_L_init', 0)
    # Parameter('R5HT2A_Gq_trimer_init', 0.37759)
    # Parameter('R5HT2A_L_Gq_trimer_init', 0)
    # Parameter('L_init', 0.01)
    # Parameter('R5HT2A_intra_init', 0.010192)

    # User R, L parameters
    Parameter('R5HT2A_init', Rtotal * 0.73)
    Parameter('R5HT2A_L_init', RL * 0.73)
    Parameter('R5HT2A_Gq_trimer_init', Rtotal*0.27)
    Parameter('R5HT2A_L_Gq_trimer_init', RL*0.27)
    Parameter('R5HT2A_intra_init',  Rtotal * 0.01)


    Parameter('Gq_a_GDP_init', 0.0027739)
    Parameter('Gq_a_GTP_init', 6.4172E-4)
    Parameter('Gq_bg_init', 0.0037173)
    Parameter('Gq_trimer_init', 0.61869)
    Parameter('RGS4_init', 0.019994)
    Parameter('RGS4_Gq_a_GTP_init', 6.4168E-6)
    Parameter('PLCb_init', 0.090022)
    Parameter('PLCb_Gq_a_GTP_init', 1.4492E-4)
    Parameter('PLCb_Ca_init', 0.0093825)
    Parameter('PLCb_Ca_Gq_a_GTP_init', 1.5038E-4)
    Parameter('PIP2_init', 2.6578)
    Parameter('IP3_init', 0.21952)
    Parameter('DAG_init', 0.055555)
    Parameter('Pi_init', 618.69)
    Parameter('Ca_init', 0.034739)
    Parameter('PA_init', 83332.0)
    Parameter('PKC_init', 0.062951)
    Parameter('PKC_Ca_init', 0.0026242)
    Parameter('PKC_DAG_init', 0.58385)
    Parameter('PKC_Ca_DAG_init', 1.3507E-7)
    Parameter('PKC_Ca_DAG_act_init', 1.3507E-6)
    Parameter('PKC_DAG_AA_init', 0.031811)
    Parameter('PKC_DAG_AA_act_init', 0.31811)
    Parameter('PKC_AA_init', 4.5732E-4)
    Parameter('PKC_Ca_AA_init', 1.9064E-4)
    #
    Parameter('APC_init',  23.496)
    Parameter('AA_init', 6.0539)
    #
    Parameter('PLA2_init', 3.5823E-4)
    Parameter('PLA2_a_init', 1.489E-5)
    Parameter('PLA2_a_ERK_PP_init', 2.5313E-7)
    Parameter('PLA2_Ca_PIP2_init', 0.0027563)
    Parameter('PLA2_PIP2_init', 0.39671)
    Parameter('PLA2_Ca_init', 1.2445E-4)
    Parameter('PLA2_a_Ca_init', 3.1036E-5)
    #
    Parameter('Ras_GEF_init', 0.072743)
    Parameter('Ras_GEF_a_init', 0.027257)
    Parameter('Ras_GAP_a_init', 0.0015858)
    #
    Parameter('ERK_PP_MKP_PP_init', 0.010263)
    Parameter('ERK_P_MKP_PP_init', 0.003791)
    Parameter('ERK_P_MKP_PP_act_init', 0.0018979)
    Parameter('Raf_aa_init', 0.0021849)
    Parameter('MEK_PP_PP2A_Init', 3.0847E-4)
    Parameter('Raf_aa_PP2A_init', 2.2131E-5)
    Parameter('ERK_PP_Raf_a_init', 2.2131E-5)
    Parameter('Raf_a_PP2A_init', 5.8373E-4)
    Parameter('ERK_P_MEK_PP_init', 6.2945E-5)
    Parameter('ERK_MEK_PP_init', 0.094895)
    Parameter('MEK_PP_init', 0.030392)
    Parameter('ERK_MKP_PP_init', 0.066479)
    Parameter('MKP_init', 0.48458)
    Parameter('ERK_PP_MKP_init', 0.0054982)
    Parameter('MKP_PP_init', 0.027491)
    Parameter('ERK_PP_init', 0.0090591)
    Parameter('ERK_P_init', 0.010356)
    Parameter('ERK_init', 0.15768)
    Parameter('MEK_P_Raf_a_init', 0.017627)
    Parameter('MEK_P_PP2A_init', 4.9691E-4)
    Parameter('PP2A_init', 0.15859)
    Parameter('MEK_P_init', 0.048958)
    Parameter('MEK_Raf_a_init', 0.028395)
    Parameter('MEK_init', 0.078866)
    Parameter('Raf_a_init', 0.057628)
    Parameter('Ras_GTP_Raf_init', 0.023412)
    Parameter('Raf_init', 0.030124)
    Parameter('Ras_GTP_init', 0.02024)
    Parameter('Ras_GAP_init', 4.1421E-4)
    Parameter('Ras_GDP_init', 0.076348)

    # Rules kinetik parameters (1/mol*s for kf and for k, 1/s for kr)

    Parameter('R5HT2A_ext_kf', 0.005)
    Parameter('R5HT2A_ext_kr', 5.0E-5)
    Parameter('R_L_kf', 1.00)
    Parameter('R_L_kr', 2.50)
    Parameter('RG_L_kf', 1.00)
    Parameter('RG_L_kr', 0.007)
    Parameter('RL_inter_k', 0.00235)
    Parameter('RL_intra_deg_k', 6.67E-4)
    Parameter('L_deg_k', 2.43E-5)

    Parameter('R5HT2A_Gq_trimer_kf', 1.00)
    Parameter('R5HT2A_Gq_trimer_kr', 1.67)
    Parameter('R5HT2A_L_Gq_trimer_kf', 1.00)
    Parameter('R5HT2A_L_Gq_trimer_kr', 0.0046)
    Parameter('Gq_trimer_split_k', 0.0001)
    Parameter('Gq_trimerization_k', 6.0)
    Parameter('R5HT2A_L_Gq_trimer_split_k', 0.04)
    Parameter('RGS4_Gq_a_GTP_kf', 20.83)
    Parameter('RGS4_Gq_a_GTP_kr', 33.32)
    Parameter('RGS4_Gq_a_GTP_diss_k', 8.33)
    Parameter('Gq_a_GTP_decay_k', 0.01)
    Parameter('Gq_a_GTP_PLCb_kf', 2.52)
    Parameter('Gq_a_GTP_PLCb_kr', 1.00)
    Parameter('Gq_a_GTP_PLCb_Ca_kf', 30.0)
    Parameter('Gq_a_GTP_PLCb_Ca_kr', 1.00)
    Parameter('Pi_gen_k', 1.00)
    Parameter('PLCb_Ca_kf', 3.00)
    Parameter('PLCb_Ca_kr', 1.00)
    Parameter('PLCb_Ca_Gq_a_GTP_kf', 25.2)
    Parameter('PLCb_Ca_Gq_a_GTP_kr', 1.00)
    Parameter('Gq_a_GTP_PLCb_Ca_diss_k', 0.013)
    Observable('obs_PIP2', PIP2(PIP2_b1=None))
    Expression('PIP2_PLCb_Ca_k', (10.0 * obs_PIP2)/(40.13 + obs_PIP2))
    Expression('PIP2_Gq_a_GTP_PLCb_Ca_k', (48.0 * obs_PIP2)/(5.00 + obs_PIP2))
    Observable('obs_IP3', IP3())
    Expression('IP3_deg_k', 0.14 - ((0.14 * 0.16)/obs_IP3))
    Parameter('DAG_deg_k', 0.15)
    Parameter('PLA2_Ca_kf', 1.00)
    Parameter('PLA2_Ca_kr', 0.10)
    Expression('PLA2_Ca_PIP2_kf', 0.10 * obs_PIP2)
    Parameter('PLA2_Ca_PIP2_kr', 0.012)
    Expression('PLA2_PIP2_kf', 0.50 * obs_PIP2)
    Parameter('PLA2_PIP2_kr', 0.0012)
    Parameter('PLA2_a_ERK_PP_kf', 3.9)
    Parameter('PLA2_a_ERK_PP_kr', 40.0)
    Parameter('PLA2_a_release_k', 10.0)
    Parameter('PLA2_a_Ca_kf', 6.00)
    Parameter('PLA2_a_Ca_kr', 0.10)
    Parameter('PLA2_a_deact_k', 0.17)
    Observable('obs_APC', APC())
    Expression('AA_from_PLA2_a_Ca_k', 120.0/(20.0 + obs_APC))
    Expression('AA_from_PLA2_Ca_k', 5.4/(20.0 + obs_APC))
    Expression('AA_from_PLA2_Ca_PIP2_k', 36.0/(20.0 + obs_APC))
    Expression('AA_from_PLA2_PIP2_k', 11.04/(20.0 + obs_APC))


    Parameter('AA_deg_k', 0.40)
    Parameter('PKC_Ca_kf', 0.60)
    Parameter('PKC_Ca_kr', 0.50)
    Parameter('PKC_Ca_DAG_kf', 0.008)
    Parameter('PKC_Ca_DAG_kr', 8.635)
    Parameter('PKC_Ca_DAG_act_kf', 1.00)
    Parameter('PKC_Ca_DAG_act_kr', 0.10)
    Parameter('PKC_DAG_kf', 0.10)
    Parameter('PKC_DAG_kr', 5.99E-4)
    Parameter('PKC_DAG_AA_kf', 0.018)
    Parameter('PKC_DAG_AA_kr', 2.00)
    Parameter('PKC_DAG_AA_act_kf', 2.00)
    Parameter('PKC_DAG_AA_act_kr', 0.20)
    Parameter('PKC_AA_kf', 1.2E-4)
    Parameter('PKC_AA_kr', 0.100)
    Parameter('PKC_Ca_AA_kf', 0.001)
    Parameter('PKC_Ca_AA_kr', 0.100)
    Observable('obs_Raf_i', Raf(Raf_b1=None, Raf_s='inact'))
    Expression('Raf_a_1_k', 4.00/(66.67 + obs_Raf_i))
    Expression('Raf_a_2_k', 4.00/(66.67 + obs_Raf_i))
    Expression('Raf_a_3_k', 4.00/(66.67 + obs_Raf_i))
    Expression('Raf_a_4_k', 4.00/(66.67 + obs_Raf_i))
    Parameter('MEK_Raf_a_kf', 3.28)
    Parameter('MEK_Raf_a_kr', 0.42)
    Parameter('MEK_phospho_k', 0.105)
    Parameter('MEK_P_PP2A_kf', 1.92)
    Parameter('MEK_P_PP2A_kr', 24.00)
    Parameter('MEK_P_dephospho_k', 6.00)
    Parameter('MEK_P_Raf_a_kf', 3.28)
    Parameter('MEK_P_Raf_a_kr', 0.42)
    Parameter('MEK_P_phospho_k', 0.105)
    Parameter('Raf_a_PP2A_kf', 1.916)
    Parameter('Raf_a_PP2A_kr', 24.00)
    Parameter('Raf_a_deact_k', 6.00)
    Parameter('ERK_PP_Raf_a_kf', 1.95)
    Parameter('ERK_PP_Raf_a_kr', 40.00)
    Parameter('Raf_a_act_k', 6.00)
    Parameter('Raf_aa_PP2A_kf', 1.916)
    Parameter('Raf_aa_PP2A_kr', 24.00)
    Parameter('Raf_aa_deact_k', 6.0)
    Parameter('ERK_MEK_PP_kf', 20.0)
    Parameter('ERK_MEK_PP_kr', 1.00)
    Parameter('ERK_phospho_k', 0.01)
    Parameter('ERK_P_MEK_PP_kf', 3.20)
    Parameter('ERK_P_MEK_PP_kr', 1.00)
    Parameter('ERK_P_phospho_k', 15.00)
    Parameter('ERK_PP_MKP_kf', 31.563)
    Parameter('ERK_PP_MKP_kr', 20.20)
    Parameter('MKP_phospho_k', 5.00)
    Parameter('ERK_PP_MKP_PP_kf', 45.00)
    Parameter('ERK_PP_MKP_PP_kr', 1.00)
    Parameter('ERK_PP_dephospho_k', 0.092)
    Parameter('ERK_P_MKP_PP_kf', 1.00)
    Parameter('ERK_P_MKP_PP_kr', 10.00)
    Parameter('MKP_PP_deg_k', 1.00)
    Parameter('MKP_PP_act_kf', 10.00)
    Parameter('MKP_PP_act_kr', 1.00)
    Parameter('ERK_P_dephospho_k', 0.50)
    Parameter('ERK_MKP_PP_kf', 0.086)
    Parameter('ERK_MKP_PP_kr', 1.100)
    Parameter('MEK_PP_PP2A_kf', 1.92)
    Parameter('MEK_PP_PP2A_kr', 24.00)
    Parameter('MEK_PP_dephospho_k', 6.00)
    Observable('obs_RasGEF_i', RasGEF(RasGEF_s='inact'))
    Expression('RasGEF_a_1_k', 4.00/(3.33 + obs_RasGEF_i))
    Expression('RasGEF_a_2_k', 4.00/(3.33 + obs_RasGEF_i))
    Expression('RasGEF_a_3_k', 4.00/(3.33 + obs_RasGEF_i))
    Expression('RasGEF_a_4_k', 4.00/(3.33 + obs_RasGEF_i))
    Parameter('RasGEF_deact_k', 1.00)
    Observable('obs_RasGAP_i', RasGAP(RasGAP_s='inact'))
    Expression('RasGAP_a_1_k', 4.00/(3.33 + obs_RasGAP_i))
    Expression('RasGAP_a_2_k', 4.00/(3.33 + obs_RasGAP_i))
    Expression('RasGAP_a_3_k', 4.00/(3.33 + obs_RasGAP_i))
    Expression('RasGAP_a_4_k', 4.00/(3.33 + obs_RasGAP_i))
    Parameter('RasGAP_deact_k', 0.1)
    Observable('obs_Ras_GDP', Ras(Ras_b1=None, Ras_s='GDP'))
    Expression('Ras_GTP_act_k', 0.020/(0.505 + obs_Ras_GDP))
    Observable('obs_Ras_GTP', Ras(Ras_b1=None, Ras_s='GTP'))
    Expression('Ras_GTP_inact_k', 10.00/(1.01 + obs_Ras_GTP))
    Parameter('Ras_GDP_Ras_GTP_kf', 1.6E-4)
    Parameter('Ras_GDP_Ras_GTP_kr', 1.2E-4)
    Parameter('Ras_GTP_Raf_kf', 24.00)
    Parameter('Ras_GTP_Raf_kr', 0.5)
    Parameter('Raf_act_k', 0.125)

# INITIAL CONDITIONS

    # Initial(L(L_b1=None), L_init)
    Initial(R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='intra'), R5HT2A_intra_init)
    Initial(R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane'), R5HT2A_init)
    Initial(R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%L(L_b1=50), R5HT2A_L_init)
    Initial(R5HT2A(R5HT2A_b1=30, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R5HT2A_Gq_trimer_init)
    Initial(R5HT2A(R5HT2A_b1=30, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50), R5HT2A_L_Gq_trimer_init)
    #
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'), Gq_a_GDP_init)
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'), Gq_a_GTP_init)
    Initial(Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), Gq_bg_init)
    Initial(Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimer_init)
    #
    Initial(RGS4(RGS4_b1=None), RGS4_init)
    Initial(RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_init)
    #
    Initial(Ca(Ca_b1=None), Ca_init)
    #
    Initial(PLCb(PLCb_b1=None, PLCb_b2=None), PLCb_init)
    Initial(PLCb(PLCb_b1=60, PLCb_b2=None)%Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP'), PLCb_Gq_a_GTP_init)
    Initial(PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70), PLCb_Ca_init)
    Initial(Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70), PLCb_Ca_Gq_a_GTP_init)
    #
    Initial(PIP2(PIP2_b1=None), PIP2_init)
    Initial(Pi(), Pi_init)
    Initial(PA(), PA_init)
    Initial(IP3(), IP3_init)
    Initial(DAG(DAG_b1=None, DAG_s='inact'), DAG_init)
    #
    Initial(PKC(PKC_b1=None, PKC_b2=None), PKC_init)
    Initial(PKC(PKC_b1=30, PKC_b2=None)%Ca(Ca_b1=30), PKC_Ca_init)
    Initial(PKC(PKC_b1=None, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact'), PKC_DAG_init)
    Initial(PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='inact'), PKC_Ca_DAG_init)
    Initial(PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'), PKC_Ca_DAG_act_init)
    Initial(PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='inact'), PKC_DAG_AA_init)
    Initial(PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='act'), PKC_DAG_AA_act_init)
    Initial(PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'), PKC_AA_init)
    Initial(PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'), PKC_Ca_AA_init)
    #
    Initial(APC(), APC_init)
    Initial(AA(AA_b1=None, AA_s='inact'), AA_init)
    #
    Initial(PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='False'), PLA2_init)
    Initial(PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='act', PIP2_b='False'), PLA2_a_init)
    Initial(PLA2(PLA2_b1=None, PLA2_b2=40,  PLA2_s='act', PIP2_b='False')%ERK(ERK_b1=40, ERK_p='p2'), PLA2_a_ERK_PP_init)
    Initial(PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='inact', PIP2_b='True')%Ca(Ca_b1=30), PLA2_Ca_PIP2_init)
    Initial(PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='True'), PLA2_PIP2_init)
    Initial(PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='inact', PIP2_b='False')%Ca(Ca_b1=30), PLA2_Ca_init)
    Initial(PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=30), PLA2_a_Ca_init)
    #
    Initial(RasGEF(RasGEF_s='inact'), Ras_GEF_init)
    Initial(RasGEF(RasGEF_s='act'), Ras_GEF_a_init)
    Initial(RasGAP(RasGAP_s='inact'), Ras_GAP_init)
    Initial(RasGAP(RasGAP_s='act'), Ras_GAP_a_init)
    #
    Initial(ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_PP_MKP_PP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_P_MKP_PP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='act'), ERK_P_MKP_PP_act_init)
    Initial(Raf(Raf_b1=None, Raf_s='act2'), Raf_aa_init)
    Initial(MEK(MEK_b1=30, MEK_p='p2')%PP2A(PP2A_b1=30), MEK_PP_PP2A_Init)
    Initial(Raf(Raf_b1=30, Raf_s='act2')%PP2A(PP2A_b1=30), Raf_aa_PP2A_init)
    Initial(ERK(ERK_b1=30, ERK_p='p2')%Raf(Raf_b1=30, Raf_s='act'), ERK_PP_Raf_a_init)
    Initial(Raf(Raf_b1=30, Raf_s='act')%PP2A(PP2A_b1=30), Raf_a_PP2A_init)
    Initial(ERK(ERK_b1=None, ERK_p='p2'), ERK_PP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p1')%MEK(MEK_b1=30, MEK_p='p2'), ERK_P_MEK_PP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p0')%MEK(MEK_b1=30, MEK_p='p2'), ERK_MEK_PP_init)
    Initial(MEK(MEK_b1=None, MEK_p='p2'), MEK_PP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p0')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_MKP_PP_init)
    Initial(MKP(MKP_b1=None, MKP_p='p0', MKP_s='inact'), MKP_init)
    Initial(ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p0', MKP_s='inact'), ERK_PP_MKP_init)
    Initial(MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact'), MKP_PP_init)
    Initial(ERK(ERK_b1=None, ERK_p='p1'), ERK_P_init)
    Initial(ERK(ERK_b1=None, ERK_p='p0'), ERK_init)
    Initial(MEK(MEK_b1=30, MEK_p='p1')%Raf(Raf_b1=30, Raf_s='act'), MEK_P_Raf_a_init)
    Initial(MEK(MEK_b1=30, MEK_p='p1')%PP2A(PP2A_b1=30), MEK_P_PP2A_init)
    Initial(PP2A(PP2A_b1=None), PP2A_init)
    Initial(MEK(MEK_b1=None, MEK_p='p1'), MEK_P_init)
    Initial(MEK(MEK_b1=30, MEK_p='p0')%Raf(Raf_b1=30, Raf_s='act'), MEK_Raf_a_init)
    Initial(MEK(MEK_b1=None, MEK_p='p0'), MEK_init)
    Initial(Raf(Raf_b1=None, Raf_s='inact'), Raf_init)
    Initial(Raf(Raf_b1=None, Raf_s='act'), Raf_a_init)
    Initial(Ras(Ras_b1=70, Ras_s='GTP')%Raf(Raf_b1=70, Raf_s='inact'), Ras_GTP_Raf_init)
    Initial(Ras(Ras_b1=None, Ras_s='GTP'), Ras_GTP_init)
    Initial(Ras(Ras_b1=None, Ras_s='GDP'), Ras_GDP_init)

# REACTIONS

# Receptor activation

    # Rule('R_externaliz', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='intra') | R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane'), R5HT2A_ext_kf, R5HT2A_ext_kr)
    # Rule('R_L', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane') + L(L_b1=None) | R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%L(L_b1=50), R_L_kf, R_L_kr)
    # Rule('RG_L', R5HT2A(R5HT2A_b1=30, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) + L(L_b1=None) | R5HT2A(R5HT2A_b1=30, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50), RG_L_kf, RG_L_kr)
    # Rule('L_deg', L(L_b1=None) >> None, L_deg_k)
    Rule('RL_internaliz', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%L(L_b1=50) >> R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='intra')%L(L_b1=50), RL_inter_k)
    Rule('RL_intra_deg', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='intra')%L(L_b1=50) >> None, RL_intra_deg_k)

# IP3, DAG production (16 Rules)

    # R+Gtrimer
    Rule('reaction1', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane') + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R5HT2A(R5HT2A_b1=30, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), R5HT2A_Gq_trimer_kf, R5HT2A_Gq_trimer_kr)
    # RL+Gtrimer
    Rule('reaction2', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%L(L_b1=50) + Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) | R5HT2A(R5HT2A_b1=30, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50), R5HT2A_L_Gq_trimer_kf, R5HT2A_L_Gq_trimer_kr)
    # Gq_trimer split
    Rule('reaction3', Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None) >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), Gq_trimer_split_k)
    # Gq_trimerization
    Rule('reaction4', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None) >> Gq_a(Gq_a_b1=None, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None), Gq_trimerization_k)
    # RL_Gq split
    Rule('reaction5', R5HT2A(R5HT2A_b1=30, R5HT2A_b2=50, R5HT2A_s='act')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50) >> R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act')%L(L_b1=50) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None), R5HT2A_L_Gq_trimer_split_k)
    # RGS4 + Gq_a_GTP
    Rule('reaction6', RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'), RGS4_Gq_a_GTP_kf, RGS4_Gq_a_GTP_kr)
    # RGS4_Gq_a_GTP dissociation
    Rule('reaction7', RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP') >> RGS4(RGS4_b1=None) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Pi(), RGS4_Gq_a_GTP_diss_k)
    # Gq_a_GTP decay
    Rule('reaction8', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + Pi(), Gq_a_GTP_decay_k)
    # Gq_a_GTP + PLCb
    Rule('reaction9', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') + PLCb(PLCb_b1=None, PLCb_b2=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None), Gq_a_GTP_PLCb_kf, Gq_a_GTP_PLCb_kr)
    # Gq_a_GTP_PLCb + Ca
    Rule('reaction10', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None) + Ca(Ca_b1=None) | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), Gq_a_GTP_PLCb_Ca_kf, Gq_a_GTP_PLCb_Ca_kr)
    # PLCb + Ca
    Rule('reaction11', PLCb(PLCb_b1=None, PLCb_b2=None) + Ca(Ca_b1=None) | PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_kf, PLCb_Ca_kr)
    # PLCb_Ca + Gq_a_GTP
    Rule('reaction12', PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1) + Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP') | Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), PLCb_Ca_Gq_a_GTP_kf, PLCb_Ca_Gq_a_GTP_kr)
    # Gq_a_GTP_PLCb_Ca dissosiation
    Rule('reaction13', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1) + Pi(), Gq_a_GTP_PLCb_Ca_diss_k)
    # PIP2_PLCb_Ca
    Rule('reaction14', PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1) >> IP3() + DAG(DAG_b1=None, DAG_s='inact') + PLCb(PLCb_b1=None, PLCb_b2=1)%Ca(Ca_b1=1), PIP2_PLCb_Ca_k )
    # PIP2_Gq_a_GTP_PLCb_Ca
    Rule('reaction15', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1) >> IP3() + DAG(DAG_b1=None, DAG_s='inact') + Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=1)%Ca(Ca_b1=1), PIP2_Gq_a_GTP_PLCb_Ca_k )
    # Pi generation
    Rule('reaction16', PA() >> PI(), Pi_gen_k)

# PLA2 activation (7 rules)

    # PLA2 + Ca
    Rule('PLA2_Ca', PLA2(PLA2_b1=None, PLA2_b2=None,  PLA2_s='inact', PIP2_b='False') + Ca(Ca_b1=None) | PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='inact', PIP2_b='False')%Ca(Ca_b1=1), PLA2_Ca_kf, PLA2_Ca_kr)
    # PLA2_Ca + PIP2
    Rule('PLA2_Ca_PIP2', PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='inact', PIP2_b='False')%Ca(Ca_b1=1)  | PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='inact', PIP2_b='True')%Ca(Ca_b1=1), PLA2_Ca_PIP2_kf, PLA2_Ca_PIP2_kr)
    # PLA2 + PIP2
    Rule('PLA2_PIP2', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='False') | PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='True'), PLA2_PIP2_kf, PLA2_PIP2_kr)
    # PLA2 + ERK_PP
    Rule('PLA2_ERK_PP', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='False') + ERK(ERK_b1=None, ERK_p='p2') | PLA2(PLA2_b1=None, PLA2_b2=40,  PLA2_s='act', PIP2_b='False')%ERK(ERK_b1=40, ERK_p='p2'), PLA2_a_ERK_PP_kf, PLA2_a_ERK_PP_kr)
    # active PLA2 relese
    Rule('PLA2_a_release', PLA2(PLA2_b1=None, PLA2_b2=40, PLA2_s='act', PIP2_b='False')%ERK(ERK_b1=40, ERK_p='p2') >> PLA2(PLA2_b1=None, PLA2_b2=None,  PLA2_s='act', PIP2_b='False') + ERK(ERK_b1=None, ERK_p='p2'), PLA2_a_release_k)
    # PLA2_act + Ca
    Rule('PLA2_a_Ca', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='act', PIP2_b='False') + Ca(Ca_b1=None) | PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=1), PLA2_a_Ca_kf, PLA2_a_Ca_kr)
    # PLA2_act deactivation
    Rule('PLA2_a_deact', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='act', PIP2_b='False') >> PLA2(PLA2_b1=None, PLA2_b2=None,  PLA2_s='inact', PIP2_b='False'), PLA2_a_deact_k)

# AA production

    # from APC to AA by PLA2_Ca
    Rule('AA_from_PLA2_Ca', APC() + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=1) >> AA(AA_b1=None, AA_s='inact') + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=1), AA_from_PLA2_Ca_k)
    # from APC to AA by PLA2_Ca_PIP2
    Rule('AA_from_PLA2_Ca_PIP2', APC() + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='inact', PIP2_b='True')%Ca(Ca_b1=1) >> AA(AA_b1=None, AA_s='inact') + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='inact', PIP2_b='True')%Ca(Ca_b1=1), AA_from_PLA2_Ca_PIP2_k)
    # from APC to AA by PLA2_PIP2
    Rule('AA_from_PLA2_PIP2', APC() + PLA2(PLA2_b1=None, PLA2_b2=None,  PLA2_s='inact', PIP2_b='True') >> AA(AA_b1=None, AA_s='inact') + PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='True') , AA_from_PLA2_PIP2_k)
    # from APC to AA by PLA2_a_Ca
    Rule('AA_from_PLA2_a_Ca', APC() + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=1) >> AA(AA_b1=None, AA_s='inact') + PLA2(PLA2_b1=1, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=1), AA_from_PLA2_a_Ca_k)
    # AA degrade
    Rule('AA_deg', AA(AA_b1=None, AA_s='inact') >> APC(), AA_deg_k)

# PKC activation (10 rules)

    #IP3 degradation
    Rule('IP3_deg', IP3() >> None, IP3_deg_k)
    #DAG degradation
    Rule('DAG_deg', DAG(DAG_b1=None, DAG_s='inact') >> PA(), DAG_deg_k)
    # PKC + Ca
    Rule('PKC_Ca', PKC(PKC_b1=None, PKC_b2=None) + Ca(Ca_b1=None) | PKC(PKC_b1=30, PKC_b2=None)%Ca(Ca_b1=30), PKC_Ca_kf, PKC_Ca_kr)
    # PKC_Ca + DAG
    Rule('PKC_Ca_DAG', PKC(PKC_b1=30, PKC_b2=None)%Ca(Ca_b1=30) + DAG(DAG_b1=None, DAG_s='inact') | PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='inact'), PKC_Ca_DAG_kf, PKC_Ca_DAG_kr)
    # PKC_Ca_DAG activation
    Rule('PKC_Ca_DAG_act', PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='inact') | PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'), PKC_Ca_DAG_act_kf, PKC_Ca_DAG_act_kr)
    # PKC + DAG
    Rule('PKC_DAG', PKC(PKC_b1=None, PKC_b2=None) + DAG(DAG_b1=None, DAG_s='inact') | PKC(PKC_b1=None, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact'), PKC_DAG_kf, PKC_DAG_kr)
    # AA + PKC_DAG
    Rule('PKC_DAG_AA', AA(AA_b1=None, AA_s='inact') + PKC(PKC_b1=None, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact') | PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='inact'), PKC_DAG_AA_kf, PKC_DAG_AA_kr)
    # PKC_DAG_AA activation
    Rule('PKC_DAG_AA_act', PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='inact') | PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='act'), PKC_DAG_AA_act_kf, PKC_DAG_AA_act_kr)
    # PKC + AA
    Rule('PKC_AA', PKC(PKC_b1=None, PKC_b2=None) + AA(AA_b1=None, AA_s='inact') | PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'), PKC_AA_kf, PKC_AA_kr )
    # PKC_Ca + AA
    Rule('PKC_Ca_AA', PKC(PKC_b1=30, PKC_b2=None)%Ca(Ca_b1=30) + AA(AA_b1=None, AA_s='inact') | PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'), PKC_Ca_AA_kf, PKC_Ca_AA_kr)

# Raf activation

    # Raf activation by PKC_Ca_DAG_act
    Rule('Raf_a_1', Raf(Raf_b1=None, Raf_s='inact') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act') >> Raf(Raf_b1=None, Raf_s='act') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'), Raf_a_1_k)
    # Raf activation by PKC_DAG_AA_act
    Rule('Raf_a_2', Raf(Raf_b1=None, Raf_s='inact') +  PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='act') >> Raf(Raf_b1=None, Raf_s='act') +  PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='act'), Raf_a_2_k)
    # Raf activation by PKC_AA
    Rule('Raf_a_3', Raf(Raf_b1=None, Raf_s='inact') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact') >> Raf(Raf_b1=None, Raf_s='act') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'), Raf_a_3_k)
    # Raf activation by PKC_Ca_AA
    Rule('Raf_a_4', Raf(Raf_b1=None, Raf_s='inact') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact') >> Raf(Raf_b1=None, Raf_s='act') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'), Raf_a_4_k)
#
# Raf -> MEK -> ERK (from ref. 32)

    # MEK + Raf_a
    Rule('MEK_Raf_a', MEK(MEK_b1=None, MEK_p='p0') + Raf(Raf_b1=None, Raf_s='act') | MEK(MEK_b1=30, MEK_p='p0')%Raf(Raf_b1=30, Raf_s='act'), MEK_Raf_a_kf, MEK_Raf_a_kr)
    # MEK phosphorilation
    Rule('MEK_phospho', MEK(MEK_b1=30, MEK_p='p0')%Raf(Raf_b1=30, Raf_s='act') >> MEK(MEK_b1=None, MEK_p='p1') + Raf(Raf_b1=None, Raf_s='act'), MEK_phospho_k)
    # MEK_P + PP2A
    Rule('MEK_P_PP2A', MEK(MEK_b1=None, MEK_p='p1') + PP2A(PP2A_b1=None) | MEK(MEK_b1=30, MEK_p='p1')%PP2A(PP2A_b1=30), MEK_P_PP2A_kf, MEK_P_PP2A_kr)
    # MEK_P dephosphorilation
    Rule('MEK_P_dephospho', MEK(MEK_b1=30, MEK_p='p1')%PP2A(PP2A_b1=30) >> MEK(MEK_b1=None, MEK_p='p0') + PP2A(PP2A_b1=None), MEK_P_dephospho_k)
    # Raf_a + MEK_P
    Rule('MEK_P_Raf_a', MEK(MEK_b1=None, MEK_p='p1') + Raf(Raf_b1=None, Raf_s='act') | MEK(MEK_b1=30, MEK_p='p1')%Raf(Raf_b1=30, Raf_s='act'), MEK_P_Raf_a_kf, MEK_P_Raf_a_kr)
    # MEK_P phosphorilation
    Rule('MEK_P_phospho', MEK(MEK_b1=30, MEK_p='p1')%Raf(Raf_b1=30, Raf_s='act') >> MEK(MEK_b1=None, MEK_p='p2') + Raf(Raf_b1=None, Raf_s='act'), MEK_P_phospho_k)
    # MEK_PP + PP2A
    Rule('MEK_PP_PP2A', MEK(MEK_b1=None, MEK_p='p2') + PP2A(PP2A_b1=None) | MEK(MEK_b1=30, MEK_p='p2')%PP2A(PP2A_b1=30), MEK_PP_PP2A_kf ,MEK_PP_PP2A_kr)
    # MEK_PP dephosphorilation
    Rule('MEK_PP_dephospho', MEK(MEK_b1=30, MEK_p='p2')%PP2A(PP2A_b1=30) >> MEK(MEK_b1=None, MEK_p='p1') + PP2A(PP2A_b1=None), MEK_PP_dephospho_k)
    # Raf_a + PP2A
    Rule('Raf_a_PP2A', Raf(Raf_b1=None, Raf_s='act') + PP2A(PP2A_b1=None) | Raf(Raf_b1=30, Raf_s='act')%PP2A(PP2A_b1=30), Raf_a_PP2A_kf, Raf_a_PP2A_kr)
    # Raf_a deactivation
    Rule('Raf_a_deact', Raf(Raf_b1=30, Raf_s='act')%PP2A(PP2A_b1=30) >> Raf(Raf_b1=None, Raf_s='inact') + PP2A(PP2A_b1=None), Raf_a_deact_k)
    # ERK_PP + Raf_a
    Rule('ERK_PP_Raf_a', ERK(ERK_b1=None, ERK_p='p2') + Raf(Raf_b1=None, Raf_s='act') | ERK(ERK_b1=30, ERK_p='p2')%Raf(Raf_b1=30, Raf_s='act'), ERK_PP_Raf_a_kf, ERK_PP_Raf_a_kr)
    # Raf_a activation
    Rule('Raf_a_act', ERK(ERK_b1=30, ERK_p='p2')%Raf(Raf_b1=30, Raf_s='act') >> ERK(ERK_b1=None, ERK_p='p2') + Raf(Raf_b1=None, Raf_s='act2'), Raf_a_act_k)
    # Raf_aa + PP2A
    Rule('Raf_aa_PP2A', Raf(Raf_b1=None, Raf_s='act2') + PP2A(PP2A_b1=None) | Raf(Raf_b1=30, Raf_s='act2')%PP2A(PP2A_b1=30), Raf_aa_PP2A_kf, Raf_aa_PP2A_kr)
    # Raf_aa deactivation
    Rule('Raf_aa_deact', Raf(Raf_b1=30, Raf_s='act2')%PP2A(PP2A_b1=30) >> Raf(Raf_b1=None, Raf_s='act') + PP2A(PP2A_b1=None), Raf_aa_deact_k)
    # ERK + MEK_PP
    Rule('ERK_MEK_PP', ERK(ERK_b1=None, ERK_p='p0') + MEK(MEK_b1=None, MEK_p='p2') | ERK(ERK_b1=30, ERK_p='p0')%MEK(MEK_b1=30, MEK_p='p2'), ERK_MEK_PP_kf, ERK_MEK_PP_kr)
    # ERK phosphorilation
    Rule('ERK_phospho', ERK(ERK_b1=30, ERK_p='p0')%MEK(MEK_b1=30, MEK_p='p2') >> ERK(ERK_b1=None, ERK_p='p1') + MEK(MEK_b1=None, MEK_p='p2'), ERK_phospho_k)
    # ERK_P + MEK_PP
    Rule('ERK_P_MEK_PP', ERK(ERK_b1=None, ERK_p='p1') + MEK(MEK_b1=None, MEK_p='p2') | ERK(ERK_b1=30, ERK_p='p1')%MEK(MEK_b1=30, MEK_p='p2'), ERK_P_MEK_PP_kf, ERK_P_MEK_PP_kr)
    # ERK_P phosphorilation
    Rule('ERK_P_phospho', ERK(ERK_b1=30, ERK_p='p1')%MEK(MEK_b1=30, MEK_p='p2') >> ERK(ERK_b1=None, ERK_p='p2') + MEK(MEK_b1=None, MEK_p='p2'), ERK_P_phospho_k)

# ERK <-> MKP

    # ERK_PP + MKP
    Rule('ERK_PP_MKP', ERK(ERK_b1=None, ERK_p='p2') + MKP(MKP_b1=None, MKP_p='p0', MKP_s='inact') | ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p0', MKP_s='inact'), ERK_PP_MKP_kf, ERK_PP_MKP_kr)
    # MKP phosphorilation
    Rule('MKP_phospho', ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p0', MKP_s='inact') >> ERK(ERK_b1=None, ERK_p='p2') + MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact'), MKP_phospho_k)
    # ERK_PP + MKP_PP
    Rule('ERK_PP_MKP_PP', ERK(ERK_b1=None, ERK_p='p2') + MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact') | ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_PP_MKP_PP_kf, ERK_PP_MKP_PP_kr)
    # ERK_PP dephosphorilation
    Rule('ERK_PP_dephospho', ERK(ERK_b1=30, ERK_p='p2')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact') >> ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_PP_dephospho_k)
    # ERK_P_MKP_PP
    Rule('ERK_P_MKP_PP', ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact') | ERK(ERK_b1=None, ERK_p='p1') + MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact'), ERK_P_MKP_PP_kf, ERK_P_MKP_PP_kr)
    # MKP_PP degrade
    Rule('MKP_PP_deg', MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact') >> MKP(MKP_b1=None, MKP_p='p0', MKP_s='inact'), MKP_PP_deg_k)
    # MKP_PP activation
    Rule('MKP_PP_act', ERK(ERK_b1=None, ERK_p='p1') + MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact') | ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='act'), MKP_PP_act_kf, MKP_PP_act_kr)
    # ERK_P dephosphorilation
    Rule('ERK_P_dephospho', ERK(ERK_b1=30, ERK_p='p1')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='act') >> ERK(ERK_b1=30, ERK_p='p0')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact'), ERK_P_dephospho_k)
    # ERK_MKP_PP
    Rule('ERK_MKP_PP', ERK(ERK_b1=30, ERK_p='p0')%MKP(MKP_b1=30, MKP_p='p2', MKP_s='inact') | ERK(ERK_b1=None, ERK_p='p0') + MKP(MKP_b1=None, MKP_p='p2', MKP_s='inact'), ERK_MKP_PP_kf, ERK_MKP_PP_kr)

# RasGEF activation

    # RasGEF activation by PKC_Ca_DAG_act
    Rule('RasGEF_a_1', RasGEF(RasGEF_s='inact') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act') >> RasGEF(RasGEF_s='act') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'), RasGEF_a_1_k)
    # RasGEF activation by PKC_DAG_AA_act
    Rule('RasGEF_a_2', RasGEF(RasGEF_s='inact') + AA(AA_b1=50, AA_s='act')%PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact') >> RasGEF(RasGEF_s='act') + AA(AA_b1=50, AA_s='act')%PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact'), RasGEF_a_2_k)
    # RasGEF activation by PKC_AA
    Rule('RasGEF_a_3', RasGEF(RasGEF_s='inact') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact') >> RasGEF(RasGEF_s='act') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'), RasGEF_a_3_k)
    # RasGEF activation by PKC_Ca_AA
    Rule('RasGEF_a_4', RasGEF(RasGEF_s='inact') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact') >> RasGEF(RasGEF_s='act') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'), RasGEF_a_4_k)
    # Ras_GEF deactivation
    Rule('RasGEF_deact', RasGEF(RasGEF_s='act') >> RasGEF(RasGEF_s='inact'), RasGEF_deact_k)

# RasGAP activation

    # RasGAP activation by PKC_Ca_DAG_act
    Rule('RasGAP_a_1', RasGAP(RasGAP_s='inact') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act') >> RasGAP(RasGAP_s='act') + PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'), RasGAP_a_1_k)
    # RasGAP activation by PKC_DAG_AA_act
    Rule('RasGAP_a_2', RasGAP(RasGAP_s='inact') + AA(AA_b1=50, AA_s='act')%PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact') >> RasGAP(RasGAP_s='act') + AA(AA_b1=50, AA_s='act')%PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact'), RasGAP_a_2_k)
    # RasGAP activation by PKC_AA
    Rule('RasGAP_a_3', RasGAP(RasGAP_s='inact') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact') >> RasGAP(RasGAP_s='act') + PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'), RasGAP_a_3_k)
    # RasGAP activation by PKC_Ca_AA
    Rule('RasGAP_a_4', RasGAP(RasGAP_s='inact') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact') >> RasGAP(RasGAP_s='act') + PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'), RasGAP_a_4_k)
    # RasGAP deactivation
    Rule('RasGAP_deact', RasGAP(RasGAP_s='act') >> RasGAP(RasGAP_s='inact'), RasGAP_deact_k)

# Ras_GDP -> Raf_act

    # Ras_GTP activation by RasGEF_act
    Rule('Ras_GDP_act', Ras(Ras_b1=None, Ras_s='GDP') + RasGEF(RasGEF_s='act') >> Ras(Ras_b1=None, Ras_s='GTP') + RasGEF(RasGEF_s='act'), Ras_GTP_act_k)
    # # Ras_GTP deactivation by Ras_GAP_inact
    Rule('Ras_GTP_inact', Ras(Ras_b1=None, Ras_s='GTP') + RasGAP(RasGAP_s='inact') >> Ras(Ras_b1=None, Ras_s='GDP') + RasGAP(RasGAP_s='inact'), Ras_GTP_inact_k)
    # Ras_GDP <-> Ras_GTP
    Rule('Ras_GDP_Ras_GTP', Ras(Ras_b1=None, Ras_s='GDP') | Ras(Ras_b1=None, Ras_s='GTP'), Ras_GDP_Ras_GTP_kf, Ras_GDP_Ras_GTP_kr)
    # Ras_GTP + Raf
    Rule('Ras_GTP_Raf', Ras(Ras_b1=None, Ras_s='GTP') + Raf(Raf_b1=None, Raf_s='inact') | Ras(Ras_b1=70, Ras_s='GTP')%Raf(Raf_b1=70, Raf_s='inact'), Ras_GTP_Raf_kf, Ras_GTP_Raf_kr )
    # Raf activation by Ras_GTP
    Rule('Raf_act', Ras(Ras_b1=70, Ras_s='GTP')%Raf(Raf_b1=70, Raf_s='inact') >> Ras(Ras_b1=None, Ras_s='GTP') + Raf(Raf_b1=None, Raf_s='act'), Raf_act_k)

# OBSERVABLES

    Observable('obs_DAG', DAG(DAG_b1=None, DAG_s='inact'))
    Observable('obs_PA', PA())
    Observable('obs_Ca', Ca(Ca_b1=None))
    Observable('obs_Gq_a_GDP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GDP'))
    Observable('obs_Gq_a_GTP', Gq_a(Gq_a_b1=None, Gq_a_b2=None, Gq_a_s='GTP'))
    Observable('obs_Gq_trimer_GDP', Gq_a(Gq_a_b1=40, Gq_a_b2=None, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))
    Observable('obs_Gq_bg', Gq_bg(Gq_bg_b1=None, Gq_bg_b2=None))
    Observable('obs_RGS4_Gq_a_GTP', RGS4(RGS4_b1=50)%Gq_a(Gq_a_b1=50, Gq_a_b2=None, Gq_a_s='GTP'))
    Observable('obs_AA', AA(AA_b1=None, AA_s='inact'))
    Observable('obs_PKC_Ca', PKC(PKC_b1=30, PKC_b2=None)%Ca(Ca_b1=30))

    Observable('obs_PKC_Ca_DAG_act', PKC(PKC_b1=30, PKC_b2=40)%Ca(Ca_b1=30)%DAG(DAG_b1=40, DAG_s='act'))
    Observable('obs_PKC_DAG_AA_a', PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='act'))
    Observable('obs_PKC_DAG_AA', PKC(PKC_b1=50, PKC_b2=40)%DAG(DAG_b1=40, DAG_s='inact')%AA(AA_b1=50, AA_s='inact'))
    Observable('obs_PKC_AA', PKC(PKC_b1=50, PKC_b2=None)%AA(AA_b1=50, AA_s='inact'))
    Observable('obs_PKC_Ca_AA',  PKC(PKC_b1=30, PKC_b2=50)%Ca(Ca_b1=30)%AA(AA_b1=50, AA_s='inact'))

    Observable('obs_PIP2_PLA2_Ca', PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='inact', PIP2_b='True')%Ca(Ca_b1=30))
    Observable('obs_PIP2_PLA2', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='True'))
    Observable('obs_PLA2', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='inact', PIP2_b='False'))
    Observable('obs_PLA2_a', PLA2(PLA2_b1=None, PLA2_b2=None, PLA2_s='act', PIP2_b='False'))

    Observable('obs_PLA2_Ca', PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='inact', PIP2_b='False')%Ca(Ca_b1=30))
    Observable('obs_PLA2_a_Ca', PLA2(PLA2_b1=30, PLA2_b2=None, PLA2_s='act', PIP2_b='False')%Ca(Ca_b1=30))
    Observable('obs_PLA2_a_ERK_PP',  PLA2(PLA2_b1=None, PLA2_b2=40,  PLA2_s='act', PIP2_b='False')%ERK(ERK_b1=40, ERK_p='p2'))

    Observable('obs_PLCb_Ca_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=70)%Ca(Ca_b1=70))
    Observable('obs_PLCb', PLCb(PLCb_b1=None, PLCb_b2=None))
    Observable('obs_PLCb_Gq_a_GTP', Gq_a(Gq_a_b1=60, Gq_a_b2=None, Gq_a_s='GTP')%PLCb(PLCb_b1=60, PLCb_b2=None))
    Observable('obs_PLCb_Ca', PLCb(PLCb_b1=None, PLCb_b2=70)%Ca(Ca_b1=70))
    Observable('obs_ERK_PP', ERK(ERK_b1=None, ERK_p='p2'))
    Observable('obs_ERK_P', ERK(ERK_b1=None, ERK_p='p1'))
    Observable('obs_ERK', ERK(ERK_b1=None, ERK_p='p0'))

    Observable('obs_RasGEF_act', RasGEF(RasGEF_s='act'))
    Observable('obs_RasGAP_act', RasGAP(RasGAP_s='act'))
    Observable('obs_Ras_GTP_Raf', Ras(Ras_b1=30, Ras_s='GTP')%Raf(Raf_b1=30, Raf_s='inact'))
    Observable('obs_MEK', MEK(MEK_b1=None, MEK_p='p0'))

    Observable('obs_MKP_PP', MKP(MKP_b1=None, MKP_p='p2'))
    Observable('obs_Raf_act', Raf(Raf_b1=None, Raf_s='act'))
    Observable('obs_Raf', Raf(Raf_b1=None, Raf_s='inact'))

    Observable('obs_PP2A', PP2A(PP2A_b1=None))

    Observable('obs_L', L(L_b1=None))
    Observable('obs_RL', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%L(L_b1=50))
    Observable('obs_R', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane'))
    Observable('obs_RG', R5HT2A(R5HT2A_b1=30, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None))


    Observable('obs_R_in', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=None, R5HT2A_s='inact', R5HT2A_l='intra'))
    Observable('obs_RL_in', R5HT2A(R5HT2A_b1=None, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='intra')%L(L_b1=50))
    Observable('obs_RGL', R5HT2A(R5HT2A_b1=30, R5HT2A_b2=50, R5HT2A_s='act', R5HT2A_l='membrane')%Gq_a(Gq_a_b1=30, Gq_a_b2=40, Gq_a_s='GDP')%Gq_bg(Gq_bg_b1=40, Gq_bg_b2=None)%L(L_b1=50))


    return model

list_of_observables=['obs_RL','obs_ERK', 'obs_ERK_P', 'obs_ERK_PP', 'obs_DAG']
