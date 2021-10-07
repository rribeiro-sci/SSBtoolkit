import math
from pysb import *
from pysb.macros import *

# define 'network()' function
def network(pKd,ligand,epsilon):

        ##### PHOTOACTIVATION       #####
    ''' The stimulus is defined as the sum  of potential backgroung illumination,
    pre- and test-flashes, and another stimulus allowing for a particular setting
    employed in certain experiments: background + preflash + testflash + otherstimulus'''
    #total receptor
    Rtotal = 6.3102E10
    #from pKd to Kd
    Kd = 10**(-(float(pKd)))*6.022E23
    #LR determination
    a = 1
    b = float(ligand)+float(Rtotal)+Kd
    c = float(Rtotal)*float(ligand)
    delta = (b**2) - (4*a*c)
    LR = (b-math.sqrt(delta))/(2*a)
    #print(Rtotal,ligand,Kd)
    #print('The solutions are {0}'.format(LR))
    #Rho_activated = stimulus*Rtotal
    #Rho_activated = float(output of equation)
    #ttotal = float(inp4)
    # create a Model
    Model()

    # MONOMERS
    Monomer('R')
    Monomer('Ga', ['b1', 'b2','b3', 's'],  {'s': ['i', 'a']})
    Monomer('Gbg', ['b'])
    Monomer('GDP', ['b2'])
    Monomer('Rho', ['b1', 'p', 's'], {'p':['O','i'], 's':['i','a']}) # p is the phospholation state and s is the activation state
    Monomer('GTP', ['b2'])
    # Assuming only one alosteric binding site (SINGLE ACTIVATION)
    Monomer('PDE', ['site', 'nGGTP', 's'], {'nGGTP':['zero','one','two'], 's': ['i', 'a', 'aa']})
    Monomer('Rec', ['b', 'b1', 's'], {'s':['R', 'T']})
    Monomer('RK', ['b'])
    Monomer('Ca', ['b', 'l', 's'], {'l':['cytos', 'ext'], 's':['free' , 'buff']}) # l (location), s (state: free or buffed)
    Monomer('Arr', ['b', 'p1', 'p2'])
    Monomer('Ops')
    Monomer('RGS', ['b'])
    Monomer('cGMP', ['b'])
    Monomer('GMP', ['b'])

    # INITIAL CONDITIONS
    #Initial(R, Parameter('R_0', R_activated))                          #artificial added by me (molecules)
    Initial(Ga(b1=None, b2=None, b3=None, s='i'), Parameter('Ga_0', 3.6E8)) # number |Total amount of Gat (total amount of G [Dell'orco et al. 2009])
    Initial(Gbg(b=None), Parameter('Gbg_0', 3.6E8))                    # number | Total amount of Gbg (total amount of G [Dell'orco et al. 2009])
    Initial(GDP(b2=None), Parameter('GDP_0', 0.00315E22))              #molecules artificial add be me
    Initial(Rho(b1=None, p='O', s='i'), Parameter('Rho_i_0', 3.6E9 - LR))   # number |Total amount of R [Dell'orco et al. 2009]
    Initial(Rho(b1=None, p='O',s='a'), Parameter('Rho_a_0', LR))
    Initial(GTP(b2=None), Parameter('GTP_0', 0.00315E22))              #molecuels artificial add by me
    Initial(PDE(site=None, nGGTP='zero', s='i'), Parameter('PDE_0', 1.335E7)) # number | Total amount of PDE6 tetramers [Dell'orco et al. 2009]
    Initial(Rec(b=None, b1=None, s='T'), Parameter('RecT_0', 35))      # µM |Total concentration of Rec [Dell'orco et al. 2009]
    Initial(RK(b=None), Parameter('RK_0', 7))                          # µM |Total concentration of RK [Dell'orco et al. 2009]
    Initial(Ca(b=None, l='cytos', s='free'), Parameter('Ca_0', 0.25))  # µM |Dark Ca concentration [invergo et al. 2014]
    Initial(Ca(b=None, l='ext', s='free'), Parameter('Ca_ext_0', 10))  # artificially added by me ; µM
    Initial(Arr(b=None,p1=None, p2=None), Parameter('Arr_0', 3.13E7))  # number |Total amount of Arr [Dell'orco et al. 2009]
    Initial(Ops(), Parameter('Ops_0', 0))                              # number |Total amount of Ops (which is 0 in the steady-state)
    Initial(RGS(b=None), Parameter('RGS_0', 3E6))                      # number |Total amount of RGS9 [Dell'orco et al. 2009]
    Initial(cGMP(b=None), Parameter('cGMP_0', 6.5))                    # µM |Dark cGMP concentration [invergo et al. 2014]
    Ca_min = 0.023                                                     # µM |Minimum intracellular Ca2+ concentration [invergo et al. 2014]

    # RATE PARAMETERS
    Parameter('kc', 1e-10)          #
    Parameter('kGrecyc', 2)         # s-1 |Binding rate for GaGDP to Gbg
    Parameter('kG1', 1e-3)          # s-1 |Binding rate of Gt to unphosphorylated R* [invergo et al. 2014]
    Parameter('kG2', 2200.0)        # s-1 |Dissociation rate of the R*.Gt complex [invergo et al. 2014]
    Parameter('kG3', 8500.0)        # s-1 |Dissociation rate of GDP from R*.Gt complex [invergo et al. 2014]
    Parameter('kG4', 400.0)         # s-1 |Association rate of the R*.Gt complex [invergo et al. 2014]
    Parameter('kG5', 3500.0)        # s-1 |Association rate of GTP to the R*.Gt complex [invergo et al. 2014]
    Parameter('kG6', 8500.0)        # s-1 |Dissociation rate of the R*.GGTP complex [invergo et al. 2014]
    Parameter('kG7', 200.0)         # s-1 |Dissosication rate of GGTP into Gbg and GaGTP
    Parameter('kP1', 0.05497)       # s-1 |Binding rate of PDE to GaGTP [invergo et al. 2014]
    Parameter('kP1prev', 0.0)       # s-1 |Dissociation rate of PDE.GaGTP, without PDE activation [invergo et al. 2014]
    Parameter('kP2',940.7)          # s-1 |Rate activation of the first PDE subunit of PDE.GaGTP [invergo et al. 2014]
    Parameter('kP3',1.4983E-9)      # s-1 |Binding rate of GaGTP to an active PDE.GaGTP complex [invergo et al. 2014]
    Parameter('kP4',21.088)         # s-1 |Rate of activation of the second PDE subunit of GaGTP.PDE.GaGTP [invergo et al. 2014]
    Parameter('kPDEprev', 0.0)      # s-1 |
    Parameter('kRec1', 0.011)       # µM-1*s-1 |Rate of Ca2+ -triggered Rec conformational change (tense to relaxed) [invergo et al. 2014]
    Parameter('kRec2', 0.05)        # s-1 |Rate of Rec conformational change (relaxed to tense) [invergo et al. 2014]
    Parameter('kRec3', 4.1081E-4)   # s-1 |Binding rate of Rec.Ca2+ to RK [invergo et al. 2014]
    Parameter('kRec4', 0.610084)    # s-1 |Dissociation rate of RK from Rec-Ca2+ [invergo et al. 2014]
    Parameter('kRK1', 0.1724)       # s-1 |Binding rate of RK to unphophorylated R* [invergo et al. 2014]
    Parameter('kRK2', 250)          # s-1 |Dissociation rate of the R* from RK prior to phosphorylation [invergo et al. 2014]
    Parameter('kRK3', 4000)         # s-1 |Binding rate of ATP to R*.RK [invergo et al. 2014]
    Parameter('kRK4', 250)          # s-1 |Dissociation rate of R* from the R*.RK complex following phosphorylation [invergo et al. 2014]
    Parameter('kA1', 9.9147E-6)     # s-1 |Binding rate of Arr to single-phophorylated R* [invergo et al. 2014]
    Parameter('kA2', 0.026)         # s-1 |Dissociation rate of R* from the Arr.R* complex prior to R* inactivation [invergo et al. 2014]
    Parameter('kA3', 1.1651)        # s-1 |Dissociation rate of R* from the Arr.R* complex following R* inactivation [invergo et al. 2014]
    Parameter('kA4', 2.9965E-7)     # s-1 |Binding rate of Arr to from homo-oligomers [invergo et al. 2014]
    Parameter('kA5', 0.424)         # s-1 |Dissociate rate of Arr from homo-oligomers [invergo et al. 2014]
    Parameter('ktherm', 0.0238)     # s-1 |Thermal decay of R* [invergo et al. 2014]
    Parameter('kRecyc', 7E-4)       # s-1 |Rate constant for regeneration from Ops [invergo et al. 2014]
    Parameter('kGRS1', 4.8182E-5)   # s-1 |Binding rate of RGS9-1 to PDE-GaGTP(one or both active subunits)[invergo et al. 2014]
    Parameter('kGRS2', 98.0)        # s-1 |Rate of hydorlysis and dissociation od one PDE subunit from GaGTP [invergo et al. 2014]
    Parameter('kGshutoff', 0.05)    # s-1 |Rate of GaGTP auto-catalytic GTPase activity [invergo et al. 2014]
    Parameter('Bdark', 3.19)        # s-1 |Dark rate of cGMP hydrolysis[invergo et al 2014]
    Parameter('Bsub', 2.1826E-3)    # s-1 |Rate constant for one catalytic PDE subunit [invergo et al 2014]
    Parameter('kPDEshutoff', 0.1)   # s-1 |Rate of PDE-induced spontaneous PDE.GaGTP shutoff [invergo et al. 2014]
    fCa = 0.12                      #     |Fraction of the circulating current carried by Ca2+ [invergo et al 2014]
    Jdark = 14.87                   # pA  |Dark circulating current [invergo et al 2014]
    F = 96485.34                    # cm-1|Faraday constant [Dell'orco et al. 2009]
    n_cg = 3.8                      #     |Hill coefficient for opening cGMP-gated ion channels [invergo et al 2014]
    V_cyto = 0.03916                # pL  |Outer segment cytoplasmatic volume [invergo et al 2014]
    YCa = 981.3558                  # s-1 |Rate of Ca2+ extrusion by Na+/Ca2+k+ ion exchanger [invergo et al 2014]
    a_max = 60                      # µMs-1 |Maximal rate of cGMP synthesis  [invergo et al 2014]
    m1 = 3.0                        #     |Hill coefficient for GCAP1-mediated Ca2+ feedback on GC activity [invergo et al 2014]
    m2 = 1.5                        #     |Hill coefficient for GCAP2-mediated Ca2+ feedback on GC activity [invergo et al 2014]
    Kc1 = 171                       # nM  |EC50 for GCAP1-mediated Ca2+ feedback on GC activity [invergo et al 2014]
    Kc2 = 59                        # nM  |EC50 for GCAP2-mediated Ca2+ feedback on GC activity [invergo et al 2014]
    eT = 400                        # µM  |Total Ca buffer molecules concentration [invergo et al 2014]
    Parameter('k1', 9.37059)        # µM-1s-1 |Binding rate of Ca2+ to cytoplasmatic buffers [invergo et al 2014]
    Parameter('k2', 46.412)         # s-1 |Dissociaiton rate of Ca2+ from cytoplasmatic buffers [invergo et al 2014]


    #####       REACTIONS       #####
    ##### CASCADE ACTIVATION #####
    ''' When rhodopsin is activated binds to the G protein inactivated forming the
    complex R:GGDP. When the later is formed there is a nucleotide exchange: GDP to GTP.
    When GTP is bonded to the G protein the complex GGTP dissosiates from the receptor
    and becomes activated, and consequently, the Gproteins subunits will dissociate.
    G monomer have to be 2 binding sites: b1 and b2, and a state site s.
    The b1 binds Rho and b2 binds GDP/GTP. The state is inactive (when GDP is bonded)
    or active (when GTP is bonded. b1 for binding Rho and b2 for binding GDP/GTP.'''


    # The Ga will bind to GDP forming GaGDP
    Rule('G_bind_GDP', Ga(b1=None, b2=None, b3=None, s='i') + GDP(b2=None)
         >>Ga(b1=None, b2=1, b3=None, s='i') % GDP(b2=1), kc)

    # The GaGDP bind to Gbg
    Rule('G_hetrotrimer', Ga(b1=None, b2=1, b3=None, s='i') % GDP(b2=1) + Gbg(b=None) >> Ga(b1=None, b2=1, b3=3, s='i') % GDP(b2=1) % Gbg(b=3), kGrecyc)

    # The Rho binds to GGDP
    bind_complex(Rho(b1=None, p='O',s='a'), 'b1',Ga(b1=None, b2=1, b3=3, s='i') % GDP(b2=1) % Gbg(b=3), 'b1', [kG1, kG2])

    # GDP dissociation from the Rho.Gt complex
    Rule('RhoGGDP_dissociation', Rho(b1=50, p='O',s='a') % Ga(b1=50, b2=1, b3=3, s='i') % GDP(b2=1) % Gbg(b=3) | Rho(b1=50,p='O', s='a') % Ga(b1=50, b2=None, b3=3, s='i') % Gbg(b=3) + GDP(b2=None), kG3, kG4)

    # GTP binding to the Rho.Gt complex
    Rule('Rho_bind_GGTP', Rho(b1=50,p='O', s='a') % Ga(b1=50, b2=None, b3=3, s='i') % Gbg(b=3) + GTP(b2=None) >> Rho(b1=50,p='O', s='a') % Ga(b1=50, b2=1,b3=3, s='i') % Gbg(b=3) % GTP(b2=1), kG5)

    # Dissociation of the Rho.G.GTP
    Rule('RhoGGTP_dissociation', Rho(b1=50,p='O', s='a') % Ga(b1=50, b2=1,b3=3, s='i') % Gbg(b=3) % GTP(b2=1) >> Rho(b1=None,p='O', s='a') + Ga(b1=None, b2=1, b3=3, s='i') % Gbg(b=3) % GTP(b2=1), kG6)

    # Activation of GaGTP
    Rule('GaGTP_act', Ga(b1=None, b2=1, b3=3, s='i') % Gbg(b=3) % GTP(b2=1) >> Ga(b1=None, b2=1, b3=None, s='a') % GTP(b2=1) + Gbg(b=None), kG7)
    #equilibrate(G(b1=None, b2=1, s='i') % GTP(b2=1), G(b1=None, b2=1, s='a') % GTP(b2=1), [kG7, kG7_1])


    ##### PDE ACTIVATION #####
    '''When GaGTP is activated, will bind to one of the two binding sites of PDE, activating it.
       However, o be fully activated, PDE as to bind to two GaGTP. In any case, single or double activated,
       PDE will convert cGMP into GMP. The PDE function is regulates by PDE itself and by regulatory proteins, like RGS.'''

    # Binding of GaGTP to one PDE inactive subunit
    Rule('Bind_PDE', PDE(site=None, nGGTP='zero', s='i') + Ga(b1=None, b2=1, b3=None, s='a') % GTP(b2=1) >> PDE(site=None, nGGTP='one', s='i'), kP1)

    # Activation of the PDE.GaGTP complex
    equilibrate(PDE(site=None, nGGTP='one', s='i'), PDE(site=None, nGGTP='one', s='a'), [kP2, kPDEprev])

    # Binding of GaGTP to one PDE active subunit
    Rule('Bind_PDE_2', PDE(site=None, nGGTP='one', s='a') + Ga(b1=None, b2=1, b3=None, s='a') % GTP(b2=1) >> PDE(site=None, nGGTP='two', s='a'), kP3)

    # Activation of both the GaGTP-bounded PDE subunits
    equilibrate(PDE(site=None, nGGTP='two', s='a'), PDE(site=None, nGGTP='two', s='aa'), [kP4, kPDEprev])

    #Hydrolysis of cGMP in dark
    Rule('cGMP_hyd_dark', cGMP(b=None) >> GMP(b=None), Bdark)

    #Hydrolysis of cGMP by PDE
    Expression('kPDE1', Bdark + Bsub)
    Expression('kPDE2', Bdark + 2 * Bsub)
    catalyze_one_step(PDE(site=None, nGGTP='one', s='a'), cGMP(b=None), GMP(b=None), kPDE1)
    catalyze_one_step(PDE(site=None, nGGTP='two', s='aa'), cGMP(b=None), GMP(b=None), kPDE2)

    #PDE deactivation by RGS
    # Binding of RGS complex to the completely activated PDE tetramer
    Rule('rule_RGS2', PDE(site=None, nGGTP='two', s='aa') + RGS(b=None) >> PDE(site=1, nGGTP='two', s='aa') % RGS(b=1), kGRS1)

    # RGS-mediated deactivation of one of the two PDE active subunits
    Rule('PDE2_deactivation', PDE(site=1, nGGTP='two', s='aa') % RGS(b=1) >> PDE(site=None, nGGTP='one', s='a') + RGS(b=None) + Ga(b1=None, b2=1, b3=None, s='a') % GDP(b2=1), kGRS2)

    # Binding of RGS complex to the PDE tetramer with one activat subunit
    Rule('rule_RGS1', PDE(site=None, nGGTP='one', s='a') + RGS(b=None) >> PDE(site=1, nGGTP='one', s='a') % RGS(b=1), kGRS1)

    # RGS-mediated deactivation of the PDE onle active subunit
    Rule('PDE1_deactivation', PDE(site=1, nGGTP='one', s='a') % RGS(b=1) >> PDE(site=None, nGGTP='zero', s='i') + RGS(b=None) + Ga(b1=None, b2=1, b3=None, s='a') % GDP(b2=1), kGRS2)

    # Inactivation of the PDE*.GaGTP complex by GaGTP's GTPase activity
    Rule('PDE_shutoff2', PDE(site=None, nGGTP='two', s='aa') >> PDE(site=None, nGGTP='one', s='a') + Ga(b1=None, b2=1, b3=None, s='a') % GDP(b2=1), kPDEshutoff)

    # Inactivation of one of the two active PDE subunits by GaGTP's GTPase activity
    Rule('PDE_shutoff1', PDE(site=None, nGGTP='one', s='a') >> PDE(site=None, nGGTP='zero', s='i') + Ga(b1=None, b2=1, b3=None, s='a') % GDP(b2=1), kPDEshutoff)

    # GaGTP auto-catalytic GTPase activity
    Rule('GaGTP_shutoff', Ga(b1=None, b2=1, b3=None, s='a') % GTP(b2=1) >> Ga(b1=None, b2=1, b3=None, s='i') % GDP(b2=1), kGshutoff)


    ##### Rec pathway #####
    '''Rec protein in two forms: R and T. When Rec binds two Ca2+ ions, it undergoes
    a conformational change, exposing a myristoyl froup- This "relaxed" (RecRca) from
    preveils in the dark, at higher intracellular Ca2+ concentrations, and increases
    the proteins affinity for the disk membrane due to an overall augmented hydrophobicity.
    In the relaxed for, Rec binds RK and prevents it from phosphorylating the activated
    photopigment rhodopsin (Rho). As Ca2+ concentration decreases during a photoresponse,
    Rec reverst to ist tense form (RecT), becoming more hydrophilic, and leading to bind RK.
    [Invergo et al. 2013]'''

    # Ca2+ -induced Rec conformation change
    Rule('RecT', Rec(b=None, b1=None, s='T') + Ca(b=None, l='cytos', s='free') | Rec(b=1, b1=None, s='R') % Ca(b=1, l='cytos', s='buff'), kRec1, kRec2)

    # Binding of RK to Rec
    bind_complex(Rec(b=1, b1=None, s='R') % Ca(b=1, l='cytos', s='free'), 'b1',RK(b=None), 'b', [kRec3, kRec4])

    ##### RK PATHWAY #####
    # Binding of Rho and Rk.
    bind(RK, 'b', Rho(b1=None, p='O',s='a'), 'b1', [kRK1, kRK2])

    #Phosphorilation of Rho by RK
    Rule('Rho_phosphorilation', RK(b=50)%Rho(b1=50, p='O', s='a') >> RK(b=50)%Rho(b1=50, p='i',s='a'), kRK3)

    #Dissociation of the complex RK.Rho
    Rule('Rho_RK_dissociation', RK(b=50)%Rho(b1=50, p='i',s='a') >> RK(b=None) + Rho(b1=None, p='i',s='a'), kRK4)


    ###### CALCIUM #####
    '''The decreasinng in Ca2+ concentrations trigger the activity of two
    guanylate cyclase activatin proteins (GCAPs), which cause Guanylate cyclases (CGs)
    to synthesize cGMP at higher rates. This leads to the re-opening of the cGMP-gated
    ion channels and a retur to the dark circulating current. The decrease on cGMP leads
    to the closure of the cGMP-gates ion channels, a subsquente drop in the intracellular
    Ca2+ concentreation due to its continued extrosion via Na+/Ca2+, K+ exchangers, and a
    consequent hyper-polarization do the cell membrane. In the cytosol some Ca2+ is buffered
    by some proteins, including the Rec protein.'''

    #Calcium association and dissociation from intracellular buffers with total concentration eT
    Observable('Ca_buff', Ca(b=None, l='cytos', s='buff'))
    Observable('Ca_free', Ca(b=None, l='cytos', s='free'))
    Expression('kbuff1', k1 * (eT - Ca_buff) * Ca_free)
    Rule('Ca_buffering', Ca(b=None, l='cytos', s='free') | Ca(b=None, l='cytos', s='buff'), kbuff1, k2)

    # Extracellular Ca2+ influx via the cGMP_gated cation channels
    Observable('cGMP_conc', cGMP(b=None))
    Expression('kCa_influx', ((10E6 * fCa * Jdark) / ((2 + fCa) * F * V_cyto)) * ((cGMP_conc / cGMP_0) ** n_cg))
    Rule('Ca_influx', Ca(b=None, l='ext', s='free') >> Ca(b=None, l='cytos', s='free') , kCa_influx)

    # Intracellular Ca2+ efflux via the Na+/Ca2+K+ exchanger
    Expression('kCa_efflux', YCa * (Ca_free - Ca_min))
    Rule('Ca_efflux', Ca(b=None, l='cytos', s='free') >> Ca(b=None, l='ext', s='free'), kCa_efflux)

    ### cGMP synthisis by guanylate cyclase
    Expression('kGC', a_max / (1 + ((Ca_free / Kc1) ** m1)) + a_max / (1 + ((Ca_free / Kc2) ** m2)))
    Rule('cGMP_synt', GMP(b=None) >> cGMP(b=None), kGC)


    ##### Arr PATHWAY #####
    '''Arrestin binds to the activated receptor and phophorylate it.
    The phospholyrated receptor will be inactivated and will loose the retinal
    becoming the so called Opsonin (Ops). The Ops is recycled to the full inactivated receptor'''
    # Binding of Arr and Rho
    bind(Rho(b1=None, p='i',s='a'), 'b1', Arr(b=None,p1=None, p2=None), 'b', [kA1, kA2])

    # Arr -mediated inactivation of Rho. Ops indicates the ligand-free receptor
    Rule('Rho_Ops', Rho(b1=50, p='i',s='a')%Arr(b=50) >> Ops() + Arr(b=None,p1=None, p2=None), kA3)

    # Thermal decay of catalyti
    Rule('Rho_decay', Rho(b1=None, p='O',s='a') >> Ops(), ktherm)

    #Recycle of Ops to R (not activated by photons)
    Rule('Rho_recycle', Ops() >> Rho(b1=None,p='O',s='i'), kRecyc) # we need to add the receptor activation reaction, otherwise the concentration of R will be artificial increased.

    ## ARR homo-tetramerization
    assemble_pore_sequential(Arr(b=None,p1=None, p2=None),'p1','p2', 4,[[kA4,kA5],[kA4,kA5],[kA4,kA5]])
    Rule('Dimer', Arr(b=None, p1=None, p2=None) + Arr(b=None, p1=None, p2=None) | Arr(b=None, p1=None, p2=1) % Arr(b=None, p1=1, p2=None), kA1, kA2)
    Rule('tetramerization', Arr(b=None, p1=None, p2=1) % Arr(b=None, p1=1, p2=None) + Arr(b=None, p1=None, p2=1) % Arr(b=None, p1=1, p2=None) | Arr(b=None, p1=1, p2=2) % Arr(b=None, p1=2, p2=3) % Arr(b=None, p1=1, p2=4) %  Arr(b=None, p1=4, p2=3), kA1, kA2)

    ######################   RESULTS ########################
    ## Look at the observables list file

    # OBSERVABLES

    #Tracking G protein activation by activated receptor
    Observable('obs_Rho_i', Rho(b1=None, p='O',s='i'))                              # Rho               | Inactivated Receptor
    Observable('obs_Rho_a', Rho(b1=None, p='O',s='a'))                              # Rho               | Activated Receptor
    #Observable('obs_G', G(b1=None, b2=None, s='i'))                                 # G                 | G protein without nucleotide
    #Observable('obs_GGDP', G(b1=None, b2=1, s='i') % GDP(b2=1))                     # G.GDP             | G protiens complexed with GDP
    #Observable('obs_RhoGGDP', Rho(b1=50, p='O',s='a')% G(b1=50, b2=1, s='i') % GDP(b2=1))   # Rho.GGDP  | Receptor binded to GGDP
    Observable('obs_GDP', GDP(b2=None))                                             # GDP
    #Observable('obs_RhoG', Rho(b1=50,p='O',s='a') % G(b1=50, b2=None, s='i'))               # Rho.G     | Receptor binded to G protein without GDP
    Observable('obs_GTP', GTP(b2=None))                                             # GTP
    #Observable('obs_RhoGGTP', Rho(b1=50,p='O',s='a') % G(b1=50, b2=1, s='i') % GTP(b2=1))   # Rho.GGTP  | Receptor binded to G protein with GTP
    #Observable('obs_GGTP', G(b1=None, b2=1, s='i') % GTP(b2=1))                     # GGTP              | G protein complexed with GTP
    #Observable('obs_GaGTP', G(b1=None, b2=1, s='a') % GTP(b2=1))                    # GaGTP             | G alpha subunit with GTP
    #Tracking PDE pathway
    Observable('obs_PDE', PDE(site=None, nGGTP='zero', s='i'))                      # PDE               | Inactivated PDE
    Observable('obs_PDE_GaGTP', PDE(site=None, nGGTP='one', s='i'))                 # PDE.GaGTP         | PDE binded to one GaGTP without activation
    Observable('obs_PDEGaGTP', PDE(site=None, nGGTP='one', s='a'))                  # PDE*.GaGTP        | PDE single activation with one GaGTP
    Observable('obs_GaGTPPDE_GaGTP', PDE(site=None, nGGTP='two', s='a'))            # GaGTP.PDE*.GaGTP  | PDE single activated binded to two GaGTP
    Observable('obs_GaGTPPDEGaGTP', PDE(site=None, nGGTP='two', s='aa'))            # GaGTP.*PDE*.GaGTP | PDE double activated without binded to two GaGTP
    Observable('obs_cGMP', cGMP(b=None))                                            # cGMP              | Free cGMP
    Observable('obs_GMP', GMP(b=None))                                              # GMP               | Free GMP
    Observable('obs_RGS', RGS(b=None))                                              # RGS               | Free RGS
    #Observable('obs_GaGDP_a', G(b1=None, b2=1, s='a') % GDP(b2=1))                  # GaGTP             | G aplha actitaved with GTP
    #Observable('obs_GaGDP_i', G(b1=None, b2=1, s='i') % GDP(b2=1))                  # GaGTP             | H alpha inactivated with GTP

    Observable('obs_PDEaa_RGS', PDE(site=1, nGGTP='two', s='aa') % RGS(b=1))        # PDE**:RGS         | PDE double activated complexed with RGS
    Observable('obs_PDEa_RGS', PDE(site=1, nGGTP='one', s='a') % RGS(b=1))          # PDE*:RGS          | PDE single actitaved complexed with RGS
    # Tracking CALCIUM
    Observable('obs_Ca_free', Ca(b=None, l='cytos', s='free'))                      # Ca++ Free         | Calcuim free in the cytosol
    Observable('obs_Ca_buff', Ca(b=None, l='cytos', s='buff'))                      # Ca++ buff         | Calcium buffered in the cytosol (except buffered by REC)
    Observable('obs_Ca_ext', Ca(b=None, l='ext', s='free'))                         # Ca++ ext          | Extracellular Calcium
    #Tracking REC protein and RK
    Observable('obs_RecT', Rec(b=None, b1=None, s='T'))                             # RecT              | REC protein in the Tense form
    Observable('obs_RecR_CA', Rec(b=1, b1=None, s='R') % Ca(b=1, l='cytos', s='buff'))        #           | Rec binds to Calium (Rec R from and caluim buffered)
    Observable('obs_Rec_Ca_RK', Rec(b=1, b1=50, s='R')%Ca(b=1, l='cytos', s='free')%RK(b=50)) # REC:Ca:RK | Rec in the R from binded to Ca and RK
    Observable('obs_RK', RK(b=None))
    #Tracking Recptor phosporilation by RK
    Observable('obs_RK_Rho',RK(b=50)%Rho(b1=50, p='O',s='a'))                       # RK:Rho            | RK binded to activated receptors
    Observable('obs_RK_Rho1', RK(b=50)%Rho(b1=50, p='i',s='a'))                     # RK:Rho 1          | RK binded to phosporilated activated receotors
    Observable('obs_Rho1', Rho(b1=None, p='i',s='a'))                               # Rho1              | phoporilated receptor
    #Tracking Receptor regulation by Arrestin
    Observable('obs_Arr1', Arr(b=None,p1=None, p2=None))                            # ARR1              | Arrestin Free
    Observable('obs_Arr2', Arr(b=None, p1=None, p2=1) % Arr(b=None, p1=1, p2=None)) # ARR2              | Dimers of arrestin
    #Observable('obs_Arr4', Arr(b=None, p1=1, p2=2)%Arr(b=None, p1=2, p2=3)%Arr(b=None, p1=1, p2=4)%Arr(b=None, p1=4, p2=3), kA1, kA2)   # Arr4 | Tetramers of Arrestin
    Observable('obs_Arr_Rho1', Rho(b1=50, p='i',s='a')%Arr(b=50,p1=None, p2=None))    # ARR:Rho1          | Arresntin bindes to the phosporilated receptor
    Observable('obs_Ops', Ops())                                                    # Ops               | Opsonin - Inactivated receptor (without retinal)

    return model

list_of_observables=['obs_Rho_i', 'obs_Rho_a', 'obs_GDP', 'obs_GTP','obs_PDE', 'obs_PDE_GaGTP', 'obs_PDEGaGTP',
'obs_GaGTPPDEGaGTP','obs_cGMP','obs_GMP','obs_RGS','obs_PDEaa_RGS','obs_PDEa_RGS','obs_Ca_free','obs_Ca_buff']
