Introduction to Structure Systems Biology
#########################################

Network and systems biology approaches have been already making important contributions to drug discovery and development and, 
despite being incomplete and error-prone, they are accurate enough to provide useful information 
(`Pujol et al., 2010 <https://doi.org/10.1016/j.tips.2009.11.006>`_). As a matter of fact, such approaches have been assisting 
in the identification and evaluation of drug targets, for instance in tumor associated diseases 
(`Wang G. et al., 2018 <https://doi.org/10.1007/978-1-4939-7456-6_11>`_, `DeWard et al., 2018 <https://doi.org/10.1007/978-1-4939-7493-1_13>`_). 

Even if network and systems biology approaches changed the way biomedical research thinks about drug action in complex diseases, 
the truth is that these approaches only focus on the interactions between its elements. This means that the quantitative element 
that is missing should be, undoubtedly, included into systems biology. In other words, not only network and systems biology 
should study the physical and functional interactions between the elements that constitute the biological system, but also should 
measure the concentrations and kinetic parameters that govern these interactions (`Pujol et al., 2010 <https://doi.org/10.1016/j.tips.2009.11.006>`_). 
Therefore, applying such concepts to pharmacology studies, the quantitative systems pharmacology rises up as a new discipline 
(`Androulakis, 2015 <https://doi.org/10.1002/wsbm.1294>`_).

Based on mathematical models, quantitative systems pharmacology provides, then, a framework that bridges spatial and temporal 
dimensions allowing a better understanding on how drugs affect complex biological systems and pathophysiological processes 
(`Krzyzanski, 2013 <https://doi.org/10.1038/psp.2013.15>`_). Such mathematical models rely on the fundamental concept of the 
receptor occupancy theory initially proposed by A. J. Clark, thereby making the central idea of quantitative systems pharmacology 
as the quantification of ligand-receptor interaction and subsequent cellular effects. The receptor occupancy theory states that on 
drug binding a cellular effect mediated by the activation of a signaling cascade may result. And the intensity of the cellular 
effect depends on drug-receptor parameters like specificity, affinity (strength and time of binding), and availability 
(drug and receptor concentration) (`Krzyzanski, 2013 <https://doi.org/10.1038/psp.2013.15>`_).

The cellular response mediated by those signaling pathways are, actually, an intricate series of molecular events, 
commonly proteins’ phosphorylation catalyzed by kinases, and each one of these events can be described by a mathematical equation. 
Typically, a mathematical model of a biochemical network consists of a set of ordinary differential equations (ODEs), 
and since ODEs depend only in one variable, they can be used to describe the change of the states of the system. For instance, 
if ODEs are integrated as a function of time, such reactions can describe how the concentration of species inside the network 
changes over time (`Stein et al., 2007 <https://doi.org/10.1016/j.sbi.2007.03.014>`_). However, the biochemical reactions underneath 
a signaling cascade are, normally, governed by kinetic parameters. While drug-receptor binding or protein-protein interactions are 
characterized by a second-order association constant, :math:`k_{on}` , and a first-order dissociation constant, :math:`k_{off}` 
(`Krzyzanski, 2013 <https://doi.org/10.1038/psp.2013.15>`_), the reactions of enzyme catalysis are characterized by :math:`K_m` , :math:`V_{max}` and 
:math:`K_{cat}` (`Stein et al., 2007 <https://doi.org/10.1016/j.sbi.2007.03.014>`_). Because under steady-state conditions these processes, from drugreceptor 
binding to downstream subcellular effects, occur at relatively short time scales equilibrium assumptions allow the explicit 
numerical calculation of concentration values of drugs, receptors, drug-receptor complexes and of all elements that make part 
of the signaling cascade (`Krzyzanski, 2013 <https://doi.org/10.1038/psp.2013.15>`_).

Mathematical models of signal-transduction cascades have been developed for decades for a variety of systems such as the signaling
mechanisms of hallucinogens towards the serotonin receptors (`Chang et al., 2009 <https://dx.doi.org/10.1016%2Fj.neuropharm.2008.07.049>`_), the 
sensing reward mechanism of dopaminergic receptors (`Nair et al., 2015 <https://dx.doi.org/10.1523%2FJNEUROSCI.0730-15.2015>`_), 
the mechanism of phototransduction mediated by rhodopsin receptors (`Dell'Orco et al., 2019 <https://doi.org/10.1039/B908123B>`_, 
`Invergo et al., 2013 <https://doi.org/10.1186/1478-811X-11-36>`_, `Invergo et al., 2014 <https://doi.org/10.1039/C3MB70584F>`_), 
or the signaling mediated by the epidermal growth factor receptor (`Kholodenko et al., 1999 <https://doi.org/10.1074/jbc.274.42.30169>`_). 
In fact, the number of mathematical models of signaling cascades have been growing so fast that many databases to store them have also been developed. 
Today, many repositories of pathways and mathematical models of biological and biomedical systems can be easily found across the internet, such as 
the BioModels database (`Glont et al., 2018 <https://doi.org/10.1093/nar/gkx1023>`_, `Malik-Sheriff et al., 2020 <https://doi.org/10.1093/nar/gkz1055>`_).

However, the dynamic modeling of biochemical networks is still hampered by the lack of kinetic parameters needed to feed the 
network. Ideally, such parameters should be determined experimentally under relevant conditions to the model; however, the truth 
is that the existing parameters are distributed in literature and relate to different experimental conditions 
(`Stein et al., 2007 <https://doi.org/10.1016/j.sbi.2007.03.014>`_, `Krzyzanski, 2013 <https://doi.org/10.1038/psp.2013.15>`_).

This is where protein structural data plays its role. Since all the kinetic parameters are encoded in the three-dimensional
macromolecular structure of proteins, such information can be used to estimate these parameters. In fact, at the same time
as the number of mathematical models grow, the number of protein structures that are being solved also rises concurrently 
(`Birch et al., 2020 <https://doi.org/10.3390/biology9110401>`_). For these reasons, the role of high-resolution three-dimensional 
protein structures in systems biology/pharmacology has become unquestionable, making many researchers claiming and advocating a new 
paradigm of structure systems pharmacology (`Xie et al., 2014 <https://dx.doi.org/10.1371%2Fjournal.pcbi.1003554>`_, 
`Duran-Frigola et al., 2013 <https://doi.org/10.1016/j.chembiol.2013.03.004>`_). The proof of this is the tremendous quantity of 
computationally methods, developed in the past years, to derive quantitative structure-kinetics relationships. For instance, 
molecular dynamics (MD) simulations in their many renditions have been used to estimate drug-receptor binding kinetics 
(`De Vivo et al., 2016 <https://doi.org/10.1021/acs.jmedchem.5b01684>`_, `Bruce N.J. et al., 2018 <https://doi.org/10.1016/j.sbi.2017.10.001>`_, 
`Nunes-Alves et al., 2020 <https://doi.org/10.1016/j.sbi.2020.06.022>`_), specially, metadynamics 
(`Capelli et al., 2020 <https://doi.org/10.1101/2020.03.30.015396>`_). Another approach is the comparison of molecular interaction fields 
by similarity indices. Assuming that the molecular interaction fields are the most relevant factors in determining the kinetic parameters 
values, such parameters can be transferable between related proteins (`Stein et al., 2007 <https://doi.org/10.1016/j.sbi.2007.03.014>`_). 
More recently, thanks to the power of machine and deep learning, many approaches using convolution neural networks have been promising 
the prediction of binding affinity data from protein-ligand complex structures (`Liu et al., 2017 <https://doi.org/10.1021/acs.accounts.6b00491>`_, 
`Ragoza et al., 2017 <https://doi.org/10.1021/acs.jcim.6b00740>`_).

Taken as a whole, structural-quantitative systems pharmacology provides an unprecedented molecular framework for understanding 
complex cell processes and molecular networks related to diseases. Bringing to light the dynamic interplay between those complex 
biological systems and drugs, systems pharmacology has been proposed as an alternative to overcome the critical steps in the drug 
discovery pipeline (`Duran-Frigola et al., 2013 <https://doi.org/10.1016/j.chembiol.2013.03.004>`_), promising then to be 
next-generation of drug-discovery and personalized medicine.
