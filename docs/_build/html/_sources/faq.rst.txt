.. _faq:

==========================
Frequently Asked Questions
==========================


General
=======

* If I sue the parameter `binding_kinetics = False` does the default pathway parameter `L_init` be considered?

    **No**. If `binding_kinetics = False` the reaction that rules the binding of ligand to the receptor will not be simulated

* Why am I obtaining same results when using `binding_kinetics = False` but changing other parameters of the pathway?

    Using `binding_kinetics = False` the reaction that rules the binding of ligand to the receptor will not be simulated, and so 
    the simulation results will be directly proportional to the parameter `LR` calculated using the Kd values.