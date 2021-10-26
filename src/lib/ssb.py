#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Rui Ribeiro
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This program is a collection of utilities necessary to perform the QSP calculations
"""
__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"

#system libraries
import sys
import os
import csv
import time
import json
import subprocess
import importlib
import urllib.request   as urllib
from glob               import glob
from xml.dom            import minidom
from urllib.request     import urlopen

#Paralleling
import multiprocessing as mp

#Scientific Libraries
import math
import numpy  as np
import pandas as pd
from scipy.spatial    import distance
from scipy.optimize   import curve_fit, minimize

#Plotting Libraries
import plotly.graph_objs as go
import plotly.io         as pio
import pylab             as pl

#PYSB
from pysb import *
from pysb.macros import *
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator

#biopython libraries
from Bio                        import pairwise2, SeqIO
from Bio.SubsMat.MatrixInfo     import blosum62
from Bio.Blast                  import NCBIWWW
from Bio.PDB.PDBParser          import PDBParser
from Bio.PDB                    import PDBIO
from Bio.PDB.vectors            import Vector
from Bio.PDB.Superimposer       import Superimposer
from Bio.Blast.Applications     import NcbiblastpCommandline
from Bio.Seq                    import Seq
from Bio.SeqRecord              import SeqRecord

#direcotires
from src.lib.directories import *
#from src.lib.autogrids.auto_grids import auto_grids


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from sklearn.preprocessing import minmax_scale

#PROGRAMS
prepare_receptor = prepare_receptor_executable
prepare_receptor_flex = prepare_receptor_flex_executable
prepare_ligand = prepare_ligand_executable
vina = vina_executable
vina_split = vina_split_executable
dlscore_path = dlscore_path
reference_files_dir = autogrids_reference_dir
bp_center_reference_files='/autogrid_reference_files/'

#Human Target Receptors DataBase directory PATH
HuTRdb_path = HuTRdb_path
swissprot_db = swissprot_db

def help():
    print('Info should be written here')


class sharedfunctions:
    def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total:
            print()

    def LR_eq_conc(receptor_conc, agonist_conc, antagonist_conc, pkd_agonist, pkd_antagonist):
        #from pKd to Kd
        kd_anatagonist = (10**(-(float(pkd_antagonist)))) * 10**6

        if pkd_antagonist == 0:
            kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)
        else:
            kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)*(1+(antagonist_conc/kd_anatagonist))

        #LR determination
        a = 1
        b = float(agonist_conc)+float(receptor_conc)+kd_agonist
        c = float(receptor_conc)*float(agonist_conc)
        delta = (b**2) - (4*a*c)
        LR = (b-math.sqrt(delta))/(2*a)
        return LR

class convert:
    def microgr2nanomolar(uniprotID, concentration):
        from scipy.constants import Avogadro, micro, nano , pico
        from Bio.Seq import Seq
        from Bio.SeqUtils import molecular_weight
        from bioservices import UniProt
        u = UniProt(verbose=False)

        def Da2gr(x):
            return x*1.6605300000013E-24
        
        #Get protein sequence from UniProt
        seq = u.get_fasta_sequence(uniprotID)
        
        #Prot MW
        prot_Da = molecular_weight(seq, "protein") #Da (daltons)
        
        #convert dalton to gram and then to picograms
        prot_gr = Da2gr(prot_Da)

        #convert grams to picograms
        prot_microgr = prot_gr/micro
        prot_microgr

        #calculate protein units (experimental)
        prot_Na = concentration/prot_microgr

        #Molar concentration of prot (experimental)
        prot_M = prot_Na/Avogadro
        prot_nM = prot_M/nano
        
        return round(prot_nM,4)

class get:
    '''tools to retrive protein information'''

    def gprotein(uniprotID):
        import sqlite3
        
        #query HuTRdb
        conn = sqlite3.connect(os.path.join(os.getcwd(),'SSBColab/src/databases/HuTRdb.sqlite3'))
        cur = conn.cursor()
        cur.execute("SELECT * FROM gpcr WHERE uniprotid=?", (uniprotID,))
        rows = cur.fetchall()
        conn.close()

        gprotein=rows[0][7]
        return gprotein

class network:

    class stimulation:

        def run(ligands=None, pkd_values=None, network_name=None, receptor_conc=None, ligand_conc_range=None, ttotal=None, nsteps=None, kinetics=True, kinetic_parameters=None):
            '''
            This function runs the network simulation and gives the raw output from it.

            ligands         list of ligands or receptors (this represent the number of simulations performed)
            pkd_data           data from dlscore (path of csv file)
            network_name       name of the netwokr (eg. G_alpha_s or G_alpha_i)
            lig_conc_range         range of the concentration of ligand.
            ttotal             total time mof simulation
            nsteps             number of time steps/windows'''

            #Check inputs
            if ligands==None: raise TypeError("ligands list undefined.")
            elif pkd_values==None: raise TypeError("pkd_values_dict undefined.")
            elif network_name==None: raise TypeError("network name undefined.")
            elif pkd_values==None: raise TypeError("pkd_values_dict undefined.")
            elif ligand_conc_range.any() == False: raise TypeError("ligand_conc_range undefined.")
            elif ttotal==None: raise TypeError("ttotal undefined.")
            elif nsteps==None: raise TypeError("nsteps undefined.")
            elif receptor_conc==None: raise TypeError("receptor_conc undefined.")
            elif kinetics==True and kinetic_parameters==None: raise TypeError("Kinetic parameters undefined.")
            else: pass


            available_networks = ['Gs', 'Gi', 'Gq']
            if network_name == 'Gz(Gi)': network_name = 'Gi'
            if network_name not in available_networks: raise Exception('Unvailable Network. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')

            #Main function
            pathway = importlib.import_module('.'+network_name, package='src.lib.pathways')
            t = pl.geomspace(0.00001, ttotal, num=nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points

            simulation_data={}
            
            for ligand in ligands:
                sim_data={}
                ligand_name = os.path.splitext(str(ligand))[0]
                data=[]
                sharedfunctions.printProgressBar(0, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                for idx in range(len(ligand_conc_range)):

                    ligand_conc = ligand_conc_range[idx]

                    if kinetics == False:
                        #get LR conc
                        parameters={'R_init':receptor_conc}
                        LR_conc_init = sharedfunctions.LR_eq_conc(receptor_conc, ligand_conc, 0, pkd_values[ligands.index(ligand)], 0)
                        
                        mymodel = pathway.network(LR=LR_conc_init, kinetics=False, **parameters)
                        simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                        yout = simres.all
                    elif kinetics == True:
                        parameters={'R_init':receptor_conc, 'L_init':ligand_conc_range[idx], **kinetic_parameters}
                        mymodel = pathway.network(kinetics=True, **parameters)
                        simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                        yout = simres.all

                    d1={'ligand_conc':ligand_conc, 'time':t }

                    for idx2 in range(len(pathway.list_of_observables)):
                        d2={pathway.list_of_observables[idx2]:yout[pathway.list_of_observables[idx2]]}
                        d1.update(d2)
                    data.append(d1)
                    sharedfunctions.printProgressBar(idx + 1, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                #sim_data['ligand_name'] = ligand_name
                #sim_data['pkd']=str(pkd_data[ligand_name])
                #sim_data['sim_data']=data
                #all_data.append(sim_data)
                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':ligand_name,
                                            'pkd':str(pkd_values[ligands.index(ligand)])}


            return simulation_data

        def analysis(simulation_data, lig_conc_range, pathway):
            '''
            This function runs the network simulation and gives the raw output from it.

            list_names         list od ligands or receptors (this represent the number of simulations performed)
            pkd_data           data from dlscore (path of csv file)
            network_name       name of the netwokr (eg. G_alpha_s or G_alpha_i)
            lig_conc_range         range of the concentration of ligand.
            ttotal             total time mof simulation
            nsteps             number of time steps/windows+
            epsilon            default is 1 (not used in this moment)
            metabolite         name of the metabolite that will be used for dose-response'''
            
            
            from sklearn.preprocessing import minmax_scale
            
            # Define all the lists and dictionaries used in this function
            raw_data=[]
            normalized_data=[]
            fitted_data=[]

            dose={}

            #defining concentration range
            lig_conc_min = lig_conc_range.min()
            lig_conc_max = lig_conc_range.max()
            
            #Main function
            for ligand in simulation_data:

                #definig and dictionaries used in this loop:
                raw_data_dict={}
                normalized_data_dict={}
                fitted_data_dict={}

                # Calculate dose-response curve
                #get metabolite concentration, rescale, and transform data if network/metabolite decrease
                # metabolite_raw is not normalized
                if pathway == 'Gi' or pathway == 'Gz(Gi)':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(lig_conc_range)):
                        n=np.amax(simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(1-np.array(metabolite_conc_raw))
                elif pathway == 'Gs':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(lig_conc_range)):
                        n=np.amax(simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                elif pathway == 'Gq':
                    metabolite='IP3'
                    metabolite_conc_raw=[]
                    for i in range(len(lig_conc_range)):
                        n=np.amax(simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                else: raise Exception('Unvailable Network. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')


                ## save results
                raw_data_dict['x']=lig_conc_range
                raw_data_dict['y']=metabolite_conc_raw
                raw_data_dict['label']=simulation_data[ligand]['label']

                normalized_data_dict['x']=lig_conc_range
                normalized_data_dict['y']=metabolite_conc_norm
                normalized_data_dict['label']=simulation_data[ligand]['label']

                ## create a list of all data
                raw_data.append(raw_data_dict)
                normalized_data.append(normalized_data_dict)

                ##Fitting curve to the data

                def equation_dose(X, Bottom, Top, EC50, p):
                    return Bottom + (Top-Bottom)/(1+np.power((EC50/X),p))

                popt_EC50, pcov = curve_fit(equation_dose, lig_conc_range, metabolite_conc_norm, bounds=([np.min(metabolite_conc_norm),-np.inf,-np.inf, 0.5],[np.inf,np.max(metabolite_conc_norm),np.inf, 2.5]))

                xfit_EC50 = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
                yfit_EC50 = equation_dose(xfit_EC50, *popt_EC50)

                fit_EC50={'x':xfit_EC50, 'y':yfit_EC50, 'label':simulation_data[ligand]['label']}

                dose[ligand] = {'raw_data': raw_data_dict,
                                'normalized_data':normalized_data_dict ,
                                'fitted_data': fit_EC50,
                                'EC50 (μM)': round(popt_EC50[2],5),
                                'pEC50': round(-np.log10(popt_EC50[2]*1E-6),2)}
                


            return dose

        def curve(simulation_data):
            import plotly
            import plotly.graph_objs as go

            colors = plotly.colors.DEFAULT_PLOTLY_COLORS

            plot_data=[]

            color_id=0
            for ligand in simulation_data:
                trace_norm = go.Scatter(x=simulation_data[ligand]['normalized_data']['x'],
                                       y=minmax_scale(simulation_data[ligand]['normalized_data']['y'])*100 ,
                                       mode='markers',
                                       showlegend=True,
                                       name=simulation_data[ligand]['normalized_data']['label'],
                                       marker=dict(color=colors[color_id]))
                plot_data.append(trace_norm)

                trace_fitted = go.Scatter(x=simulation_data[ligand]['fitted_data']['x'],
                                    y=minmax_scale(simulation_data[ligand]['fitted_data']['y'])*100,
                                    mode='lines',
                                    showlegend=False,
                                    name=simulation_data[ligand]['fitted_data']['label'],
                                    line=dict(color=colors[color_id]))
                plot_data.append(trace_fitted)
                color_id +=1

            layout = dict(title = '',
                          xaxis = dict(
                              title = '[ligand] μM',
                              type ='log',
                              range = [-3, 2],
                              titlefont=dict(
                                  size=20
                              ),
                              tickfont=dict(
                                  size=20
                              )),
                          yaxis = dict(
                              title = '% Response',
                              range = [0, 100],
                              titlefont=dict(
                                  size=20),
                              tickfont=dict(
                              size=20)

                          ),
                          legend=dict(font=dict(size=15)),
                          autosize=False,
                          width=850,
                          height=650,
                         )

            fig = go.Figure(data=plot_data, layout=layout)
            return fig

        def constants(simulation_data):
            kvalues={}
            for ligand in simulation_data:
                EC50 = list(simulation_data[ligand].keys())[-2]
                EC50_value = simulation_data[ligand][EC50]
                pEC50 = list(simulation_data[ligand].keys())[-1]
                pEC50_value = simulation_data[ligand][pEC50]
                kvalues[ligand]={EC50:EC50_value, pEC50:pEC50_value}

            return kvalues

