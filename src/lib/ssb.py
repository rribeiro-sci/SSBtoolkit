#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2019 Rui Ribeiro
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

#MODELLER
#from modeller                   import *
#from modeller.optimizers        import molecular_dynamics, conjugate_gradients
#from modeller.automodel         import autosched

#direcotires
from src.lib.directories import *
from src.lib.autogrids.auto_grids import auto_grids
#from Tools.paths import *

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


# Print iterations progress
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


def LR_eq_conc_test(receptor_conc, agonist_conc, antagonist_conc, pkd_agonist, pkd_antagonist):
    #from pKd to Kd
    kd_anatagonist = (10**(-(float(pkd_antagonist)))) * 10**6

    if pkd_antagonist == 0:
        kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)
    else:
        kd_agonist = ((10**(-(float(pkd_agonist)))) * 10**6)*(1+(antagonist_conc/kd_anatagonist))

    #LR determination
    LR = (receptor_conc*agonist_conc)/(kd_agonist+agonist_conc)
    return LR

def equilibrium_binding_curve(receptor_conc, lig_conc_range, antagonist_conc, pkd_agonist, pkd_antagonist):
    curve_data=[]
    for conc in lig_conc_range:
        curve_data.append(LR_eq_conc(receptor_conc, conc, antagonist_conc, pkd_agonist, pkd_antagonist))
    return curve_data

def equilibrium_binding_curve_test(receptor_conc, lig_conc_range, antagonist_conc, pkd_agonist, pkd_antagonist):
    curve_data=[]
    for conc in lig_conc_range:
        curve_data.append(LR_eq_conc_test(receptor_conc, conc, antagonist_conc, pkd_agonist, pkd_antagonist))
    return curve_data

def maxbend(drug_receptor, lig_conc_range):
    from scipy.optimize   import curve_fit, minimize


    def equation_binding(X, Bottom, Top, Kd, p):
        return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

    xfit = np.geomspace(1E-3, 1E4, 50000)
    popt, pcov = curve_fit(equation_binding, lig_conc_range, drug_receptor, bounds=([np.min(drug_receptor),-np.inf,-np.inf, 0.5],[np.inf,np.max(drug_receptor),np.inf, 2.5]))


    def equation_binding_deriv_b(x, a,d,c,b):
        return (x/c)**b*(a - d)*np.log(x/c)/((x/c)**b + 1)**2

    min_value = minimize(equation_binding_deriv_b, np.max(xfit), args=(popt[0],popt[1],popt[2],popt[3]), method = 'Nelder-Mead')

    submaximal = round(min_value.x[0],3)
    return submaximal

def import_dlscore(filename,ligands_names):
    root, ext = os.path.splitext(filename)

    try:
        if 'csv' in ext:
            # Assume that the user uploaded a CSV file
            pkd_list=[]

            with open(filename, "r") as infile:
                next(infile)
                for i in infile:
                    pkd_list.append(i.split(',')[2].rstrip()) # this take the dlscore from csv

            dict_kd = dict(zip(ligands_names, pkd_list))

    except:
        print('There was an error processing this file.')

    return dict_kd


class get:
    '''tools to retrive protein information'''

    def get_seq(structure_filename, chain):
        '''This function extracts the aminoacid sequence of pdb file and return the sequence in one letter code format'''
        #return the one-letter-code sequence of a structure
        #NB: specified chain must be in the FIRST model of the structure

        letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
                   'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
                   'TYR':'Y','VAL':'V'}
        seq = ""

        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure('query', structure_filename)

        chain = structure[0][chain]
        for residue in chain:
            if residue.id[0] == ' ':
                seq = seq+letters[residue.get_resname()]

        return seq

    def get_uniprot_entry_online(structure_filename, chain):
        '''This function returns the Uniprot id of the pdb structure given'''
        seq = get.get_seq(structure_filename, chain)
        results = NCBIWWW.qblast("blastp", "swissprot", seq, format_type="Text",hitlist_size=5 )

        with open("my_blast.txt", "w") as out:
            out.write(results.read())

        results.close()

        with open("my_blast.txt", 'r') as out:
            for i,line in enumerate(out):
                if i == 34:
                    refseq_entry = line.split(sep='.')[0]
                    out.close()
                    break

        os.remove("my_blast.txt")
        return refseq_entry

    def get_uniprot_entry_local(structure_filename, chain, path):
        '''This function returns the Uniprot id of the pdb structure given'''

        seq = get.get_seq(structure_filename, chain)

        fasta_seq = SeqRecord(Seq(seq), id='query', description='')

        SeqIO.write(fasta_seq, path+'seq.fasta', 'fasta')
        blast_command = NcbiblastpCommandline(db=swissprot_db, out=path+'my_blast.txt', query=path+'seq.fasta')
        stdout, stder = blast_command()

        with open(path+"my_blast.txt", 'r') as out:
            for i,line in enumerate(out):
                if i == 29:
                    refseq_entry = line.split(sep='.')[0]
                    out.close()
                    break

        os.remove(path+"my_blast.txt")
        os.remove(path+'seq.fasta')
        return refseq_entry

    def get_geneid(uniprot_entry):
        '''uniprot id as input and give us the gene id/name of the protein'''
        uniprot_url_query = "https://www.ebi.ac.uk/proteins/api/proteins/"+uniprot_entry
        html_data = urllib.urlopen(uniprot_url_query)
        datas = json.loads(html_data.read())

        organism = datas["organism"]["names"][0]["value"]
        if organism != "Homo sapiens":
            print("\nWARNING:   This is not a Human gene!\n")

        gene_id = datas["gene"][0]["name"]["value"]

        return gene_id

    def get_uniprot_seq(uniprot_entry):
        '''uniport id as input and give us the original uniprot sequence of protein (download FASTA but return as string)'''
        uniprot_url = "https://www.uniprot.org/uniprot/"
        query_url = uniprot_entry+".fasta"
        html_data = urllib.urlopen(uniprot_url+query_url)

        seq_fasta = html_data.read().decode("utf-8")

        #transform fasta sequence into straight sequence string
        fasta_splitted = seq_fasta.split(sep='\n')
        seq_str = ''.join(fasta_splitted[1:])

        return seq_str

    def search_approved_drugs(uniprot_entry, ligands_directory):
        try:
            mychem_url_query = "http://mychem.info/v1/query?q=drugbank.targets.uniprot:"+uniprot_entry+"%20AND%20drugbank.groups:approved&fields=drugbank.id,drugbank.name"
            html_data = urllib.urlopen(mychem_url_query)
            datas = json.loads(html_data.read())
            list_of_ligands=[]
            #retrive SDF structures from DrugBank
            for hit in datas["hits"]:

                DrugBank_accession = hit["drugbank"]["id"]
                ligand_name = hit["drugbank"]["name"]
                drugbank_ligand_structure_url = "https://www.drugbank.ca/structures/small_molecule_drugs/"+DrugBank_accession+".sdf?type=3d"
                html_data = urllib.urlopen(drugbank_ligand_structure_url)

                #check if there are spaces in the name of the ligand
                if len(ligand_name.split()) > 1:
                    ligand_name = ''.join(ligand_name.split())

                with open(ligands_directory+'/'+ligand_name+'.sdf', "w") as f:
                    f.write(html_data.read().decode("utf-8"))

                list_of_ligands.append(ligand_name)
            #if no drug-bank approved drugs were founded
            # if not datas["hits"]:
            #
            #     print("\nWARNIG: NO APPROVED DRUGS FOUNDED\n")

        except:
            print("\nWARNING: SOME STRUCTURES OF APPROVED DRUGS NOT FOUND!\n")
        return list_of_ligands

    def search_exac_variants(gene_id):
        '''gene id as input. we can call the function retrive_geneid to get the gene id. Give us a list of all of annotated missense variant
        for the gene.    p.Arg300Lys'''
        #retrive all the missense variations annotated in the ExAc browser for the passed gene

        exac_url = "http://exac.hms.harvard.edu/"
        api_query_url = "/rest/awesome?query="+gene_id+"&service=variants_in_gene"

        html_data_query = urllib.urlopen(exac_url+api_query_url)

        variants = json.loads(html_data_query.read())
        variants_list = []

        for variant in variants:
            if variant["category"] == "missense_variant" and variant["major_consequence"] == "missense_variant":
                variants_list.append(variant["HGVSp"])


        return variants_list

    def filter_bindingsite_varitants(receptor_filename, chain, family, variants_list, max_distance, path):
        '''pdb file as input + chain of our subunit + family + python list of HGVSp variants (from previous function).
        give us a python list of all the input variants che stano nel biding site.

        family = GPCR or IGR (ionotropic glutamatergic receptor)

        chain in capital'''
        #given a receptor leading to the specified family and a list of variants (in the HGVSp notation)
        #annotated for the specified chain, return only the variants wich lay on the Ligand Binding Domain of the passed protein's structure.

        #N.B. : a preventive renumeration of the structure is recommended

        def evaluate_distance(center_coordinates, atoms_coordinates, max_distance):

            #if the distance from the passed center of (at least) one of the passed atoms is lower than max_distance return true, false otherwise
            for atom_coordinates in atoms_coordinates:

                d = distance.euclidean(atom_coordinates, center_coordinates)
                if d < max_distance:
                    return True

            return False

        #find the center of the Binding Pocket with bp_center
        center_coordinates = bp_center.bp_center(path+receptor_filename, chain, family)
        #center_coordinates = [32.03, 8.71, -16.1]
        #print(center_coordinates)
        #set maximum allowed distance (Amstrong)
        #max_distance = 22

        #parse the structure
        parser = PDBParser(PERMISSIVE=1)
        structure =  parser.get_structure("", path+receptor_filename)

        #NB: passed chain must be in the first model of the structure
        chain = structure[0][chain]
        structure_res_numbers = []
        filtered_variants_list = []

        for residue in chain:
            structure_res_numbers.append(list(residue.id)[1])


        for variant in variants_list:
            variant_res_number = int(variant[5:len(variant)-3])

            if variant_res_number in structure_res_numbers:
                atoms_coordinates = []

                for atom in chain[variant_res_number]:
                    atoms_coordinates.append(atom.get_coord())

                if evaluate_distance(center_coordinates, atoms_coordinates, max_distance):
                    filtered_variants_list.append(variant)

        return filtered_variants_list

    def get_prot_info(receptor, chain, dirpath):
        import sqlite3
        #get uniprotid in local
        uniprotid = get.get_uniprot_entry_local(os.path.join(dirpath, receptor), chain, dirpath)
        #query our db
        conn = sqlite3.connect(HuTRdb_path)
        cur = conn.cursor()
        cur.execute("SELECT * FROM gpcr WHERE uniprotid=?", (uniprotid,))
        rows = cur.fetchall()
        conn.close()

        info={'class':rows[0][1], 'class':rows[0][2], 'type':rows[0][3],'name':rows[0][4],'uniprotid':rows[0][5],'geneid':rows[0][6],'gprotein':rows[0][7]}
        return info

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

class handle_pdb:

    def clean_pdb(file_pdb, out_pdb):

        with open(file_pdb, "r") as input:
            lines = []
            for line in input:
                columns = line.split()
                if columns[0] in ['ATOM', 'TER']:
                    lines.append(line)

        with open("temp.tmp", "w") as output:
            output.write("MODEL1\n")
            output.writelines(lines)
            output.write("ENDMDL\n")
            output.write("END")

        if file_pdb == out_pdb:
            os.remove(file_pdb)

        os.rename("temp.tmp", out_pdb)

        return out_pdb

    def pdbqt2pdb(file_pdbqt, file_pdb):

        #convert pdbqt to pdb
        #subprocess.call('babel -i pdbqt '+file_pdbqt+' -o pdb '+file_pdb)
        os.system('babel -i pdbqt '+file_pdbqt+' -o pdb '+file_pdb)
        print('Converting '+file_pdbqt+' to '+file_pdb)
        return file_pdb

    def append_pdb(file1_pdb, file2_pdb, out_filename):

        #append file1_pdb to file2_pdb
        ref = open(file2_pdb,'r')
        app = open(file1_pdb,'r')
        out = open(out_filename,'w+')

        for line in ref:
            columns=line.split()
            if columns[0] in ['ATOM','HETATM']:
                out.write(line)

        ref.close()

        for line in app:
            columns = line.split()
            if columns[0] in ['ATOM','CONNECT', 'HETATM']:
                out.write(line)

        app.close()

        out.write('END')
        out.close()
        return out_filename

    def merge_pdbqts(pdbqt_filename1, pdbqt_filename2, pdb_out_filename):

        #merge two pdbqt files in a pdb file
        pdbqt_basename1 = (os.path.splitext(os.path.basename(pdbqt_filename1)))[0]
        pdbqt_basename2 = (os.path.splitext(os.path.basename(pdbqt_filename2)))[0]
        converted_filename1 = pdbqt_basename1+'.pdb'
        converted_filename2 = pdbqt_basename2+'.pdb'
        pdbqt2pdb(pdbqt_filename1, converted_filename1)
        pdbqt2pdb(pdbqt_filename2, converted_filename2)
        append_pdb(converted_filename1, converted_filename2, out_filename)

        os.remove(converted_filename1)
        os.remove(converted_filename2)

    def vina_best_result(pdbqt_result_filename, pdbqt_best_result_filename):

        #extract the first results from vina pdbqt result file
        with open(pdbqt_result_filename, "r") as file:
            lines = []
            for line in file:
                lines.append(line)
                columns=line.split()
                if columns[0] == 'ENDMDL': break

        with open(pdbqt_best_result_filename, "w") as out:
            out.writelines(lines)

        return pdbqt_best_result_filename

    def renumerate_structure(structure_pdb, chain, uniprot_entry, structure_pdb_out):

        #Renumerate the residues of the passed structure's chain according with the original uniprot sequence
        #Sequences must be passed as string objects
        #NB: specified chain must be in the FIRST model of the structure

        def recursive_renumbering(new_id, new_ids, chain, dict):

            residue = chain[new_id]
            current_id = list(residue.id)
            index = dict[current_id[1]]
            new_id = new_ids[index]
            if int(current_id[1]) > int(new_id) and new_id in chain:
                recursive_renumbering(new_id, new_ids, chain, dict)

            current_id[1] = new_id
            residue.id = tuple(current_id)

        structure_seq = get.get_seq(structure_pdb, chain)
        original_seq_fasta_file = get.get_uniprot_seq(uniprot_entry)

        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure('query', structure_pdb)

        alignment = pairwise2.align.localdd(original_seq_fasta_file, structure_seq, blosum62, -1000, -1000, -10, -0.5)[0]

        new_ids = []
        for i in range(len(alignment[1])):

            if alignment[1][i] != '-':
                new_ids.append(i+1)

        chain = structure[0][chain]

        #create a residue : structure-relative-position dict
        dict = {}
        for index, residue in enumerate(chain):
            dict.update({residue.id[1] : index})

        #start the re-enumeration from the end of the structure
        for index, residue in reversed(list(enumerate(chain))):

            current_id = list(residue.id)
            new_id = new_ids[index]

            #check if we're trying to use an already exixting residue identifier, if so, first re-enumerate that residue
            if int(current_id[1]) > int(new_id) and new_id in chain:
                recursive_renumbering(new_id, new_ids, chain, dict)

            current_id[1] = new_id
            residue.id = tuple(current_id)


        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(structure_pdb_out)

        return structure_pdb_out

class modelling:
    def mutate(structure_pdb, chain, variant, out_pdb, path):
        '''pdb file as input + chain + variant in HGVSp notation + output name (inlcuding pdb extension)'''
        #perform the computation af a mutatant model for the passed structure
        #the mutation will be performed over the specified chain, according to the passed HGVSp variant
        #NB : variant must be passed in the HGVSp notation
        '''modeller script'''
        def optimize(atmsel, sched):
            #conjugate gradient
            for step in sched:
                step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
            #md
            refine(atmsel)
            cg = conjugate_gradients()
            cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)

        #molecular dynamics
        def refine(atmsel):
            # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
            md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                                    md_return='FINAL')
            init_vel = True
            for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                        (200, 600,
                                         (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
                for temp in temps:
                    md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                                 max_iterations=its, equilibrate=equil)
                    init_vel = False

        #use homologs and dihedral library for dihedral angle restraints
        def make_restraints(mdl1, aln):
           rsr = mdl1.restraints
           rsr.clear()
           s = selection(mdl1)
           for typ in ('stereo', 'phi-psi_binormal'):
               rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
           for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
               rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                        spline_dx=0.3, spline_min_points = 5, aln=aln,
                        spline_on_site=True)

        residues = {'Ala':'ALA','Arg':'ARG','Asn':'ASN','Asp':'ASP','Cys':'CYS','Glu':'GLU','Gln':'GLN','Gly':'GLY','His':'HIS',
                   'Ile':'ILE','Leu':'LEU','Lys':'LYS','Met':'MET','Phe':'PHE','Pro':'PRO','Ser':'SER','Thr':'THR','Trp':'TRP',
                   'Tyr':'TYR','Val':'VAL'}

        # Set a different value for rand_seed to get a different final model
        env = environ(rand_seed=-49837)

        env.io.hetatm = True
        #soft sphere potential
        env.edat.dynamic_sphere=False
        #lennard-jones potential (more accurate)
        env.edat.dynamic_lennard=True
        env.edat.contact_shell = 4.0
        env.edat.update_dynamic = 0.39

        # Read customized topology file with phosphoserines (or standard one)
        env.libs.topology.read(file='$(LIB)/top_heav.lib')

        # Read customized CHARMM parameter library with phosphoserines (or standard one)
        env.libs.parameters.read(file='$(LIB)/par.lib')

        # Read the original PDB file and copy its sequence to the alignment array:
        mdl1 = model(env, file=structure_pdb)
        ali = alignment(env)
        ali.append_model(mdl1, atom_files=structure_pdb, align_codes=structure_pdb)

        '''pietro start'''
        #extract the position and the substitute residue from the passed HGVSp-notated variant
        respos = variant[5:len(variant)-3]
        restyp = residues[variant[len(variant)-3:]]
        '''pietro stop'''

        #set up the mutate residue selection segment
        s = selection(mdl1.chains[chain].residues[respos])

        #perform the mutate residue operation
        s.mutate(residue_type=restyp)
        #get two copies of the sequence.  A modeller trick to get things set up
        ali.append_model(mdl1, align_codes=structure_pdb)

        # Generate molecular topology for mutant
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])

        # Transfer all the coordinates you can from the template native structure
        # to the mutant (this works even if the order of atoms in the native PDB
        # file is not standard):
        #here we are generating the model by reading the template coordinates
        mdl1.transfer_xyz(ali)

        # Build the remaining unknown coordinates
        mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

        #yes model2 is the same file as model1.  It's a modeller trick.
        mdl2 = model(env, file=structure_pdb)

        #required to do a transfer_res_numb
        #ali.append_model(mdl2, atom_files=structure_pdb, align_codes=structure_pdb)
        #transfers from "model 2" to "model 1"
        mdl1.res_num_from(mdl2,ali)

        #It is usually necessary to write the mutated sequence out and read it in
        #before proceeding, because not all sequence related information about MODEL
        #is changed by this command (e.g., internal coordinates, charges, and atom
        #types and radii are not updated).

        mdl1.write(file='temp.tmp')
        mdl1.read(file='temp.tmp')

        #set up restraints before computing energy
        #we do this a second time because the model has been written out and read in,
        #clearing the previously set restraints
        make_restraints(mdl1, ali)

        #a non-bonded pair has to have at least as many selected atoms
        mdl1.env.edat.nonbonded_sel_atoms=1

        sched = autosched.loop.make_for_model(mdl1)

        #only optimize the selected residue (in first pass, just atoms in selected
        #residue, in second pass, include nonbonded neighboring atoms)
        #set up the mutate residue selection segment
        s = selection(mdl1.chains[chain].residues[respos])

        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)

        s.energy()

        s.randomize_xyz(deviation=4.0)

        mdl1.env.edat.nonbonded_sel_atoms=2
        optimize(s, sched)

        #feels environment (energy computed on pairs that have at least one member
        #in the selected)
        mdl1.env.edat.nonbonded_sel_atoms=1
        optimize(s, sched)

        s.energy()

        #give a proper name
        mdl1.write(file=path+out_pdb)

        #delete the temporary file
        os.remove('temp.tmp')

        return out_pdb

    def variants(receptor, list_of_variants, receptor_basename, path, modeller_pwd):
        models=[]
        if modeller_pwd == 'MODELIRANJE':
            for variant in list_of_variants:
                out = receptor_basename+'_'+variant.split('.')[1]+'.pdb'
                mut = modelling.mutate(path+receptor, 'A', variant, out, path)
                handle_pdb.clean_pdb(path+out, path+out)
                models.append(out)
        else: raise Exception('Wrong MODELLER license key. Please obtain a license key from https://salilab.org/modeller/registration.html') 
        
        return models

class docking:

    def grid(receptor_name, family):
        '''perform an automatic detection of the binding pocket center coordinates, for specific proteins families (also works for heteromultimer)'''
        coordinates_dict = auto_grids(receptor_name, family, reference_files_dir)
        return coordinates_dict

    def vina(receptors, ligands, path, family, size_x, size_y, size_z, exhaustiveness, energy_range, num_modes, ncpu):
        import ast
        '''This function execute AUTODOCK VINA for virtual screening. It receives as input
        a list of receptors' names, a list of ligands' names, and a path where the files a located.
        It is only working with pdb files as input because the pdb extension if automatically added.
        The function prepares the receptors and ligands  according the default scripts of MGLTOOLS.
        The output is a pandas dataframe of vina results for each ligand docked, organized in a Pandas Panel.'''
        print('Calculating grid...\n')
        mygrid = auto_grids(path+receptors[0], family, reference_files_dir)
        print('Vina is running...\n')
        # Initial call to print 0% progress
        printProgressBar(0, len(ligands), prefix = 'Progress:', suffix = 'Complete', length = 50)
        i=0
        for receptor in receptors:
            receptor_name, receptor_ext = os.path.splitext(receptor)
            #print('\nLoading '+receptor_name+'...')
            prepare_receptor_command = prepare_receptor + ' -r '+path+receptor+' -o '+path+receptor_name+'.pdbqt'+' -A hydrogens'
            #print(prepare_receptor_command) #for debugging
            os.system(prepare_receptor_command)
            ii=0
            for ligand in ligands:
                ligand_name, ligand_ext = os.path.splitext(ligand)
                #print('    Docking '+ligand_name+'...')
                prepare_ligand_command = prepare_ligand + ' -l '+path+ligand+' -o '+path+ligand_name+'.pdbqt'
                #print(prepare_ligand_command) #for debugging
                os.system(prepare_ligand_command)
                vina_command = vina+' --receptor '+path+receptor_name+'.pdbqt --ligand '+path+ligand_name+'.pdbqt --center_x '+str(mygrid['A'][0])+' --center_y '+str(mygrid['A'][1])+' --center_z '+str(mygrid['A'][2])+ ' --size_x '+str(size_x)+ ' --size_y '+str(size_y)+ ' --size_z '+str(size_z)+ ' --exhaustiveness '+str(exhaustiveness)+ ' --energy_range '+str(energy_range)+ ' --num_modes '+str(num_modes)+ ' --cpu '+str(ncpu)+ ' --out '+path+receptor_name+'_'+ligand_name+'_vina_out.pdbqt --log '+path+receptor_name+'_'+ligand_name+'_vina_out.log'
                #print(vina_command) #for debugging
                os.system(vina_command)
                # Update Progress Bar
                if len(receptors)<2:
                  printProgressBar(ii + 1, len(ligands), prefix = 'Progress:', suffix = 'Complete', length = 50)
                  ii +=1
            
            # Update Progress Bar
            if len(receptors)>1:
                
                printProgressBar(i + 1, len(receptors), prefix = 'Progress:', suffix = 'Complete', length = 50)
                i +=1
        
        print('\nDone!')
        #########################################

        vina_dict={}
        vina_files=glob(path+'*vina_out.pdbqt')

        for docking in vina_files:
            root, ext = os.path.splitext(docking)
            #print(root)
            head, tail =os.path.split(root)
            #parse vina output
            with open(root+'.log') as vina_log:
                lines = vina_log.read()
                lines_splitted = lines.split('\n')
            vina_log_data = []
            for line in lines_splitted[24:-2]:
                mode = line.strip().split()
                vina_log_data.append(mode)
            vina_log_df = pd.DataFrame(vina_log_data, columns=['mode','affinity (kcal/mol)', 'rmsd l.b.','rmsd u.b.'])
            vina_log_df['mode'] = pd.to_numeric(vina_log_df['mode'])

            results = vina_log_df[vina_log_df['mode'] <= 10]
            vina_dict[tail] = ast.literal_eval(results.to_json(orient='records'))
        #########################################
        #     #Split modes
        #     #file_pdbqt = path+ligand_name+'_out.pdbqt'
        #     #split = vina_split+' --input '+docking+' --ligand '+path+tail
        #     #print(split)
        #     #os.system(split)

            handle_pdb.vina_best_result(docking, root+'_best.pdbqt')

        return vina_dict

    def vina_flex(receptors, ligands, flex_residues, path, family, size_x, size_y, size_z, exhaustiveness, energy_range, num_modes, ncpu):
        import ast
        '''This function execute AUTODOCK VINA for virtual screening. It receives as input
        a list of receptors' names, a list of ligands' names, and a path where the files a located.
        It is only working with pdb files as input because the pdb extension if automatically added.
        The function prepares the receptors and ligands  according the default scripts of MGLTOOLS.
        The output is a pandas dataframe of vina results for each ligand docked, organized in a Pandas Panel.'''
        print('Calculating grid...\n')
        mygrid = auto_grids(path+receptors[0], family, reference_files_dir)
        print('Vina is running...\n')
        # Initial call to print 0% progress
        printProgressBar(0, len(ligands), prefix = 'Progress:', suffix = 'Complete', length = 50)
        i=0
        for receptor in receptors:
            receptor_name, receptor_ext = os.path.splitext(receptor)
            #print('\nLoading '+receptor_name+'...')
            prepare_receptor_command = prepare_receptor + ' -r '+path+receptor+' -o '+path+receptor_name+'.pdbqt -A hydrogens'
            os.system(prepare_receptor_command)
            #print(prepare_receptor_command) #for debugging

            prepare_receptor_flex_command = prepare_receptor_flex + ' -r '+path+receptor_name+'.pdbqt -s '+flex_residues+' -g '+path+receptor_name+'_rigid.pdbqt -x '+path+receptor_name+'_flex.pdbqt > '+path+'error2.log'
            #print(prepare_receptor_flex_command) #for debugging
            os.system(prepare_receptor_flex_command)
            
            
            ii=0
            for ligand in ligands:
                ligand_name, ligand_ext = os.path.splitext(ligand)
                #print('    Docking '+ligand_name+'...')
                prepare_ligand_command = prepare_ligand + ' -l '+path+ligand+' -o '+path+ligand_name+'.pdbqt '
                #print(prepare_ligand_command) #for debugging
                os.system(prepare_ligand_command)
                vina_command = vina+' --receptor '+path+receptor_name+'_rigid.pdbqt --flex '+path+receptor_name+'_flex.pdbqt --ligand '+path+ligand_name+'.pdbqt --center_x '+str(mygrid['A'][0])+' --center_y '+str(mygrid['A'][1])+' --center_z '+str(mygrid['A'][2])+ ' --size_x '+str(size_x)+ ' --size_y '+str(size_y)+ ' --size_z '+str(size_z)+ ' --exhaustiveness '+str(exhaustiveness)+ ' --energy_range '+str(energy_range)+ ' --num_modes '+str(num_modes)+ ' --cpu '+str(ncpu)+ ' --out '+path+receptor_name+'_'+ligand_name+'_vina_out.pdbqt --log '+path+receptor_name+'_'+ligand_name+'_vina_out.log'
                #print(vina_command) #for debugging
                os.system(vina_command)
                # Update Progress Bar
                if len(receptors)<2:
                    printProgressBar(ii + 1, len(ligands), prefix = 'Progress:', suffix = 'Complete', length = 50)
                    ii +=1
            
            # Update Progress Bar
            if len(receptors)>1:
                
                printProgressBar(i + 1, len(receptors), prefix = 'Progress:', suffix = 'Complete', length = 50)
                i +=1
        
        print('\nDone!')
        #########################################

        vina_dict={}
        vina_files=glob(path+'*vina_out.pdbqt')

        for docking in vina_files:
            root, ext = os.path.splitext(docking)
            #print(root)
            head, tail =os.path.split(root)
            #parse vina output
            with open(root+'.log') as vina_log:
                lines = vina_log.read()
                lines_splitted = lines.split('\n')
            vina_log_data = []
            for line in lines_splitted[24:-2]:
                mode = line.strip().split()
                vina_log_data.append(mode)
            vina_log_df = pd.DataFrame(vina_log_data, columns=['mode','affinity (kcal/mol)', 'rmsd l.b.','rmsd u.b.'])
            vina_log_df['mode'] = pd.to_numeric(vina_log_df['mode'])

            results = vina_log_df[vina_log_df['mode'] <= 10]
            vina_dict[tail] = ast.literal_eval(results.to_json(orient='records'))
        #########################################
        #     #Split modes
        #     #file_pdbqt = path+ligand_name+'_out.pdbqt'
        #     #split = vina_split+' --input '+docking+' --ligand '+path+tail
        #     #print(split)
        #     #os.system(split)

            handle_pdb.vina_best_result(docking, root+'_best.pdbqt')

        return vina_dict

    def prepare_dlscore(ligands,docked_files, path):

        file = open(path+'dlscore_input.pdbqt','w')
        number = 1
        for ligand in ligands:
            ligand_name, ligand_extension = os.path.splitext(ligand)
            for i in range(len(docked_files)):
                if ligand_name in docked_files[i]:
                    modelfile=open(docked_files[i])
                    lines=modelfile.readlines()
                    modelfile.close()
                    firt_line='MODEL '+str(number)+'\n'
                    last_line='ENDMDL\n'
                    file.write(firt_line)
                    for n in lines:
                        if n.startswith('MODEL 1'):pass
                        else: file.write(n)
                    file.write(last_line)
                    number +=1
        file.close()
        dlscore_file = path+'dlscore_input.pdbqt'
        return dlscore_file

    def prepare_dlscore_input(docked_files, path):
        dlscore_input_file = open(path+'dlscore_input.pdbqt', 'w')
        number = 1
        for file in docked_files:
            modelfile = open(file)
            lines=modelfile.readlines()
            modelfile.close()
            first_line='MODEL '+str(number)+'\n'
            last_line='ENDMDL\n'
            dlscore_input_file.write(first_line)
            for line in lines:
                if line.startswith('MODEL 1'): pass
                else: dlscore_input_file.write(line)
            dlscore_input_file.write(last_line)
            number +=1
        return dlscore_input_file.name

    def dlscore_local(receptor_file, dlscore_input_file, path):
        '''DLSCORE: A deep learning based scoring function for predicting protein-ligand binding affinity
            DLSCORE recieves as input the docking structure but in different files for ligands and receptor.'''


        receptor_name, receptor_ext = os.path.splitext(receptor_file)

        #Defining dlscore command
        print('DLSCORE is running...')
        dlscore_command='python3 '+dlscore_path+'dlscore.py --receptor '+path+receptor_name+'.pdbqt --ligand '+dlscore_input_file+' --network_type refined --vina_executable '+dlscore_path+'autodock_vina_1_1_2/bin/vina --output ' +path+'dlscore_raw_output --verbose 1 > ' +path+'dlscore_raw_output.log'
        #print (dlscore_command) #for debbuging

        #Running DLSCORE
        os.system(dlscore_command)

        #dlscore_output_path = job_directory+job_name+'.csv'
        out_file=path+'dlscore_raw_output.csv'
        #print('DLSCORE is done in '+str(end-start)+' seconds! Results saved in '+job_directory+job_name+'.csv')

        print('Done!')
        return out_file

    def dlscore(receptors, ligands, path, network_type,num_networks):
        results={}
        
        #we need to diferenciate multiple receptors and multiple ligands. this is important for naming
        if len(receptors) >1: 
            dlscore_rtype = 'multiplereceptor'
        else: 
            dlscore_rtype = "singlereceptor" 
        if len(ligands)>1: 
            dlscore_ltype = 'multipleligand'
        else: 
            dlscore_ltype = 'singleligand'
    
        for receptor in receptors:
            receptor_name = os.path.splitext(receptor)[0]
            receptor_file = path+receptor_name+'.pdbqt'
            for ligand in ligands:
                ligand_name = os.path.splitext(ligand)[0]
                docked_ligand_file= path+receptor_name+'_'+ligand_name+'_vina_out_best.pdbqt'
                dlscore_out = path+receptor_name+'_'+ligand_name+'_dlscore_out'
                docking.run_dlscore(receptor_file, docked_ligand_file, dlscore_out, path, network_type, num_networks)

                #naming
                if dlscore_rtype == 'singlereceptor' and dlscore_ltype == 'singleligand': 
                    dlscore_name = ligand_name
                elif dlscore_rtype == 'singlereceptor' and dlscore_ltype == 'multipleligand': 
                    dlscore_name = ligand_name
                elif dlscore_rtype == 'multiplereceptor' and dlscore_ltype == 'singleligand':
                    dlscore_name = receptor_name
                elif dlscore_rtype == 'multiplereceptor' and dlscore_ltype == 'multipleligand':
                    dlscore_name = receptor_name+'_'+ligand_name
                    
                    
                    
                with open(dlscore_out+'.csv', 'r') as f:
                    headline = f.readline()
                    score_values = f.readline().rstrip().split(',')
                    results[dlscore_name] = {'vina_score':round(float(score_values[0]),3),
                                                            'nnscore':round(float(score_values[1]),3),
                                                            'dlscore': round(float(score_values[2]),3)}

        # Serializing json
        json_object = json.dumps(results, indent = 4)

        # Writing to sample.json
        with open(os.path.join(path,"dlscore_output.json"), "w") as outfile:
            outfile.write(json_object)
        return results

    ######
    def receptor2pdbqt(receptors_list, path):
        for receptor in receptors_list:
            receptor_name, receptor_ext = os.path.splitext(receptor)
            prepare_receptor_command = prepare_receptor + ' -r '+path+receptor+' -o '+path+receptor_name+'.pdbqt'+' -A hydrogens'
            os.system(prepare_receptor_command)
        return

    def ligand2pdbqt(ligands_list, path):
        for ligand in ligands_list:
            ligand_name, ligand_ext = os.path.splitext(ligand)
            prepare_ligand_command = prepare_ligand + ' -l '+path+ligand+' -o '+path+ligand_name+'.pdbqt'
            os.system(prepare_ligand_command)
        return
            

    def vina_HPC(receptors, ligands, path, family, size_x, size_y, size_z, exhaustiveness, energy_range, num_modes, ncpu):
        import ast
        '''This function execute AUTODOCK VINA for virtual screening. It receives as input
        a list of receptors' names, a list of ligands' names, and a path where the files a located.
        It is only working with pdb files as input because the pdb extension if automatically added.
        The function prepares the receptors and ligands  according the default scripts of MGLTOOLS.
        The output is a pandas dataframe of vina results for each ligand docked, organized in a Pandas Panel.'''
        print('Calculating grid...\n')
        mygrid = auto_grids(path+receptors[0], family, reference_files_dir)
        print('Vina is running...\n')
        
        for receptor in receptors:
            receptor_name, receptor_ext = os.path.splitext(receptor)
            #print('\nLoading '+receptor_name+'...')
            prepare_receptor_command = prepare_receptor + ' -r '+path+receptor+' -o '+path+receptor_name+'.pdbqt'+' -A hydrogens'
            #print(prepare_receptor_command) #for debugging
            os.system(prepare_receptor_command)
            
            for ligand in ligands:
                ligand_name, ligand_ext = os.path.splitext(ligand)
                #print('    Docking '+ligand_name+'...')
                prepare_ligand_command = prepare_ligand + ' -l '+path+ligand+' -o '+path+ligand_name+'.pdbqt'
                #print(prepare_ligand_command) #for debugging
                os.system(prepare_ligand_command)
                vina_command = vina+' --receptor '+path+receptor_name+'.pdbqt --ligand '+path+ligand_name+'.pdbqt --center_x '+str(mygrid['A'][0])+' --center_y '+str(mygrid['A'][1])+' --center_z '+str(mygrid['A'][2])+ ' --size_x '+str(size_x)+ ' --size_y '+str(size_y)+ ' --size_z '+str(size_z)+ ' --exhaustiveness '+str(exhaustiveness)+ ' --energy_range '+str(energy_range)+ ' --num_modes '+str(num_modes)+ ' --cpu '+str(ncpu)+ ' --out '+path+receptor_name+'_'+ligand_name+'_vina_out.pdbqt --log '+path+receptor_name+'_'+ligand_name+'_vina_out.log'
                #print(vina_command) #for debugging
                os.system(vina_command)
        
        print('\nDone!')

    def parse_vina(path):
        import ast
        vina_dict={}
        vina_files=glob(path+'*vina_out.pdbqt')

        for docking in vina_files:
            root, ext = os.path.splitext(docking)
            #print(root)
            head, tail =os.path.split(root)
            #parse vina output
            with open(root+'.log') as vina_log:
                lines = vina_log.read()
                lines_splitted = lines.split('\n')
            vina_log_data = []
            for line in lines_splitted[24:-2]:
                mode = line.strip().split()
                vina_log_data.append(mode)
            vina_log_df = pd.DataFrame(vina_log_data, columns=['mode','affinity (kcal/mol)', 'rmsd l.b.','rmsd u.b.'])
            vina_log_df['mode'] = pd.to_numeric(vina_log_df['mode'])

            results = vina_log_df[vina_log_df['mode'] <= 10]
            vina_dict[tail] = ast.literal_eval(results.to_json(orient='records'))
        #########################################
        #     #Split modes
        #     #file_pdbqt = path+ligand_name+'_out.pdbqt'
        #     #split = vina_split+' --input '+docking+' --ligand '+path+tail
        #     #print(split)
        #     #os.system(split)

            handle_pdb.vina_best_result(docking, root+'_best.pdbqt')

        return vina_dict


    ######

    def parse_dlscore_input(dlscore_out, docked_files, ):
        '''This function will match the ligands to each row in the dlscore output'''
        results={}
        ligand_names=[]
        for file in docked_files:
            ligand_name = os.path.splitext(os.path.basename(file))[0].split('_vina_out_best')[0].split('_')[-1]
            ligand_names.append(ligand_name)

        dlscore_lines=[]
        with open(dlscore_out, 'r') as f:
            headline = f.readline()
            for line in f:
                dlscore_lines.append(line.rstrip().split(','))

        for i in range(len(ligand_names)):
            vina_score=(round(float(dlscore_lines[i][0]),3))
            vina_score=(round(float(dlscore_lines[i][0]),3))
            vina_score=(round(float(dlscore_lines[i][0]),3))

            results[ligand_names[i]] = {'vina_score':round(float(dlscore_lines[i][0]),3),
                                        'nnscore':round(float(dlscore_lines[i][1]),3),
                                       'dlscore': round(float(dlscore_lines[i][2]),3)}

        # Serializing json
        json_object = json.dumps(results, indent = 4)

        # Writing to sample.json
        with open(os.path.join(os.path.dirname(dlscore_out),"dlscore_output.json"), "w") as outfile:
            outfile.write(json_object)
        return results

    def run_dlscore(receptor_file, docked_ligand_file, dlscore_out, path, network_type, num_networks):
        #Defining dlscore command
        dlscore_command='python3 '+dlscore_path+'dlscore.py --receptor '+receptor_file+' --ligand '+docked_ligand_file+' --network_type '+network_type+' --num_networks '+num_networks+' --vina_executable '+vina+' --output ' +dlscore_out+' --verbose 1 > ' +dlscore_out+'.log'
        #Running DLSCORE
        os.system(dlscore_command)
        return

#change class name (danger!!)
class bp_center:
    '''trova il centro del ligand binding domain'''
    #reference_files = bp_center_reference_files

    def get_seq(res_list):

        #return the one-letter-code sequence of a list of residues
        letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
                   'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
                   'TYR':'Y','VAL':'V'}
        residues = []
        seq = ""

        for res in res_list:
            if res.id[0] == ' ':
                seq = seq + letters[res.get_resname()]

        return seq

    def superimpose_get_rotranrms(fixed, moving):

        # 'fixed' and 'moving' are lists of Atom objects.
        # The moving atoms will be put on the fixed atoms.
        # return the transformation matrices

        sup = Superimposer()
        sup.set_atoms(fixed, moving)

        # calculate rot and tran matrices.

        rot,tran = sup.rotran
        rms = sup.rms

        return rms, rot, tran

    def allign_3D(ref_list_atoms, ref_list_res, targ_list_res):

        seq1 = bp_center.get_seq(targ_list_res)
        seq2 = bp_center.get_seq(ref_list_res)

        alignment = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -.5)[0]
        score = alignment[2]

        #set the two atoms lists to use to find the transformation matrix
        ref_list = []
        targ_list = []

        #create a list of indices iterating over alignment, where each element contains a couple of value
        #(start and stop of the continuous match in the alignment)

        start_tar = 0
        stop_tar = 0
        start_ref = 0
        stop_ref = 0
        indexes_tar =[]
        indexes_ref = []
        continuous = False
        count_tar = 0
        count_ref = 0

        max = 0

        for i in range(len(alignment[1])):

            if (alignment[0])[i] != '-' :
                count_tar += 1

            if (alignment[1])[i] != '-' :
                count_ref += 1

            if (alignment[0])[i] != '-' and (alignment[1])[i] != '-' :
                if continuous == False:
                    start_tar = count_tar
                    start_ref = count_ref
                    continuous = True
            else:
                if continuous == True:
                    stop_tar = count_tar
                    indexes_tar.append([start_tar,stop_tar])
                    stop_ref = count_ref
                    indexes_ref.append([start_ref,stop_ref])
                    continuous = False
                    if (stop_tar - start_tar) > max :
                        max = stop_tar-start_tar

        #check if the alignment produced a perfect match. if true perform the superimposition on the whole atom lists,
        #otherwise perform the superimposition on the atoms leading to continuous-matching subsequences of the alignment:
        #set k as the minimum length of the alignment's continuous subsequences you want to consider.
        #NB! max is used to fix the maximum value that k can assume, and refears to the maximum continuous subsequence length.

        if len(indexes_tar) == 0 and count_ref == len(ref_list_res):

            for res in targ_list_res:
                for atom in res:
                    targ_list.append(atom)

            ref_list = ref_list_atoms

        else:
            k = max -10

            #extract from the target structure the atoms leading to continuous match
            #subsequences of the alignment with length >= keys

            for index_tar, index_ref in zip(indexes_tar, indexes_ref):

                if index_tar[1]-index_tar[0] >= k :
                    for i in range(index_tar[0],index_tar[1]):

                        for atom in targ_list_res[i]:
                            targ_list.append(atom)

                    for i in range(index_ref[0],index_ref[1]):

                        for atom in ref_list_res[i]:
                            ref_list.append(atom)

        #resize the two lists to perform superimposition

        if len(targ_list)>=len(ref_list):
            targ_list = targ_list[:len(ref_list)]

        else:

            ref_list = ref_list[:len(targ_list)]

        #try superimpose
        rms,rot,tran = bp_center.superimpose_get_rotranrms(targ_list,ref_list)

        return rms, rot, tran, score

    def bp_center(receptor_pdb, chain, family, verbose=False):

        #perform an automatic detection of the orthosteric binding pocket's center for the passed chain of the passed structure

        #set the reference files directory
        #reference_files = bp_center.reference_files
        reference_files = '/home/rribeiro/Projects/pharmacodynamics/Tools/autogrid_reference_files'


        #NEW RECEPTORS FAMILIES SHOULD BE ADDED HERE!

        families = {
                        'IGR' : [reference_files+'/IGr_Ref.pdb', reference_files+'/IGr_Ref_grids.txt'],
                        'GPCR' : [reference_files+'/GPCr_Ref.pdb', reference_files+'/GPCr_Ref_grids.txt']
        }


        if family not in families:
            print("\nArgoument Error: -family argoument should be ", families.keys(),'\n\n')

        else:
            reference_structures_filename = (families[family])[0]
            reference_grids_filename = (families[family])[1]

        parser = PDBParser(PERMISSIVE=1)

        s1 = parser.get_structure("tar", receptor_pdb)
        s2 = parser.get_structure("ref", reference_structures_filename)

        tar_chain = s1[0][chain]
        tar_atoms_list = []
        tar_res_list = []

        for residue in tar_chain:
            tar_res_list.append(residue)
            for atom in residue:
                tar_atoms_list.append(atom)

        ref_atoms_dict = {}
        ref_res_dict = {}

        for ref_chain in s2[0]:
            residues = []
            atoms = []
            ref_res_dict.update({ref_chain.id : residues})
            ref_atoms_dict.update({ref_chain.id : atoms})
            for residue in ref_chain:
                ref_res_dict[ref_chain.id].append(residue)
                for atom in residue:
                    ref_atoms_dict[ref_chain.id].append(atom)

        grids_file_ref = open(reference_grids_filename, "r")

        #create a dictionary containing the coordinates of the reference grids (as vectors) which can be
        #modified during the execution in order to not to modify the original reference grids file,
        #and a dictionary that will contain the transformed coordinates of the grids.

        ref_grids_dic = {}
        for line in grids_file_ref :

            toks = line.split()
            chain_id = toks[0]
            coor = Vector(toks[2],toks[3],toks[4])
            ref_grids_dic.update({chain_id : coor})

        #Start the detection

        #set a reasonable high value to minimize the rmsd
        score_opt = 0
        rms_opt = 1000

        for ref_chain in s2[0]:

            #create copies in order to preserve the integrity of the original chains
            ref_atoms_list = ref_atoms_dict[ref_chain.id]
            ref_res_list = ref_res_dict[ref_chain.id]

            rms,rot,tran,score = bp_center.allign_3D(ref_atoms_list,ref_res_list,tar_res_list)

            if rms < rms_opt:
                rms_opt = rms
                rot_opt, tran_opt = rot, tran
                score_opt = score
                opt = ref_chain.id

        #set a threshold for the alignment score to judge the goodness of the calculated alignment
        if score_opt < 80 and rms_opt > 7:

            print('Error: no good structural alignment found for chain ', chain)
            print(score_opt,rms_opt)
        else:

            #read the reference coordinates file and transform the relative coordinates currently stored in the dictionary,
            #then write them to the output target's grids coordinates file.
            #NB: this is important , because all the transformations performed of the reference chains change the
            #position of the chains themselves, but not the relative grids coordinates!

            ref_grid_coor = ref_grids_dic[opt]
            targ_grid_coor = list(Vector.right_multiply(ref_grid_coor, rot_opt) + tran_opt)

            if verbose :

                #print summary
                print("###############################################################")
                print("                                                               ")
                print("Target chain '"+chain+"' Summary :                             ")
                print("                                                               ")
                print("reference chain used : '"+opt+"'"                               )
                print("calculated rmsd :",rms_opt,                                     )
                print("calculated score alignment :",score_opt,                        )
                print("grid center coordinates : ",targ_grid_coor,                     )
                print("                                                               ")
                print("###############################################################")
                print("                                                               ")

            #return the coordinates of the binding pocket's center for the passed chain
            return targ_grid_coor

class binding:

    def agonist(receptor_conc, lig_conc_range, pkd_agonist):
        data=[]
        for conc in lig_conc_range:
            data.append(LR_eq_conc(receptor_conc, conc, 0, pkd_agonist, 0))
        return data

    def agonist_test(receptor_conc, lig_conc_range, pkd_agonist):
        data=[]
        for conc in lig_conc_range:
            data.append(LR_eq_conc_test(receptor_conc, conc, 0, pkd_agonist, 0))
        return data

    def direct(LR,receptor_conc, lig_conc_min, lig_conc_max, dlscore_results):
        lig_conc_range = np.geomspace(lig_conc_min, lig_conc_max, 20)
        binding_data = {}
        for entry in LR:
            data = binding.agonist(receptor_conc, lig_conc_range, dlscore_results[entry]['dlscore'])
            
            raw_data = {'x':lig_conc_range, 'y':data, 'label':entry}
        
            #fit sigmoid to data
            def sigmoid_binding(X, Bottom, Top, Kd, p):
                return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

            popt, pcov = curve_fit(sigmoid_binding, lig_conc_range, data, bounds=([np.min(data),-np.inf,-np.inf, 0.5],[np.inf,np.max(data),np.inf, 2.5]))

            xfit = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
            yfit = sigmoid_binding(xfit, *popt)

            fitted_data ={'x':xfit, 'y':yfit, 'label':entry + ' (fitted)'}
            
            slope = popt[1]
            #if slope is negative the constant is ki. If it is positive the constant is kd
            if slope <0:
                binding_data[entry] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'ki (μM)': round(popt[2],2)}
            elif slope >0:
                binding_data[entry] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'kd (μM)': round(popt[2],2)}
            else:
                raise Exception('Error. Slope = 0')
        return binding_data

    def direct_test(LR,receptor_conc, lig_conc_min, lig_conc_max, dlscore_results):
        lig_conc_range = np.geomspace(lig_conc_min, lig_conc_max, 20)
        binding_data = {}
        for entry in LR:
            data = binding.agonist_test(receptor_conc, lig_conc_range, dlscore_results[entry]['dlscore'])
            
            raw_data = {'x':lig_conc_range, 'y':data, 'label':entry}
        
            #fit sigmoid to data
            def sigmoid_binding(X, Bottom, Top, Kd, p):
                return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

            popt, pcov = curve_fit(sigmoid_binding, lig_conc_range, data, bounds=([np.min(data),-np.inf,-np.inf, 0.5],[np.inf,np.max(data),np.inf, 2.5]))

            xfit = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
            yfit = sigmoid_binding(xfit, *popt)

            fitted_data ={'x':xfit, 'y':yfit, 'label':entry + ' (fitted)'}
            
            slope = popt[1]
            #if slope is negative the constant is ki. If it is positive the constant is kd
            if slope <0:
                binding_data[entry] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'ki (μM)': round(popt[2],2)}
            elif slope >0:
                binding_data[entry] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'kd (μM)': round(popt[2],2)}
            else:
                raise Exception('Error. Slope = 0')
        return binding_data

    def binding_assay(agonist, antagonists, receptor_conc, agonist_submaximal, ant_conc_min, ant_conc_max, dlscore_results):
        ant_conc_range = np.geomspace(ant_conc_min, ant_conc_max, 40)

        def agonist_antagonist(receptor_conc, lig_bound_conc,ant_conc_range, pkd_ligand_bound, pkd_ligand):
            binding_data = []
            for conc in ant_conc_range:
                binding_data.append(LR_eq_conc(receptor_conc, lig_bound_conc, conc, pkd_ligand_bound, pkd_ligand))
            return binding_data


        binding_data={}
        for antagonist in antagonists:
            data = agonist_antagonist(receptor_conc,agonist_submaximal,ant_conc_range,
                                               dlscore_results[agonist[0]]['dlscore'],
                                               dlscore_results[antagonist]['dlscore'])


            raw_data = {'x':ant_conc_range, 'y':data, 'label':agonist[0] +' + '+ antagonist}
            #fit sigmoid to data
            def sigmoid_binding(X, Bottom, Top, Kd, p):
                return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

            popt, pcov = curve_fit(sigmoid_binding, ant_conc_range, data, bounds=([np.min(data),-np.inf,-np.inf, 0.5],[np.inf,np.max(data),np.inf, 2.5]))

            xfit = np.geomspace(ant_conc_min, ant_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
            yfit = sigmoid_binding(xfit, *popt)

            fitted_data ={'x':xfit, 'y':yfit, 'label':agonist[0] +' + '+ antagonist + ' (fitted)'}

            #append results
            slope = popt[1]
            #if slope is negative the constant is ki. If it is positive the constant is kd
            if slope <0:
                binding_data[antagonist] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'ki (μM)': round(popt[2],5)}
            elif slope >0:
                binding_data[antagonist] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'kd (μM)': round(popt[2],5)}
            else:
                raise Exception('Error. Slope = 0')

        return binding_data
    
    def radioligand_binding_assay(radioligand_name, radioligand_submaximal_conc, ligands_names_list, ligand_conc_range, receptor_conc, ligands_pkd_list):
        from scipy.optimize   import curve_fit
        
        def LR_eq_calc(radioligand_submaximal_conc, pKd_radioligand, ligand_conc, pKd_ligand, receptor_conc):
            import math
            #from pKd to Kd
            Kd_ligand = (10**(-(float(pKd_ligand)))) * 10**6

            Kd_radioligand = ((10**(-(float(pKd_radioligand)))) * 10**6)*(1+(ligand_conc/Kd_ligand))

            #LR determination
            a = 1
            b = float(radioligand_submaximal_conc)+float(receptor_conc)+Kd_radioligand
            c = float(receptor_conc)*float(radioligand_submaximal_conc)
            delta = (b**2) - (4*a*c)
            LR = (b-math.sqrt(delta))/(2*a)
            return LR
        
        def radioligand_occupancy_calc(radioligand_submaximal_conc, pKd_radioligand, ligand_conc_range, pKd_ligand, receptor_conc):
            binding_data = []
            for ligand_conc in ligand_conc_range:
                LR_conc_init = LR_eq_calc(radioligand_submaximal_conc, pKd_radioligand, ligand_conc, pKd_ligand, receptor_conc)
                binding_data.append(LR_conc_init)
            
            return binding_data


        binding_data={}     
        for ligand in ligands_names_list:
            
            data = radioligand_occupancy_calc(radioligand_submaximal_conc, ligands_pkd_list[radioligand_name]['pKd'],
                                            ligand_conc_range, ligands_pkd_list[ligand]['pKd'],
                                            receptor_conc)


            raw_data = {'x':ligand_conc_range, 'y':data, 'label':ligand}
            #fit sigmoid to data
            def sigmoid_binding(X, Bottom, Top, Kd, p):
                return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

            popt, pcov = curve_fit(sigmoid_binding, ligand_conc_range, data, bounds=([np.min(data),-np.inf,-np.inf, 0.5],[np.inf,np.max(data),np.inf, 2.5]))

            xfit = np.geomspace(ligand_conc_range.min(), ligand_conc_range.max(), 50000) # These values are the same as the values for the simulation time and not ligand concentration
            yfit = sigmoid_binding(xfit, *popt)

            fitted_data ={'x':xfit, 'y':yfit, 'label':ligand + ' (fitted)'}

            #append results
            slope = popt[1]
            #if slope is negative the constant is ki. If it is positive the constant is kd
            if slope <0:
                binding_data[ligand] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'ki': round(popt[2],5)}
            elif slope >0:
                binding_data[ligand] = {'raw_data': raw_data, 'fitted_data':fitted_data, 'kd': round(popt[2],5)}
            else:
                raise Exception('Error. Slope = 0')

        return binding_data

    def constants(binding_data):
        kvalues={}
        for ligand in binding_data:
            constant = list(binding_data[ligand].keys())[-1]
            value = binding_data[ligand][constant]
            kvalues[ligand]={constant:value}
        return kvalues

    def maxbend(binding_data, l_conc_range):
        from scipy.optimize   import curve_fit, minimize


        def sigmoid(X, Bottom, Top, Kd, p):
            return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

        xfit = np.geomspace(1E-3, 1E4, 50000) #warning: shoud this be the minimum and maximum of concentration
        popt, pcov = curve_fit(sigmoid, l_conc_range, binding_data, bounds=([np.min(binding_data),-np.inf,-np.inf, 0.5],[np.inf,np.max(binding_data),np.inf, 2.5]))


        def sigmoid_deriv_b(x, a,d,c,b):
            return (x/c)**b*(a - d)*np.log(x/c)/((x/c)**b + 1)**2

        min_value = minimize(sigmoid_deriv_b, np.max(xfit), args=(popt[0],popt[1],popt[2],popt[3]), method = 'Nelder-Mead')

        submaximal = round(min_value.x[0],3)
        return submaximal

    def curve(binding_data):
        import plotly
        import plotly.graph_objs as go

        colors = plotly.colors.DEFAULT_PLOTLY_COLORS

        plot_data=[]

        color_id=0

        for ligand in binding_data:
            trace_raw = go.Scatter(x=binding_data[ligand]['raw_data']['x'],
                                y=minmax_scale(binding_data[ligand]['raw_data']['y']),
                                mode='markers',
                                showlegend=True,
                                name=binding_data[ligand]['raw_data']['label'],
                                marker=dict(color=colors[color_id]))

            plot_data.append(trace_raw)

            trace_fitted = go.Scatter(x=binding_data[ligand]['fitted_data']['x'],
                                y=minmax_scale(binding_data[ligand]['fitted_data']['y']),
                                mode='lines',
                                showlegend=False,
                                name=binding_data[ligand]['fitted_data']['label'],
                                line=dict(color=colors[color_id]))

            plot_data.append(trace_fitted)

            color_id +=1

        layout = dict(title = '',
                      xaxis = dict(
                          title = '[Ligand] μM',
                          range = [-3, 6],
                          type='log',
                          titlefont=dict(
                              size=20
                          ),
                          tickfont=dict(
                              size=20
                          )),
                      yaxis = dict(
                          title = '[LR]/[R0]',
                          range = [0, 1],
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

        fig = go.FigureWidget(data=plot_data, layout=layout)
        return fig

 

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
                printProgressBar(0, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                for idx in range(len(ligand_conc_range)):

                    ligand_conc = ligand_conc_range[idx]

                    if kinetics == False:
                        #get LR conc
                        parameters={'R_init':receptor_conc}
                        LR_conc_init = LR_eq_conc(receptor_conc, ligand_conc, 0, pkd_values[ligands.index(ligand)], 0)
                        
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
                    printProgressBar(idx + 1, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                #sim_data['ligand_name'] = ligand_name
                #sim_data['pkd']=str(pkd_data[ligand_name])
                #sim_data['sim_data']=data
                #all_data.append(sim_data)
                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':ligand_name,
                                            'pkd':str(pkd_values[ligands.index(ligand)])}


            return simulation_data
        
        def run_test(ligands, dlscore_results, network_name, receptor_conc, ligand_conc_range, ttotal, nsteps):
            '''
            This function runs the network simulation and gives the raw output from it.

            ligands         list od ligands or receptors (this represent the number of simulations performed)
            pkd_data           data from dlscore (path of csv file)
            network_name       name of the netwokr (eg. G_alpha_s or G_alpha_i)
            lig_conc_range         range of the concentration of ligand.
            ttotal             total time mof simulation
            nsteps             number of time steps/windows'''

            available_networks = ['Gs', 'Gi', 'Gq']
            if network_name == 'Gz(Gi)': network_name = 'Gi'
            if network_name not in available_networks: raise Exception('Unvailable Network. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')

            # parse pkd_data
            pkd_data={}
            for key in dlscore_results.keys():
                pkd_data[key]=dlscore_results[key]['dlscore']

            #Main function
            pathway = importlib.import_module('.'+network_name, package='src.metnetlib.pathways')
            t = pl.geomspace(0.00001, ttotal, num=nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points

            simulation_data={}
            
            for ligand in ligands:
                sim_data={}
                ligand_name = os.path.splitext(ligand)[0]
                data=[]
                printProgressBar(0, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                for i_conc in range(len(ligand_conc_range)):

                    ligand_conc = ligand_conc_range[i_conc]

                    #get LR conc
                    LR_conc_init = LR_eq_conc_test(receptor_conc, ligand_conc, 0, pkd_data[ligand_name], 0)
                    mymodel = pathway.network(receptor_conc, LR_conc_init)
                    simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                    yout = simres.all

                    d1={'ligand_conc':ligand_conc, 'time':t }

                    for i in range(len(pathway.list_of_observables)):
                        d2={pathway.list_of_observables[i]:yout[pathway.list_of_observables[i]]}
                        d1.update(d2)
                    data.append(d1)
                    printProgressBar(i_conc + 1, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                #sim_data['ligand_name'] = ligand_name
                #sim_data['pkd']=str(pkd_data[ligand_name])
                #sim_data['sim_data']=data
                #all_data.append(sim_data)
                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':ligand_name,
                                            'pkd':str(pkd_data[ligand_name])}


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

    class inhibition:
        def run(agonist, antagonists, dlscore_results, network_name, receptor_conc, ligand_conc_range, agonist_submaximal_conc, ttotal, nsteps):

            #check network
            available_networks = ['Gs', 'Gi', 'Gq']
            if network_name == 'Gz(Gi)': network_name = 'Gi'
            if network_name not in available_networks: raise Exception('Unvailable Network. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')

            # parse pkd_data
            pkd_data={}
            for key in dlscore_results.keys():
                pkd_data[key]=dlscore_results[key]['dlscore']

            #Main function
            pathway = importlib.import_module('.'+network_name, package='src.metnetlib.pathways')
            t = pl.geomspace(0.00001, ttotal, num=nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points

            simulation_data={}
            agonist_name=os.path.splitext(agonist[0])[0]
            for ligand in antagonists:
                sim_data={}
                ligand_name = os.path.splitext(ligand)[0]
                data=[]
                printProgressBar(0, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                for i_conc in range(len(ligand_conc_range)):

                    ligand_conc = ligand_conc_range[i_conc]

                    #get LR conc
                    LR_conc_init = LR_eq_conc(receptor_conc, agonist_submaximal_conc, ligand_conc, pkd_data[agonist_name], pkd_data[ligand_name])
                    mymodel = pathway.network(receptor_conc, LR_conc_init)
                    simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                    yout = simres.all

                    d1={'ligand_conc':ligand_conc, 'time':t }

                    for i in range(len(pathway.list_of_observables)):
                        d2={pathway.list_of_observables[i]:yout[pathway.list_of_observables[i]]}
                        d1.update(d2)
                    data.append(d1)
                    printProgressBar(i_conc + 1, len(ligand_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                #sim_data['ligand_name'] = ligand_name
                #sim_data['pkd']=str(pkd_data[ligand_name])
                #sim_data['sim_data']=data
                #all_data.append(sim_data)
                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':agonist_name+' + ' + ligand_name,
                                            'pkd':str(pkd_data[ligand_name])}

            return simulation_data

        def analysis(simulation_data, lig_conc_range, pathway):
            from sklearn.preprocessing import minmax_scale

            # Define all the lists and dictionaries used in this function
            raw_data=[]
            normalized_data=[]
            fitted_data=[]

            dose={}

            #defining concentration range
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

                popt_IC50, pcov = curve_fit(equation_dose, lig_conc_range, metabolite_conc_norm, bounds=([np.min(metabolite_conc_norm),-np.inf,-np.inf, 0.5],[np.inf,np.max(metabolite_conc_norm),np.inf, 2.5]))

                xfit_IC50 = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
                yfit_IC50 = equation_dose(xfit_IC50, *popt_IC50)

                fit_IC50={'x':xfit_IC50, 'y':yfit_IC50, 'label':simulation_data[ligand]['label']}




                dose[ligand] = {'raw_data': raw_data_dict,
                                'normalized_data':normalized_data_dict ,
                                'fitted_data': fit_IC50,
                                'IC50 (μM)': round(popt_IC50[2],2),
                                'pA2': round(-np.log10(popt_IC50[2]*1E-6),2)}

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

            fig = go.FigureWidget(data=plot_data, layout=layout)
            return fig

        def constants(simulation_data):
            kvalues={}
            for ligand in simulation_data:
                IC50 = list(simulation_data[ligand].keys())[-2]
                IC50_value = simulation_data[ligand][IC50]
                pA2 = list(simulation_data[ligand].keys())[-1]
                pA2_value = simulation_data[ligand][pA2]
                kvalues[ligand]={IC50:IC50_value, pA2:pA2_value}

            return kvalues

    class agonist_antagonist:
        def run(agonists_list, antagonists_list, ligands_pKd_list, network_name, receptor_conc, ligand_conc_range, ligand_submaximal_conc, ttotal, nsteps, radioligand_action):
            
            #function
            def LR_eq_calc(radioligand_submaximal_conc, pKd_radioligand, ligand_conc, pKd_ligand, receptor_conc):
                import math
                #from pKd to Kd
                Kd_ligand = (10**(-(float(pKd_ligand)))) * 10**6

                Kd_radioligand = ((10**(-(float(pKd_radioligand)))) * 10**6)*(1+(ligand_conc/Kd_ligand))

                #LR determination
                a = 1
                b = float(radioligand_submaximal_conc)+float(receptor_conc)+Kd_radioligand
                c = float(receptor_conc)*float(radioligand_submaximal_conc)
                delta = (b**2) - (4*a*c)
                LR = (b-math.sqrt(delta))/(2*a)
                return LR

            #check network
            available_networks = ['Gs', 'Gi', 'Gq']
            if network_name == 'Gz(Gi)': network_name = 'Gi'
            if network_name not in available_networks: raise Exception('Unvailable Network. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')

            # parse pkd_data
            pkd_data={}
            for key in ligands_pKd_list.keys():
                pkd_data[key]=ligands_pKd_list[key]['pKd']
            print(pkd_data)

            #Main function
            pathway = importlib.import_module('.'+network_name, package='src.metnetlib.pathways')
            t = pl.geomspace(0.00001, ttotal, num=nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points

            simulation_data={}
            for agonist in agonists_list:
                #agonist_name=os.path.splitext(agonist[0])[0]
                for antagonist in antagonists_list:
                    sim_data={}
                    #antagonist_name = os.path.splitext(antagonist)[0]
                    interaction_name = agonist+'_'+antagonist
                    data=[]
                    printProgressBar(0, len(ligand_conc_range), prefix = "{:<15}".format(interaction_name[:15]), suffix = 'Complete', length = 50)

                    for i_conc in range(len(ligand_conc_range)):

                        ligand_conc = ligand_conc_range[i_conc]

                        #get LR conc
                        if radioligand_action == 'agonist':
                            LR_conc_init = LR_eq_calc(ligand_submaximal_conc, pkd_data[agonist], ligand_conc, pkd_data[antagonist], receptor_conc)
                        elif radioligand_action == 'antagonist':
                            LR_conc_init = LR_eq_calc(ligand_conc, pkd_data[agonist], ligand_submaximal_conc, pkd_data[antagonist], receptor_conc)
                        else: 
                            raise TypeError("set radioligand action: agonist or antagonist") 

                        mymodel = pathway.network(receptor_conc, LR_conc_init)
                        simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                        yout = simres.all

                        d1={'ligand_conc':ligand_conc, 'time':t }

                        for i in range(len(pathway.list_of_observables)):
                            d2={pathway.list_of_observables[i]:yout[pathway.list_of_observables[i]]}
                            d1.update(d2)
                        data.append(d1)
                        printProgressBar(i_conc + 1, len(ligand_conc_range), prefix = "{:<15}".format(interaction_name[:15]), suffix = 'Complete', length = 50)

                simulation_data[interaction_name] = {'sim_data':data,
                                            'label':interaction_name,
                                            'pkd_agonist':str(pkd_data[agonist]),
                                            'pkd_antagonist':str(pkd_data[antagonist])}

            return simulation_data
            
        def analysis(simulation_data, lig_conc_range, pathway):
            from sklearn.preprocessing import minmax_scale

            # Define all the lists and dictionaries used in this function
            raw_data=[]
            normalized_data=[]
            fitted_data=[]

            dose={}

            #defining concentration range
            #defining concentration range
            lig_conc_min = lig_conc_range.min()
            lig_conc_max = lig_conc_range.max()

            #Main function
            for interaction in simulation_data:

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
                raw_data_dict['label']=simulation_data[interaction]['label']

                normalized_data_dict['x']=lig_conc_range
                normalized_data_dict['y']=metabolite_conc_norm
                normalized_data_dict['label']=simulation_data[interaction]['label']

                ## create a list of all data
                raw_data.append(raw_data_dict)
                normalized_data.append(normalized_data_dict)

                ##Fitting curve to the data

                def equation_dose(X, Bottom, Top, EC50, p):
                    return Bottom + (Top-Bottom)/(1+np.power((EC50/X),p))

                popt_IC50, pcov = curve_fit(equation_dose, lig_conc_range, metabolite_conc_norm, bounds=([np.min(metabolite_conc_norm),-np.inf,-np.inf, 0.5],[np.inf,np.max(metabolite_conc_norm),np.inf, 2.5]))

                xfit_IC50 = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
                yfit_IC50 = equation_dose(xfit_IC50, *popt_IC50)

                fit_IC50={'x':xfit_IC50, 'y':yfit_IC50, 'label':simulation_data[ligand]['label']}




                dose[ligand] = {'raw_data': raw_data_dict,
                                'normalized_data':normalized_data_dict ,
                                'fitted_data': fit_IC50,
                                'IC50': round(popt_IC50[2],4),
                                'pIC50': round(-np.log10(popt_IC50[2]*1E-6),4)}

            return dose

        def curve(simulation_data):
            import plotly
            import plotly.graph_objs as go

            colors = plotly.colors.DEFAULT_PLOTLY_COLORS

            plot_data=[]

            color_id=0
            for interaction in simulation_data:
                trace_norm = go.Scatter(x=simulation_data[interaction]['normalized_data']['x'],
                                    y=minmax_scale(simulation_data[interaction]['normalized_data']['y'])*100 ,
                                    mode='markers',
                                    showlegend=True,
                                    name=simulation_data[interaction]['normalized_data']['label'],
                                    marker=dict(color=colors[color_id]))
                plot_data.append(trace_norm)

                trace_fitted = go.Scatter(x=simulation_data[interaction]['fitted_data']['x'],
                                    y=minmax_scale(simulation_data[interaction]['fitted_data']['y'])*100,
                                    mode='lines',
                                    showlegend=False,
                                    name=simulation_data[interaction]['fitted_data']['label'],
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

            fig = go.FigureWidget(data=plot_data, layout=layout)
            return fig

        def constants(simulation_data):
            kvalues={}
            for interaction in simulation_data:
                IC50 = list(simulation_data[interaction].keys())[-2]
                IC50_value = simulation_data[interaction][IC50]
                pIC50 = list(simulation_data[interaction].keys())[-1]
                pIC50_value = simulation_data[interaction][pIC50]
                kvalues[interaction]={IC50:IC50_value, pIC50:pIC50_value}

            return kvalues
