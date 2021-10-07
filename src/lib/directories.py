'''Here we have to define the programs paths'''
import os

cwd = os.getcwd()

mainpath = '/usr/local/opt/'

vina_executable = mainpath+'vina1.1.2/bin/vina'
vina_split_executable = mainpath+'vina1.1.2/bin/vina_split'
prepare_receptor_executable = mainpath+'MGLTools1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
prepare_receptor_flex_executable = mainpath+'MGLTools1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_flexreceptor4.py'
prepare_ligand_executable = mainpath+'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
dlscore_path = mainpath+'DLSCORE/'
autogrids_reference_dir = cwd+'/src/metnetlib/autogrids/reference_files/'

## Set here the databases paths
#Human Target Receptors DataBase directory PATH
HuTRdb_path = os.path.join(os.getcwd(),'SSBColab/src/databases/HuTRdb.sqlite3')
swissprot_db = os.path.join(os.getcwd(),'src/databases/blastdb/swissprot.00')