#!/usr/bin/env python
"""Generate 3D structures of mutant based on the 3D of wildtype
    
"""
# Original source code: https://codeocean.com/capsule/2202829/tree/v1
# Modified by:  Duc Nguyen (@NguyenLab)

import os
import sys 
import shutil
import pandas as pd
import myconfig
from pqdm.processes import pqdm


mypath = myconfig.mypath()

os.environ["JACKALDIR"] = mypath['JACKALDIR']
profix = mypath['profix']
scap = mypath['scap']
vmd = mypath['vmd']


class gen_3D:
    """
    Generate 3D structures for mutants and
    create the mutant sites and wildtype sites around
    the mutant location
    
    Parameters
    ----
    name: str
        Name of the dataset. E.g.: S1131, AB_Bind, etc.
    mutants: list(dict)
        Example: [{'pdbid' : '1AK4',
                   'chainid': 'D',
                   'resid': '488',
                   'wildname': 'A',
                   'mutname': 'G',
                   'ddG': 1.010373}]          
                   
    wildtype_path: str
        Location of all widltype PDB files
        
    output_path: str
        Location for storing the generated structures
        
    cutoff: float
        Define the region around the mutant position
    """
    def __init__(self, name, mutants, wildtype_path, output_path, cutoff=12):
        self.name = name
        self.mutants = mutants
        self.wildtype_path = wildtype_path
        self.output_path = output_path
        self.cutoff = cutoff

        
    @staticmethod
    def modifypdb(name1, name2):
        """Simplify the original pdb file <name1>
            by extracting ATOM lines only.
        """
        pdbfile = open(name1)
        outfile = open(name2, "w")
        lines = pdbfile.read().splitlines()
        PRO = []
        for line in lines:
            if line[0:4] == 'ATOM':
                outfile.write(line+'\n')

        outfile.close()                    
        return 


    def prepare_dir(self, mut_info):
        """
        Create folder based on the mutant inforation, and
        copy the original wildtype pdb to that folder
            
        Parameters
        ----------
        mut_info: dict
            Example: {'pdbid' : '1AK4',
                   'chainid': 'D',
                   'resid': '488',
                   'wildname': 'A',
                   'mutname': 'G',
                   'ID': '1AK4_D_A488G',
                   }
        """
        dir_name = f"{self.output_path}/{mut_info['ID']}"
        if not os.path.isdir(f'{dir_name}'):
            os.mkdir(f'{dir_name}')
        self.modifypdb(f"{self.wildtype_path}/{mut_info['pdbid']}.pdb", 
                  f"{dir_name}/{mut_info['pdbid']}.pdb")
        return dir_name
    
    
    @staticmethod
    def gen_mut(workdir, mut_info, dir_name):
        """Generate 3D structure of a mutant
        """
        os.chdir(workdir)
        chainid, resid = mut_info['chainid'], mut_info['resid']
        mutname = mut_info['mutname']
        pdbid = mut_info['pdbid']
        ID = mut_info['ID']
        pdb_file = pdbid+'.pdb'
        os.system(f'{profix} -fix 0 {dir_name}/{pdb_file}')
        os.system(f'mv {pdbid}_fix.pdb {ID}_wild.pdb')

        scapfile = open(f'{workdir}/tmp_scap.list', 'w')
        scapfile.write(f'{chainid},{resid},{mutname}')
        scapfile.close()   
        os.system(f'{scap} -ini 20 -min 4  {pdb_file} ./tmp_scap.list')
        os.system(f'mv {pdbid}_scap.pdb {ID}_mut.pdb')
        os.system('rm ./tmp_scap.list')
        
        
    @staticmethod
    def def_site(workdir, mut_info, wildtype, cutoff):
        """
        Generate the local sites within <cutoff> around mutant position
        """
        os.chdir(workdir)
        chainid = mut_info['chainid']
        resid = mut_info['resid']
        if wildtype:
            pdbid = mut_info['ID'] + '_wild'
            #pdbid = f'{pdbid}_wild'
        else:
            pdbid = mut_info['ID'] + '_mut'
            #pdbid = f'{pdbid}_mut'
        tclname = 'def_site.tcl'
        filename = f'{pdbid}.pdb'
        b_name = pdbid +'_bindingsite.pdb'
        m_name = pdbid +'_mutationsite.pdb'
        tclfile = open(tclname,'w')
        print(chainid)
        tclfile.write("mol new {" + filename +"} type {pdb} first 0 last 0 step 1 waitfor 1\n")
        tclfile.write(f'set prot [atomselect top "within {cutoff} of chain {chainid}"]\n')
        tclfile.write('$prot writepdb '+ b_name +'\n')
        tclfile.write(f'set prot2 [atomselect top "within {cutoff} of resid {resid} and chain {chainid}"]\n')
        tclfile.write('$prot2 writepdb '+ m_name +'\n')
        tclfile.write('exit')
        tclfile.close()
        os.system(f"{vmd} -dispdev text -e {tclname}")

        
    def gen_3D_single(self, mut_info):
        print(mut_info)
        if mut_info['pdbid']:
            dir_name = self.prepare_dir(mut_info)
            workdir = f"{self.output_path}/{mut_info['ID']}"
            self.gen_mut(dir_name, mut_info, dir_name)
            print(mut_info)
            self.def_site(workdir , mut_info, True, self.cutoff) # wild type
            self.def_site(workdir , mut_info, False, self.cutoff) # mutant type
        

    def gen_3D_all(self, n_jobs=1):
        results = pqdm(self.mutants, self.gen_3D_single, n_jobs=n_jobs)        
        print('All Done!')
        
        
    def get_labels(self):
        PDBIDs = []
        ddGs = []
        for mut_info in self.mutants:
            if mut_info['ID']:
                PDBIDs.append(mut_info['ID'])
                ddGs.append(mut_info['ddG'])
        df = pd.DataFrame({'ID': PDBIDs, 'ddG': ddGs})
        df.to_csv(f'{self.output_path}/{self.name}.csv', float_format='%.3f', index=None)