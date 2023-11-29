#!/usr/bin/env python
"""Generate 3D structures of mutant based on the 3D of wildtype
    
"""
# Authors:  Duc Nguyen (@NguyenLab)  and Masud Rana (@NguyenLab)
# Date last modified: Nov 28, 2023

import _gen_3D
import argparse
import pandas as pd
from pqdm.processes import pqdm


class Get3dStructures():

    cutoff = 12.0

    def __init__(self, name, csv_path, wildtype_path, output_path):
        self.name = name
        self.csv_path = csv_path
        self.wildtype_path = wildtype_path
        self.output_path = output_path

    def read_csv(self):
        df = pd.read_csv(self.csv_path)
        pdbids = df['pdbid'].values
        chain_id = df['chain'].values
        res_id = df['residue_id'].values
        wildtype = df['wild_type'].values
        muttype = df['mutant'].values
        ddGs = df['ddg'].values
        self.mutants = []
        
        for _index, _pdbid in enumerate(pdbids):
            #mutant = mutations[_index]
            chainid = chain_id[_index]
            wildname = wildtype[_index]
            mutname = muttype[_index]
            resid = res_id[_index]
            ddG = ddGs[_index]
            _dict = {'pdbid': _pdbid, 'chainid': chainid,
                    'resid': resid, 'wildname': wildname,
                    'mutname': mutname,
                    'ddG': ddG,
                    'ID': f'{_pdbid}_{chainid}_{wildname}{resid}{mutname}'
                    }
            self.mutants.append(_dict)
        
        return self.mutants
    
    
    def gen_3D(self, n_jobs):
        self.read_csv()        
        data = _gen_3D.gen_3D(self.name, self.mutants, 
                              self.wildtype_path, self.output_path, self.cutoff)
        data.gen_3D_all(n_jobs)
        data.get_labels()
    
               

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate 3D structures and collect labels for different datasets",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--data-name', help='name of data',
                        type=str)
    parser.add_argument('-f', '--csv_path', help='path for data file in csv',
                        type=str)
    parser.add_argument('-w', '--wildtype_path', help='path for wild type PDB structures',
                        type=str)
    parser.add_argument('-o', '--output_path', help='path for output directory',
                        type=str)
    parser.add_argument('-n', '--n-jobs', help='number of cores',
                        type=int, default=1)
       
    args = parser.parse_args()
    
    #main(args.data_name, args.n_jobs)

    Get3D = Get3dStructures(name = args.data_name, csv_path = args.csv_path,
                            wildtype_path = args.wildtype_path, output_path = args.output_path)
    Get3D.gen_3D(args.n_jobs)
    
