#!/usr/bin/env python
"""Generate 3D structures of mutant based on the 3D of wildtype
    
"""
# Authors:  Duc Nguyen (@NguyenLab)

import _gen_3D
import argparse
import pandas as pd
from pqdm.processes import pqdm


class do_S1131():
    name = 'S1131'
    csv_path = '/home/ddng223/RI-Score-Mutation/data/S1131/labels/S1131.csv'
    wildtype_path = '/home/ddng223/RI-Score-Mutation/data/S1131/wild_types/XBindProfPaperData/pdbs'
    output_path = '/home/ddng223/RI-Score-Mutation/data/S1131'
    cutoff = 12.0

    def read_csv(self):
        df = pd.read_csv(self.csv_path, dtype={'protein':str})
        pdbids = df['protein']
        mutations = df['mutation']
        ddGs = df['DDG']
        self.mutants = []
        
        for _index, _pdbid in enumerate(pdbids):
            mutant = mutations[_index]
            chainid = mutant[0]
            wildname = mutant[2]
            mutname = mutant[-1]
            resid = mutant[3:-1]
            ddG = ddGs[_index]
            _dict = {'pdbid': _pdbid, 'chainid': chainid,
                    'resid': resid, 'wildname': wildname,
                    'mutname': mutname,
                    'ddG': ddG,
                    'ID': f'{_pdbid}_{chainid}_{wildname}{resid}{mutname}'}
            self.mutants.append(_dict)
        
        return self.mutants
    
    
    def gen_3D(self, n_jobs):
        self.read_csv()        
        data = _gen_3D.gen_3D(self.name, self.mutants, 
                              self.wildtype_path, self.output_path, self.cutoff)
        data.gen_3D_all(n_jobs)
        data.get_labels()

    
def main(data_name='S1131', n_jobs=1):
    if data_name == 'S1131':
        data = do_S1131()
        data.gen_3D(n_jobs)
        
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate 3D structures and collect labels for different datasets",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--data-name', help='name of data',
                        type=str)
    parser.add_argument('-n', '--n-jobs', help='number of cores',
                        type=int, default=1)
       
    args = parser.parse_args()
    
    main(args.data_name, args.n_jobs)
    
