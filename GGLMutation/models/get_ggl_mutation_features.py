#!/usr/bin/env python
"""
Introduction:
    Geometric Graph Learning Features for Protein Mutation Data

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    Oct 4, 2022

"""

import sys
# import pandas as pd
# from itertools import product
import ntpath
import argparse
import time

from ggl_mutation import *


class GGLMutationFeatures:
    abs_path = os.path.dirname(__file__)
    rel_path = '../kernels/kernels.csv'
    full_path = os.path.join(abs_path, rel_path)
    df_kernels = pd.read_csv(full_path)

    def __init__(self, args):

        self.kernel_index = args.kernel_index

        self.dataset = args.dataset

        self.data_folder = args.data_folder
        self.feature_folder = args.feature_folder

    @staticmethod
    def get_mutation_info(mutation):
        chain = mutation[1]
        wildtype = mutation[2][0]
        mutanttype = mutation[2][-1]
        resid = mutation[2][1:-1]
        resid = int(resid)

        return chain, wildtype, mutanttype, resid

    # def ggl_mutation_score(self):


    def get_ggl_mutation_features(self, parameters):

        Kernel = KernelFunction(kernel_type=parameters['type'],
                                kappa=parameters['power'], tau=parameters['tau'])

        data_df = pd.read_csv(self.dataset)

        pdbids = data_df['ID'].values
        # mutant_list = data_df['Mutation'].values
        ddG = data_df['ddG'].tolist()

        for index in range(data_df.shape[0]):
            pdbid = pdbids[index]
            mutation = pdbid.split('_')
            chain, wildtype, mutanttype, residue_id = self.get_mutation_info(mutation)

            gglMut = GGLMutation(Kernel, chain, residue_id)

            # file_folder = f'{pdbid}_{chain}_{wildtype}{residue_id}{mutanttype}'

            sites = ['bindingsite', 'mutantsite']
            prot_types = ['wild', 'mut']

            count = 0

            for site in sites:
                for prot in prot_types:
                    if prot == 'wild':
                        pdb_file = f'{self.data_folder}/{pdbid}/{pdbid}_{site}.pdb'
                    else:
                        pdb_file = f'{self.data_folder}/{pdbid}/{pdbid}_{prot}_{site}.pdb'

                    ggl_score = gglMut.get_site_specific_ggl_score(
                        pdb_file, site)

                    atom_pairs = ggl_score['ATOM_PAIR'].tolist()
                    features = ggl_score.columns[1:].tolist()

                    pairwise_features = [i[0] + '_' + i[1] + '_' + prot[0] + site[0]
                                         for i in product(atom_pairs, features)]

                    feature_values = ggl_score.drop(
                        ['ATOM_PAIR'], axis=1).values.flatten()

                    ggl_score_df = pd.DataFrame(columns=[pairwise_features])
                    ggl_score_df.loc[0] = feature_values

                    if count == 0:
                        ggl_score_df_all = ggl_score_df
                    else:
                        ggl_score_df_all = pd.concat([ggl_score_df_all, ggl_score_df], axis=1)

                    count += 1

            if index == 0:
                ggl_features_df = ggl_score_df_all
            else:
                ggl_features_df.loc[index] = ggl_score_df_all.values.flatten()

        ggl_features_df.insert(0, 'ID', pdbids)
        ggl_features_df.insert(1, 'ddG', ddG)

        return ggl_features_df

    def main(self):

        kernel_params = {
            'type': self.df_kernels.loc[self.kernel_index, 'type'],
            'power': self.df_kernels.loc[self.kernel_index, 'power'],
            'tau': self.df_kernels.loc[self.kernel_index, 'tau']
        }

        df_features = self.get_ggl_mutation_features(kernel_params)

        dataset_file_name_only = ntpath.basename(self.dataset).split('.')[0]

        output_file_name = f'{dataset_file_name_only}_ker{self.kernel_index}.csv'

        df_features.to_csv(f'{self.feature_folder}/{output_file_name}', index=False, float_format='%.5f')


def get_args(args):
    parser = argparse.ArgumentParser(description="Get GGL Mutation Features")

    parser.add_argument('-k', '--kernel-index', help='Kernel Index (see kernels/kernels.csv)',
                        type=int)

    parser.add_argument('-f', '--dataset',
                        help='path to CSV dataset')

    parser.add_argument('-dd', '--data_folder', type=str,
                        help='path to data folder directory')

    parser.add_argument('-fd', '--feature_folder', type=str,
                        help='path to the directory where features will be saved')

    args = parser.parse_args()

    return args


def cli_main():
    args = get_args(sys.argv[1:])

    ggl_mut_features = GGLMutationFeatures(args)

    ggl_mut_features.main()


if __name__ == "__main__":
    t0 = time.time()

    cli_main()

    print('Done!')
    print('Elapsed time: ', time.time() - t0)
