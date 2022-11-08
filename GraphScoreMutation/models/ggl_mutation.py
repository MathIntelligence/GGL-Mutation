#!/usr/bin/env python

"""
Introduction:
    GGL-Mutation: Geometric Graph Learning for Protein Mutation

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    Oct 4, 2022

"""

import numpy as np
import pandas as pd
# 4import os
# from os import listdir
# from rdkit import Chem
from scipy.spatial.distance import cdist
from itertools import product
# from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from biopandas.pdb import PandasPdb


# from biopandas.mol2 import PandasMol2


class KernelFunction:

    def __init__(self, kernel_type='exponential_kernel',
                 kappa=2.0, tau=1.0):
        self.kernel_type = kernel_type
        self.kappa = kappa
        self.tau = tau

        self.kernel_function = self.build_kernel_function(kernel_type)

    def build_kernel_function(self, kernel_type):
        if kernel_type[0] in ['E', 'e']:
            return self.exponential_kernel
        elif kernel_type[0] in ['L', 'l']:
            return self.lorentz_kernel

    def exponential_kernel(self, d, vdw_radii):
        eta = self.tau * vdw_radii

        return np.exp(-(d / eta) ** self.kappa)

    def lorentz_kernel(self, d, vdw_radii):
        eta = self.tau * vdw_radii

        return 1 / (1 + (d / eta) ** self.kappa)


class GGLMutation:
    protein_atom_types_df = pd.read_csv(
        '../utils/protein_atom_types.csv')

    protein_atom_types = protein_atom_types_df['AtomType'].tolist()
    protein_atom_radii = protein_atom_types_df['Radius'].tolist()

    protein_atom_type_pair = [
        i[0] + "-" + i[1] for i in product(protein_atom_types, protein_atom_types)]

    def __init__(self, Kernel, chain, residue_id):
        self.Kernel = Kernel

        self.chain = chain
        self.residue_id = residue_id

        self.pairwise_atom_type_radii = self.get_pairwise_atom_type_radii()

    def get_pairwise_atom_type_radii(self):

        protein_atom_radii_dict = {a: r for (a, r) in zip(
            self.protein_atom_types, self.protein_atom_radii)}

        pairwise_atom_type_radii = {i[0] + "-" + i[1]: protein_atom_radii_dict[i[0]] +
                                                       protein_atom_radii_dict[i[1]] for i in
                                    product(self.protein_atom_types, self.protein_atom_types)}

        return pairwise_atom_type_radii

    def pdb_to_df(self, pdb_file: str) -> pd.DataFrame:
        ppdb = PandasPdb()
        ppdb.read_pdb(pdb_file)
        ppdb_all_df = ppdb.df['ATOM']
        ppdb_df = ppdb_all_df[ppdb_all_df['atom_name'].isin(
            self.protein_atom_types)]
        atom_index = ppdb_df['atom_number']
        atom_name = ppdb_df['atom_name']
        chain_id = ppdb_df['chain_id']
        residue_name = ppdb_df['residue_name']
        residue_id = ppdb_df['residue_number']
        x, y, z = ppdb_df['x_coord'], ppdb_df['y_coord'], ppdb_df['z_coord']
        df = pd.DataFrame.from_dict({'ATOM_INDEX': atom_index, 'ATOM_NAME': atom_name,
                                     'RESIDUE': residue_name, 'CHAIN': chain_id,
                                     'RESIDUE_ID': residue_id,
                                     'X': x, 'Y': y, 'Z': z})

        return df

    def binding_site_atoms(self, protein_df):
        protein_1 = protein_df[protein_df['CHAIN'] == self.chain]

        protein_2 = protein_df[protein_df['CHAIN'] != self.chain]

        return protein_1, protein_2

    def mutation_site_atoms(self, protein_df):
        protein_1 = protein_df[(protein_df['CHAIN'] == self.chain) & (
                protein_df['RESIDUE_ID'] == self.residue_id)]

        protein_2 = protein_df.loc[~protein_df.index.isin(protein_1.index)]

        return protein_1, protein_2

    def get_mwcg_rep(self, protein_1: pd.DataFrame, protein_2: pd.DataFrame) -> pd.DataFrame:
        atom_pairs = list(
            product(protein_1['ATOM_NAME'], protein_2['ATOM_NAME']))
        atom_pairs = [x[0] + "-" + x[1] for x in atom_pairs]
        pairwise_radii = [self.pairwise_atom_type_radii[x]
                          for x in atom_pairs]
        pairwise_radii = np.asarray(pairwise_radii)

        pairwise_mwcg = pd.DataFrame(atom_pairs, columns=["ATOM_PAIR"])
        distances = cdist(protein_1[["X", "Y", "Z"]],
                          protein_2[["X", "Y", "Z"]], metric="euclidean")
        pairwise_radii = pairwise_radii.reshape(
            distances.shape[0], distances.shape[1])
        mwcg_distances = self.Kernel.kernel_function(distances, pairwise_radii)

        distances = distances.ravel()
        mwcg_distances = mwcg_distances.ravel()
        mwcg_distances = pd.DataFrame(
            data={"DISTANCE": distances, "MWCG_DISTANCE": mwcg_distances})
        pairwise_mwcg = pd.concat([pairwise_mwcg, mwcg_distances], axis=1)
        # pairwise_mwcg = pairwise_mwcg[pairwise_mwcg["DISTANCE"] <= self.cutoff].reset_index(
        #     drop=True)

        return pairwise_mwcg

    def get_site_specific_mwcg(self, pdb_file, site):
        prot_df = self.pdb_to_df(pdb_file)

        if site[0] in ['b', 'B']:
            prot1, prot2 = self.binding_site_atoms(prot_df)
        elif site[0] in ['m', 'M']:
            prot1, prot2 = self.mutation_site_atoms(prot_df)
        else:
            # Need a better control here
            return False

        ssMWCG = self.get_mwcg_rep(prot1, prot2)

        return ssMWCG

    def get_site_specific_ggl_score(self, pdb_file, site):
        features = ['COUNTS', 'SUM', 'MEAN', 'STD', 'MIN', 'MAX']

        ss_mwcg = self.get_site_specific_mwcg(pdb_file, site)

        mwcg_temp_grouped = ss_mwcg.groupby('ATOM_PAIR')
        mwcg_temp_grouped.agg(['sum', 'mean', 'std', 'min', 'max'])
        mwcg_temp = mwcg_temp_grouped.size().to_frame(name='COUNTS')

        mwcg_temp = (mwcg_temp
                     .join(mwcg_temp_grouped.agg({'MWCG_DISTANCE': 'sum'}).rename(columns={'MWCG_DISTANCE': 'SUM'}))
                     .join(mwcg_temp_grouped.agg({'MWCG_DISTANCE': 'mean'}).rename(columns={'MWCG_DISTANCE': 'MEAN'}))
                     .join(mwcg_temp_grouped.agg({'MWCG_DISTANCE': 'std'}).rename(columns={'MWCG_DISTANCE': 'STD'}))
                     .join(mwcg_temp_grouped.agg({'MWCG_DISTANCE': 'min'}).rename(columns={'MWCG_DISTANCE': 'MIN'}))
                     .join(mwcg_temp_grouped.agg({'MWCG_DISTANCE': 'max'}).rename(columns={'MWCG_DISTANCE': 'MAX'}))
                     )

        mwcg_columns = {'ATOM_PAIR': self.protein_atom_type_pair}

        for _f in features:
            mwcg_columns[_f] = np.zeros(len(self.protein_atom_type_pair))
        ggl_score = pd.DataFrame(data=mwcg_columns)
        ggl_score = ggl_score.set_index('ATOM_PAIR').add(
            mwcg_temp, fill_value=0).reindex(self.protein_atom_type_pair).reset_index()

        return ggl_score

    # def get_ggl_mutation_score(self, prot_wild, prot_mutant):

    #     #get site-specific MWCG for both wild and mutant protein
    #     ssmwcg = self.get_site_specific_mwcg(prot_df)
