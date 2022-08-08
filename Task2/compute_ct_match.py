#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np

def id_celltype_match(solution_file, datafile):
    # Returns a dict of cell type labels to list of ids in the object
    from collections import defaultdict

    ct_id_dict = defaultdict(list)

    ad_sol = sc.read(solution_file)
    adata = sc.read(datafile)

    obs_scrambler = ad_sol.uns['pairing_ix'].astype(int)

    #subset to same observations
    adata_sub = adata[ad_sol.obs_names]
    
    # Get mapping of cell type to IDs in the matching object
    for idx,ct in enumerate(adata_sub.obs['cell_type'].iloc[obs_scrambler]):
        ct_id_dict[ct].append(idx)

    return ct_id_dict, obs_scrambler


def prep_res_df(metafile):
    metadata = pd.read_csv(metafile)
    
    res = metadata[['ATAC2GEX', 'Team Name', 'id']].set_index('id')
    res['CT_match'] = np.zeros(res.shape[0]) - 1

    return res


def CT_score(filename, ct_id_dict, obs_scrambler):
    # Read in file output
    print(f'Reading file: {filename}')
    adata = sc.read(filename)
    
    ct_sum = 0

    # Sum weights over all cell type matches
    for ct in ct_id_dict.keys():
        # print(f'Cell type: {ct}')
        # print(f'type of obs_scrambler:{type(obs_scrambler)}')
        # print(f'with units: {type(obs_scrambler[0])}')
        # print(f'in the dict are: {type(ct_id_dict[ct])}')
        # print(f'dict vals are lists of: {type(ct_id_dict[ct][0])}')
        data = adata.X[obs_scrambler]
        ct_sum += data[ct_id_dict[ct],:][:,ct_id_dict[ct]].sum()
    
    val = ct_sum/adata.X.sum()
    
    return val


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Compute celltype similarity')

    parser.add_argument('-i', '--input_files', required=True, nargs='+')
    parser.add_argument('-d', '--data_file', required=True)
    parser.add_argument('-m', '--metadata_file', required=True)
    parser.add_argument('-s', '--solution_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)

    args = parser.parse_args()
    filenames = args.input_files
    datafile = args.data_file
    metafile = args.metadata_file
    outfile = args.output_file
    solfile = args.solution_file

    # Find mapping of ids to cell types
    ct_id_matches, obs_unscrambler = id_celltype_match(solfile, datafile)

    # Prep results object
    res = prep_res_df(metafile)

    res_ct = dict()
    # Loop through input files
    for fn in filenames:
        sub_id = int(fn.split('submission_')[1].split('/')[0])

        # Compute cell_type_matching score
        res_ct[sub_id] = CT_score(fn, ct_id_matches, obs_unscrambler)

    res['CT_match'] = pd.Series(res_ct)

    # Save output
    res.to_csv(outfile)
