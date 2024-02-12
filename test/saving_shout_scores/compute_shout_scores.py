import run_shout
import argparse
import pickle
import pandas as pd

def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_path', type=str, help='Please specify path to anndata pickle data.')
    parser.add_argument('cluster_key', type=str, help='Please specify cluster key.')
    parser.add_argument('radii', type=int, nargs='+', help='Please specify list of radii within which to calculate local heterogeneity scores.')
    parser.add_argument('--normalize', type=bool, default=True, help='Optionally specify "normalize", default value is True.')
    parser.add_argument('--num_cell_types', type=int, default=None, help='Optionally specify no. of celltypes in sample.')
    parser.add_argument('--coord_type', type=str, choices=['generic', 'grid'], default='GENE_SYMBOL')
    parser.add_argument('--copy', type=bool, default=False, help='Optionally specify whether heterogenety scores to be written in anndata object (copy=False), or in a csv (copy=True). Default is copy=False.')
    return parser

def save_het_scores_and_times_as_csv(het_scores, time_elapsed):
    for i in het_scores:
        df=pd.DataFrame(het_scores[i])
        df['overall_time_taken_allHetScores']=time_elapsed[i]
        df.to_csv(f'{i}.csv')
        


if __name__ == '__main__':
    args = _get_parser().parse_args()
    patients_good, patients_bad, het_scores, time_elapsed=run_shout.run(args.adata_path, args.cluster_key, args.radii, args.normalize, args.num_cell_types, args.coord_type, args.copy)
    
    if args.copy==True:
        save_het_scores_and_times_as_csv(het_scores, time_elapsed)
    
    print('Heterogeneity scores saved successfully!')


















