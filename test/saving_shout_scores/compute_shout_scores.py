import run_shout
import argparse
import pickle
import pandas as pd

def validate_list_of_radii(str_):
    radii_=[]
    for i in str_:
        try:
            radii_.append(int(i))
        except:
            print(f'{i} is not a valid value for radius, hence dropped.')
            continue
    if len(radii_)==0:
        return None
    else:
        _len_valid_=0
        _radii_valid_=[]
        for i in radii_:
            if i>0:
                _len_valid_+=1
                _radii_valid_.append(i)
        if len(_radii_valid_)==0:
            return None
        else:
            return _radii_valid_

def list_of_ints(arg):
    radii=validate_list_of_radii(arg.split(','))
    return radii


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_path', type=str, help='Please specify path to anndata pickle data.')
    parser.add_argument('cluster_key', type=str, help='Please specify cluster key.')
    parser.add_argument('radii', type=list_of_ints, help='Please specify list of radii within which to calculate local heterogeneity scores.')
    parser.add_argument('--normalize', type=bool, default=True, help='Optionally specify "normalize", default value is True.')
    parser.add_argument('--num_cell_types', type=int, default=None, help='Optionally specify no. of celltypes in sample.')
    parser.add_argument('--coord_type', type=str, choices=['generic', 'grid'], default='GENE_SYMBOL')
    parser.add_argument('--copy', type=bool, default=False, help='Optionally specify whether heterogenety scores to be written in anndata object (copy=False), or in a csv (copy=True). Default is copy=False.')
    return parser

def save_het_scores_and_times_as_csv(Het_scores, time_elapsed):
    for i in Het_scores:
        df=pd.DataFrame(Het_scores[i])
        df['overall_time_taken_allHetScores']=time_elapsed[i]
        df.to_csv(f'{i}.csv')
        


if __name__ == '__main__':
    args = _get_parser().parse_args()
    if args.radii==None:
        raise ValueError('List of radii needs to contain at least one positive integer.')
    else:
        print(args.radii)
        radii=(args.radii).copy()
    patients_good, patients_bad, Het_scores, time_elapsed=run_shout.run(args.adata_path, args.cluster_key, radii, args.normalize, args.num_cell_types, args.coord_type, args.copy)
    # ---
    if args.copy==True:
        save_het_scores_and_times_as_csv(Het_scores, time_elapsed)
    # ---
    print('Heterogeneity scores saved successfully!')


















