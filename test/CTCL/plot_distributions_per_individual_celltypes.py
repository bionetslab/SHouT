
import pandas as pd
import scipy.stats as stats
import os
import itertools as itt
from IPython.display import set_matplotlib_formats
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
import itertools
from statannotations.Annotator import Annotator
from decimal import Decimal


if __name__ == '__main__':

    p_values_cell_type=pd.read_csv(os.path.join('results', 'p_values_cell_type.csv'))
    cell_results=pd.read_csv(os.path.join('results', 'cell_results.csv'))
    sample_results=pd.read_csv(os.path.join('results', 'sample_results.csv'))
    celltypes_=list(np.unique(p_values_cell_type['cell_type']))
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    conditions_=list(np.unique(cell_results.condition))
    conditions_abbreviations_dict=dict(zip(conditions_, ['AD', 'PSO', 'CTCL']))
    disease_combinations=list(itertools.combinations(conditions_, 2))
    radii=[1, 2, 3, 4, 5]
    
    for i in celltypes_:
        for heterogeneity_measure in heterogeneity_measures:
            if heterogeneity_measure=='egophily':
                local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
                var_name=f'{heterogeneity_measure}_measure'
                value_name=f'{heterogeneity_measure}_score'
            else:
                local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
                var_name=f'local_{heterogeneity_measure}_measure'
                value_name=f'local_{heterogeneity_measure}_score'
            
            columns=local_heterogeneity_measure+['condition']
            data_=cell_results[cell_results['cell_type']==i][columns]
            data=data_.melt('condition', var_name=var_name, value_name=value_name)
            args = dict(x=var_name, y=value_name, data=data, hue="condition", hue_order=list(conditions_), order=local_heterogeneity_measure)
            pairs=[]
            pvals=[]
            for j in local_heterogeneity_measure:
                for k in disease_combinations:
                    pairs.append(((j, k[0]), (j, k[1])))
                    
                    if len(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'])==0:
                        p_=str(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==k[1]][p_values_cell_type['condition_2']==k[0]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
                    else:
                        p_=str(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'].values[0])
                    if float(p_)<0.05:
                        pvals.append('%.2E' % Decimal(p_))
                    else:
                        pvals.append('$X$')
            
            pairs_pvals_dict=dict(zip(pairs, pvals))
            
            fig, axes = plt.subplots(figsize=(20,10))
            sns.set(font_scale = 1.2)
            sns.set_style("white")
            # ax = sns.boxplot(**args)
            ax = sns.violinplot(**args, cut=0)
            annot = Annotator(ax, pairs, **args)
            annot.set_custom_annotations(pvals)
            annot.annotate()
            
            pairs_corrected=[]
            pvals_corrected=[]
            for j in range(len(annot.__dict__['annotations'])):
                pair_1=annot.__dict__['annotations'][j].__dict__['structs'][0]['group']
                pair_2=annot.__dict__['annotations'][j].__dict__['structs'][1]['group']   
                pairs_corrected.append((pair_1, pair_2))
                pvals_corrected.append(pairs_pvals_dict[(pair_1, pair_2)])
            
            fig, axes = plt.subplots(figsize=(20,10))
            sns.set(font_scale = 1.2)
            sns.set_style("white")
            # ax = sns.boxplot(**args)
            ax = sns.violinplot(**args, cut=0)
            
            text = "$X$: non-significant"
            txt=plt.text(.885, .7, text, ha='left', va='top', transform=plt.gcf().transFigure, fontsize=22, bbox=dict(boxstyle='square', fc='0.9', ec='0.9'))
            ax.set_xticklabels(radii, size=20)
            
            # if heterogeneity_measure=='entropy':
            #     ax.axhline(y=np.mean(sample_results.global_entropy), color='red')
            # elif heterogeneity_measure=='homophily':
            #     ax.axhline(y=np.mean(sample_results.global_homophily), color='red')
            
            annot = Annotator(ax, pairs, **args)
            annot.configure(text_format='simple', loc='inside', verbose=2, fontsize=25)
            annot.set_custom_annotations(pvals_corrected)
            annot.annotate()
            plt.xlabel('radius', fontsize=25, labelpad=20)
            plt.ylabel(f'{heterogeneity_measure} score', fontsize=25, labelpad=20)
            plt.title(f'Local {heterogeneity_measure} ({i})', fontsize=25, pad=20)
            plt.savefig(f'{i}_{heterogeneity_measure}_with_varying_radius.jpg', format='jpg', bbox_inches='tight')
            
            plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
            plt.setp(ax.get_legend().get_texts(), fontsize='25') # for legend text
            plt.setp(ax.get_legend().get_title(), fontsize='25') # for legend title
            new_labels=list(conditions_abbreviations_dict.values())
            for t, l in zip(ax.legend_.texts, new_labels):
                t.set_text(l)
            xticks=[]
            for j in ax.get_yticks():
                xticks.append(round(j,1))
            ax.set_yticklabels(xticks, size = 20)
            fig.tight_layout()
            plt.show()
    
    