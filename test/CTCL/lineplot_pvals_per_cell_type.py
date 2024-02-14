
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
    
    for heterogeneity_measure in heterogeneity_measures:
        if heterogeneity_measure=='egophily':
            local_heterogeneity_measure=[f'{heterogeneity_measure}_{radius}' for radius in radii]
            var_name=f'{heterogeneity_measure}_measure'
            value_name=f'{heterogeneity_measure}_score'
        else:
            local_heterogeneity_measure=[f'local_{heterogeneity_measure}_{radius}' for radius in radii]
            var_name=f'local_{heterogeneity_measure}_measure'
            value_name=f'local_{heterogeneity_measure}_score'
                    
        for k in disease_combinations:
            DF=[]
            for i in celltypes_:
                df=pd.DataFrame()
                p=[]
                for j in local_heterogeneity_measure:
                    
                    if len(p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['score']==j][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'])==0:
                        p_=p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==k[1]][p_values_cell_type['condition_2']==k[0]][p_values_cell_type['score']==j]['p_value_adj'].values[0]
                    else:
                        p_=p_values_cell_type[p_values_cell_type['cell_type']==i][p_values_cell_type['condition_1']==k[0]][p_values_cell_type['condition_2']==k[1]][p_values_cell_type['score']==j]['p_value_adj'].values[0]
                    p.append(p_)
                df[i]=p
                DF.append(df)
            DF=pd.concat(DF, axis=1)
            DF['radii']=radii
            
            fig, axes = plt.subplots(figsize=(20,10))
            sns.set(font_scale = 1.2)
            sns.set_style("white")
            sns.color_palette("tab10")
            ax=sns.lineplot(data=DF[celltypes_], markers=True, dashes=False, lw=1)
            sns.color_palette("tab10")
            
            for i in range(int(len(ax.lines.__dict__['_axes'].__dict__['_children'])/2)):
                ax.lines[i].set_markersize(25)
            
            plt.yscale('log')
            ax.axhline(y=0.05, color='black', linestyle='dashed')
           
            plt.xlabel('radius', fontsize=25, labelpad=20)
            plt.ylabel('p-values', fontsize=25, labelpad=20)
            if heterogeneity_measure=='egophily':
                plt.title(f'{heterogeneity_measure} ({k[0]} vs {k[1]})', fontsize=25, pad=20)
            else:
                plt.title(f'local {heterogeneity_measure} ({k[0]} vs {k[1]})', fontsize=25, pad=20)
            
            plt.savefig(f'lineplot_pvals_{heterogeneity_measure}_{k[0]}vs{k[1]}.jpg', format='jpg', bbox_inches='tight')
                
            plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=25, markerscale=3)
            plt.setp(ax.get_legend().get_texts(), fontsize='25') # for legend text
            plt.setp(ax.get_legend().get_title(), fontsize='25') # for legend title
            # ax.legend(fontsize=25)
            yticks=[]
            for j in ax.get_yticks():
                yticks.append(round(j,1))
            ax.set_yticklabels(yticks, size = 20)
            xticks=[]
            for j in ax.get_xticks():
                xticks.append(round(j,1)+1)
            ax.set_xticklabels(xticks, size = 20)
            ax.set_xticklabels(xticks, size = 20)
            fig.tight_layout()
            plt.show()
    
    