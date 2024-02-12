
import os
import scanpy as sc
import pandas as pd
import numpy as np
import itertools
import matlab
import matlab.engine
myEngine = matlab.engine.start_matlab()
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from itertools import permutations
from statannotations.Annotator import Annotator

def plots_across_all_celltypes(DF, conditions_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_significant, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index):
    from IPython.display import set_matplotlib_formats
    import seaborn as sns
    sns.set_theme(style="whitegrid")
    import matplotlib.pyplot as plt
    
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    
    for heterogeneity_measure in heterogeneity_measures:
        if heterogeneity_measure=='entropy':
            local_heterogeneity_measure=['local_entropy_1', 'local_entropy_2', 'local_entropy_3', 'local_entropy_4', 'local_entropy_5']
            var_name='local_entropy_measure'
            value_name='local_entropy_score'
        elif heterogeneity_measure=='homophily':
            local_heterogeneity_measure=['local_homophily_1', 'local_homophily_2', 'local_homophily_3', 'local_homophily_4', 'local_homophily_5']
            var_name='local_homophily_measure'
            value_name='local_homophily_score'
        elif heterogeneity_measure=='egophily':
            local_heterogeneity_measure=['egophily_1', 'egophily_2', 'egophily_3', 'egophily_4', 'egophily_5']
            var_name='egophily_measure'
            value_name='egophily_score'
            
        columns=local_heterogeneity_measure.copy()
        columns.append('Group')
        DF_het=DF[columns]
        data = DF_het.melt('Group', var_name=var_name, value_name=value_name)
        fig, axes = plt.subplots(figsize=(15,15))
        sns.set(font_scale = 1.2)
        sns.set_style("white")
        args = dict(x=var_name, y=value_name, data=data, hue="Group", hue_order=list(conditions_), order=local_heterogeneity_measure)
        # ax = sns.boxplot(**args)
        ax = sns.violinplot(**args)
            
        if heterogeneity_measure=='entropy':
            ax.axhline(y=np.mean(DF.global_entropy), color='red')
        elif heterogeneity_measure=='homophily':
            ax.axhline(y=np.mean(DF.global_homophily), color='red')
            
        disease_combinations=list(itertools.combinations(conditions_, 2))
            
        pairs=[]
        pvals=[]
        for j in local_heterogeneity_measure:
            for k in disease_combinations:
                pairs.append([(j, k[0]), (j, k[1])])
                pvals.append(str('$p_{MWU}=$')+str(round(mwu_pvals_all_celltypes[(k[0], k[1])][j],2)))
            
        annot = Annotator(ax, pairs, **args)
        annot.configure(text_format='simple', loc='inside', verbose=2)
        annot.set_custom_annotations(pvals)
        annot.annotate()
        plt.xlabel(f'{heterogeneity_measure} measure (with varying radius)')
        plt.ylabel(f'{heterogeneity_measure} score')
        plt.title(f'All celltypes\n({heterogeneity_measure} scores across conditions with varying radius)')
        # plt.legend(title=f'Celltypes/ conditions', title_fontsize = 13, labels=list(conditions_))
        plt.savefig(f'allCelltypes_{heterogeneity_measure}_with_varying_radius.jpg', format='jpg', bbox_inches='tight')
        fig.tight_layout()
        plt.show()

def plots_per_celltype(DF, conditions_, celltypes_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_significant, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index):
    from IPython.display import set_matplotlib_formats
    import seaborn as sns
    sns.set_theme(style="whitegrid")
    import matplotlib.pyplot as plt
    
    heterogeneity_measures=['entropy', 'homophily', 'egophily']
    
    for i in celltypes_:
        for heterogeneity_measure in heterogeneity_measures:
            if heterogeneity_measure=='entropy':
                local_heterogeneity_measure=['local_entropy_1', 'local_entropy_2', 'local_entropy_3', 'local_entropy_4', 'local_entropy_5']
                var_name='local_entropy_measure'
                value_name='local_entropy_score'
            elif heterogeneity_measure=='homophily':
                local_heterogeneity_measure=['local_homophily_1', 'local_homophily_2', 'local_homophily_3', 'local_homophily_4', 'local_homophily_5']
                var_name='local_homophily_measure'
                value_name='local_homophily_score'
            elif heterogeneity_measure=='egophily':
                local_heterogeneity_measure=['egophily_1', 'egophily_2', 'egophily_3', 'egophily_4', 'egophily_5']
                var_name='egophily_measure'
                value_name='egophily_score'
            
            columns=local_heterogeneity_measure.copy()
            columns.append('Group')
            DF_=DF.copy()
            DF_=DF_[DF_['celltype']==i]
            DF_het=DF_[columns]
            data = DF_het.melt('Group', var_name=var_name, value_name=value_name)
            fig, axes = plt.subplots(figsize=(15,15))
            sns.set(font_scale = 1.2)
            sns.set_style("white")
            args = dict(x=var_name, y=value_name, data=data, hue="Group", hue_order=list(conditions_), order=local_heterogeneity_measure)
            # ax = sns.boxplot(**args)
            ax = sns.violinplot(**args)
            
            if heterogeneity_measure=='entropy':
                ax.axhline(y=np.mean(DF.global_entropy), color='red')
            elif heterogeneity_measure=='homophily':
                ax.axhline(y=np.mean(DF.global_homophily), color='red')
            
            disease_combinations=list(itertools.combinations(conditions_, 2))
            
            pairs=[]
            pvals=[]
            for j in local_heterogeneity_measure:
                for k in disease_combinations:
                    pairs.append([(j, k[0]), (j, k[1])])
                    pvals.append(str('$p_{MWU}=$')+str(round(mwu_pvals_per_celltype[i][(k[0], k[1])][j],2)))
            
            annot = Annotator(ax, pairs, **args)
            annot.configure(text_format='simple', loc='inside', verbose=2)
            annot.set_custom_annotations(pvals)
            annot.annotate()
            plt.xlabel(f'{heterogeneity_measure} measure (with varying radius)')
            plt.ylabel(f'{heterogeneity_measure} score')
            plt.title(f'{i}\n({heterogeneity_measure} scores across conditions with varying radius)')
            plt.savefig(f'{i}_{heterogeneity_measure}_with_varying_radius.jpg', format='jpg', bbox_inches='tight')
            fig.tight_layout()
            plt.show()
    plots_across_all_celltypes(DF, conditions_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_significant, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index)

def retain_significant_pvalues(DF, conditions_, celltypes_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_insignificant_index, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index):
    mwu_pvals_per_celltype_significant={}
    for i in mwu_pvals_per_celltype:
        mwu_pvals_per_celltype_significant[i]={}
        for j in mwu_pvals_per_celltype[i]:
            mwu_pvals_per_celltype_significant[i][j]={}
            for k in mwu_pvals_per_celltype[i][j]:
                if not([i, j, k] in mwu_pvals_per_celltype_insignificant_index):
                    mwu_pvals_per_celltype_significant[i][j][k]=mwu_pvals_per_celltype[i][j][k]
    plots_per_celltype(DF, conditions_, celltypes_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_significant, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index)

def mwu_computation(DF, celltypes_, conditions_, condition_celltype_DF_dict, condition_celltype_dict_dict):
    local_heterogeneity_measures=['local_entropy_1','local_entropy_2','local_entropy_3','local_entropy_4','local_entropy_5',
                                  'local_homophily_1','local_homophily_2','local_homophily_3','local_homophily_4','local_homophily_5',
                                  'egophily_1','egophily_2', 'egophily_3', 'egophily_4', 'egophily_5']
    condition_combinations=list(itertools.combinations(conditions_, 2))
    # ======================== p-values across all celltypes: ===============================
    mwu_pvals_all_celltypes={}
    mwu_pvals_all_celltypes_insignificant_index=[]
    for i1, i2 in condition_combinations:
        mwu_pvals_all_celltypes[(i1, i2)]={}
        
        for k in local_heterogeneity_measures:
            print('i1, i2:'+str(i1)+str(', ')+str(i2)+'k:'+str(k))
            l1=list(DF[DF['Group']==i1][k])
            matList1=matlab.double(l1)
            l2=list(DF[DF['Group']==i2][k])
            matList2=matlab.double(l2)
            p=myEngine.ranksum(matList1,matList2)
            mwu_pvals_all_celltypes[(i1, i2)][k]=p
            if p>=0.05:
                mwu_pvals_all_celltypes_insignificant_index.append([(i1, i2), k])
    
    # ======================= p-values across individual celltypes: ================================
    
    mwu_pvals_per_celltype={}
    mwu_pvals_per_celltype_insignificant_index=[]
    for j in celltypes_:
        mwu_pvals_per_celltype[j]={}
        for i1, i2 in condition_combinations:
            mwu_pvals_per_celltype[j][(i1, i2)]={}
            for k in local_heterogeneity_measures:
                print('j:'+j+', '+'i1, i2:'+str(i1)+str(', ')+str(i2)+'k:'+str(k))
                l1=list(condition_celltype_dict_dict[i1][j][k].values())
                matList1=matlab.double(l1)
                l2=list(condition_celltype_dict_dict[i2][j][k].values())
                matList2=matlab.double(l2)
                p=myEngine.ranksum(matList1,matList2)
                mwu_pvals_per_celltype[j][(i1, i2)][k]=p
                if p>=0.05:
                    mwu_pvals_per_celltype_insignificant_index.append([j, (i1, i2), k])
    retain_significant_pvalues(DF, conditions_, celltypes_, mwu_pvals_per_celltype, mwu_pvals_per_celltype_insignificant_index, mwu_pvals_all_celltypes, mwu_pvals_all_celltypes_insignificant_index)

def data_preparation_for_mwu_computation(DF):
    DF_reset=DF.copy()
    DF_reset=DF_reset.reset_index(drop=True)
    celltypes_=np.unique(DF_reset.celltype)
    conditions_=np.unique(DF_reset.Group)
    DF_celltypes_dict={}
    DF_celltypeIndices_dict={}
    for i in celltypes_:
        df_=DF_reset[DF_reset['celltype']==i]
        DF_celltypes_dict[i]=df_
        DF_celltypeIndices_dict[i]=list(df_.index)
    DF_conditions_dict={}
    DF_conditionIndices_dict={}
    for i in conditions_:
        df_=DF_reset[DF_reset['Group']==i]
        DF_conditions_dict[i]=df_
        DF_conditionIndices_dict[i]=list(df_.index)
    condition_celltype_indices_dict={}
    condition_celltype_DF_dict={}
    condition_celltype_dict_dict={}
    for i in conditions_:
        condition_celltype_indices_dict[i]={}
        condition_celltype_DF_dict[i]={}
        condition_celltype_dict_dict[i]={}
        for j in celltypes_:
            indices_=list(set(DF_conditionIndices_dict[i]).intersection(set(DF_celltypeIndices_dict[j])))
            condition_celltype_indices_dict[i][j]=indices_
            condition_celltype_DF_dict[i][j]=DF_reset.loc[indices_]
            condition_celltype_dict_dict[i][j]=DF_reset.loc[indices_].copy().to_dict()
    mwu_computation(DF, celltypes_, conditions_, condition_celltype_DF_dict, condition_celltype_dict_dict)

def convert_dict_to_df(all_samples):
    LOCAL_ENTROPY_2={}
    LOCAL_HOMOPHILY_2={}
    EGOPHILY_2={}
    DF=[]
    for sample in all_samples:
        df=pd.DataFrame()
        celltype=all_samples[sample].obs['celltype']
        df['celltype']=celltype
        local_entropy_1=all_samples[sample].obs['local_entropy_1']
        df['local_entropy_1']=local_entropy_1
        local_entropy_2=all_samples[sample].obs['local_entropy_2']
        df['local_entropy_2']=local_entropy_2
        local_entropy_3=all_samples[sample].obs['local_entropy_3']
        df['local_entropy_3']=local_entropy_3
        local_entropy_4=all_samples[sample].obs['local_entropy_4']
        df['local_entropy_4']=local_entropy_4
        local_entropy_5=all_samples[sample].obs['local_entropy_5']
        df['local_entropy_5']=local_entropy_5
        local_homophily_1=all_samples[sample].obs['local_homophily_1']
        df['local_homophily_1']=local_homophily_1
        local_homophily_2=all_samples[sample].obs['local_homophily_2']
        df['local_homophily_2']=local_homophily_2
        local_homophily_3=all_samples[sample].obs['local_homophily_3']
        df['local_homophily_3']=local_homophily_3
        local_homophily_4=all_samples[sample].obs['local_homophily_4']
        df['local_homophily_4']=local_homophily_4
        local_homophily_5=all_samples[sample].obs['local_homophily_5']
        df['local_homophily_5']=local_homophily_5
        egophily_1=all_samples[sample].obs['egophily_1']
        df['egophily_1']=egophily_1
        egophily_2=all_samples[sample].obs['egophily_2']
        df['egophily_2']=egophily_2
        egophily_3=all_samples[sample].obs['egophily_3']
        df['egophily_3']=egophily_3
        egophily_4=all_samples[sample].obs['egophily_4']
        df['egophily_4']=egophily_4
        egophily_5=all_samples[sample].obs['egophily_5']
        df['egophily_5']=egophily_5
        df['cell_id']=all_samples[sample].obs['cell_id']
        df['celltype']=all_samples[sample].obs['celltype']
        df['cell_index_in_patient']=all_samples[sample].to_df().index
        len_=len(all_samples[sample].to_df().index)
        df['Group']=[all_samples[sample].uns['Group'][0]]*len_
        df['SHouT_execution_time']=[all_samples[sample].uns['SHouT_execution_time']]*len_
        df['patient_id']=[all_samples[sample].uns['patient_id']]*len_
        df['patient_index']=[sample]*len_
        df['global_entropy']=[all_samples[sample].uns['global_entropy']]*len_
        df['global_homophily']=[all_samples[sample].uns['global_homophily']]*len_
        DF.append(df)
    DF=pd.concat(DF, axis=0)
    conditions=np.unique(DF.Group)
    data_preparation_for_mwu_computation(DF)

def read_h5ad_to_dict(het_score_path):
    all_samples={} # Make empty dict.
    for file in os.listdir(het_score_path): # Go through all files in adata_path
        filename = os.fsdecode(file)
        if filename.endswith(".h5ad"): # Check to make sure we are only reading .h5ad files.
            all_samples[filename.split('.')[0]]=sc.read_h5ad(het_score_path+filename) # or, sc.read()
            continue
        else:
            continue
    convert_dict_to_df(all_samples)


def run(het_score_path):
    read_h5ad_to_dict(het_score_path)


