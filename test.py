#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:00:52 2024

@author: surya
"""

import numpy as np
import pandas as pd
import pickle
import squidpy as sq
import seaborn as sns
import shout
from IPython.display import set_matplotlib_formats
import random
import scanpy as sc
import time
import pickle
import matplotlib.pyplot as plt
import itertools
# ---
from statannotations.Annotator import Annotator
# ---
with open('MELC_TCL_caseStudy_notebook_data2.pkl', 'rb') as f:
    DF_reset, celltypes_, conditions_, DF_celltypes_dict, DF_celltypeIndices_dict, DF_conditions_dict, DF_conditionIndices_dict, condition_celltype_indices_dict, condition_celltype_DF_dict, condition_celltype_dict_dict, mwu_pvals_per_celltype, local_heterogeneity_measures, condition_combinations, mwu_pvals_per_celltype_insignificant_index, mwu_pvals_per_celltype_significant=pickle.load(f)
# ---
def createList(r1, r2):
 if (r1 == r2):
  return r1
 else:
  res = []
  while(r1 < r2+1 ):
   res.append(r1)
   r1 += 1
  return res
# ---
from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('retina', quality=100)
import seaborn as sns
# sns.set_theme()
sns.set_theme(style="whitegrid")
# plt.rcdefaults()
import matplotlib.pyplot as plt
# ---
def annot_stat(star, x1, x2, y, h, col='k', ax=None):
    ax = plt.gca() if ax is None else ax
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, star, ha='center', va='bottom', color=col)
# ---
def plot_local_heterogeneity_based_differential_analyses_conditionPairWise(condition1, condition2, heterogeneity_score_under_review):
    conditions_under_review=[condition1, condition2]
    cond1, cond2=sorted(conditions_under_review)
    # =====================
    count = 0
    _list_=[]
    _p_=[]
    not_present_list=[]
    for j in mwu_pvals_per_celltype_significant:
        if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
            count+=1
            _list_.append(j)
            _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
        else:
            not_present_list.append(j)
    # ---
    _dict_=dict(zip(_list_, _p_))
    _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
    _list_=_dict_.keys()
    _len_=len(_list_)
    no_of_rows=int(count/2)-1
    if no_of_rows==0:
        no_of_rows=1
    # =====================
    if len(not_present_list)>0:
        str_='['
        for i in not_present_list:
            str_=str_+str(i)+', '
        str_ = str_[:-2 or None]
        str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
    # =====================
    fig = plt.figure(layout="constrained", figsize=(15, 15))
    _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
    fig.suptitle(_title_, fontsize=12)
    ax_array = fig.subplots(no_of_rows,3, squeeze=False)
    ax_arrays=[]
    # ---
    count=-1
    count_outer=-1
    count_inner=-1
    for i in _list_:
        DF=[]
        _DF1_=pd.DataFrame()
        _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
        _DF1_['celltype']=i
        _DF1_['condition']=cond1
        DF.append(_DF1_)
        _DF2_=pd.DataFrame()
        _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
        _DF2_['celltype']=i
        _DF2_['condition']=cond2
        DF.append(_DF2_)
        DF=pd.concat(DF, axis=0)
        # ---
        if((float(count+1)%float(3))==0):
            count_outer+=1
            count_inner=-1
        # ---
        count+=1
        count_inner+=1
        print(count, count_outer, count_inner)
        # ---
        xlabel='\nCondition'
        _str_=f'n({i})'
        ylabel='Neighborhood enrichment score\n'
        # ---
        sns.boxplot(x="celltype", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
        # ---
        _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
        title=f'{i}('+'${{p}_{MWU}=}$'+f'{_x_}'+')'
        ax_array[count_outer, count_inner].set_title(title)
        # ---
        ax_array[count_outer, count_inner].set_xlabel('')
        # axes.set_ylabel(ylabel,fontsize=50)
        ax_array[count_outer, count_inner].set_ylabel('')
        # ---
        handles, labels = ax_array[count_outer, count_inner].get_legend_handles_labels()
        ax_array[count_outer, count_inner].legend(handles=handles[::-1], labels=[cond1, cond2])
        # ---
        ax_arrays.append(ax_array[count_outer, count_inner])
        # ==================================
    fig.tight_layout()
    plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
    plt.show()
    plt.close()
    # ---
# ---
condition_celltype_DF_dict

DF_local_heterogeneity_scores=[]
for i in condition_celltype_DF_dict:
    df_=[]
    for j in condition_celltype_DF_dict[i]:
        df_.append(condition_celltype_DF_dict[i][j])
    df_=pd.concat(df_, axis=0)
    df_['condition']=i
    DF_local_heterogeneity_scores.append(df_)
DF_local_heterogeneity_scores=pd.concat(DF_local_heterogeneity_scores, axis=0)

DF_local_heterogeneity_scores
# ---


#### LOCAL_ENTROPY_2: #####
# =========================================================
heterogeneity_score_under_review='local_entropy_2'
# ----------------------------------------------------
# ---
df=DF_local_heterogeneity_scores.copy()
# ---
celltypes=list(np.unique(df.celltype))
print(celltypes)
no_of_celltypes=len(celltypes)
no_of_rows=int(no_of_celltypes/2)-1
if no_of_rows==0:
    no_of_rows=1
# ---

    
# ---
fig = plt.figure(layout="constrained", figsize=(25, 25))
# # _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# # fig.suptitle(_title_, fontsize=12)
ax_array = fig.subplots(no_of_rows,3, squeeze=False)
ax_arrays=[]
# ---
count=-1
count_outer=-1
count_inner=-1
for i in celltypes:
    print(i)
    df_=df[df['celltype']==i]
    # ==================================
    pvals_=[]
    for j in condition_combinations:
        pvals_.append(str(mwu_pvals_per_celltype[i][j][heterogeneity_score_under_review]))
    condComb_pval_dict=dict(zip(condition_combinations, pvals_))
    # ==================================
    pairs=condition_combinations
    r1, r2 = 0, len(conditions_)-1
    conditions_idx=createList(r1, r2)
    conditions_idx_dict=dict(zip(conditions_, conditions_idx))
    condition_col=list(df_['condition'])
    _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
    df_['_idx_']=_idx_
    # ==================================
    # ---
    if((float(count+1)%float(3))==0):
        count_outer+=1
        count_inner=-1
    # ---
    count+=1
    count_inner+=1
    print(count, count_outer, count_inner)
    # ---
    xlabel='\nCondition'
    _str_=f'n({i})'
    ylabel='Neighborhood enrichment score\n'
    # ---
    ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=df_, ax=ax_array[count_outer, count_inner])
    # ax=sns.boxplot(x="Number of Research Years", y="Total Publications", hue="Number of Research Years", data=zero, ax=ax_array[count_outer, count_inner])
    # ---
    # from statannotations.Annotator import Annotator
    # annot = Annotator(ax, pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
    # annot.new_plot(ax, pairs=pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
    # annot.configure(loc='outside', comparisons_correction=None,)
    # annot.set_custom_annotations(pvals_)
    # annot.annotate()
    # ---
    annot_y_list=[1.40, 1.20, 1.00]
    annot_indexing=[(0,1), (0,2), (1,2)]
    count_annot=0
    idx_combinations=list(itertools.combinations(conditions_idx, 2))
    for k in idx_combinations:
        count_annot+=1
        if count_annot==4:
            count_annot=1
        star_=str(round(float(condComb_pval_dict[(conditions_[k[0]], conditions_[k[1]])]),5))
        annot_stat(star_, annot_indexing[count_annot-1][0], annot_indexing[count_annot-1][1], annot_y_list[count_annot-1], 0.05, ax=ax)
        print(conditions_[k[0]], conditions_[k[1]])
    # # ---
    # annot_stat('*', 0, 1, 1.40, 0.05, ax=ax)
    # annot_stat('**', 0, 2, 1.20, 0.05, ax=ax)
    # annot_stat('***', 1, 2, 1, 0.05, ax=ax)
    # # ---
    
# =============================================================================
#     for i in range(len(pairs)):
#         if i==2:
#             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 2, 1, ax=ax)
#             # annot_stat(pvals_[i], 0, 1, 1, 1, ax=ax)
#             annot_stat(pvals_[i], i, i+1, 1, 0.5, ax=ax)
# 
#         else:
#             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 1, 1, ax=ax)
#             annot_stat(pvals_[i], i, i+1, 0, 0.5, ax=ax)
# =============================================================================
    # ---
    # # _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
    # # title=f'{i}('+'${{p}_{MWU}=}$'+f'{_x_}'+')'
    # # ax_array[count_outer, count_inner].set_title(title)
    # ---
    _title_=f'\n\n\n\n{i}'
    ax_array[count_outer, count_inner].set_title(_title_)
    # ---
    ax_array[count_outer, count_inner].set_xlabel('')
    ax_array[count_outer, count_inner].set_ylabel('')
    # ---
    ax.set_xticklabels(list(conditions_))
    # ---
    ax.legend(title='Condition', loc='upper right')
    # ---
    handles, labels = ax_array[count_outer, count_inner].get_legend_handles_labels()
    # # ax_array[count_outer, count_inner].legend(handles=handles[::-1], labels=[cond1, cond2])
    # ---
    ax_arrays.append(ax_array[count_outer, count_inner])
    # ---
    # sns.move_legend(ax.figure.axes[0], "upper left", bbox_to_anchor=(1.25, 1))
    # ================================== 
fig.tight_layout()
plt.suptitle(f'{heterogeneity_score_under_review}\n\n\n')
plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_allCellTypes_allConditions.jpg', format='jpg', bbox_inches='tight')
plt.show()
plt.close()    
# ---
# # plot_local_heterogeneity_based_differential_analyses_cumulative('local_entropy_2')

# # =================================================================================
# # =================================================================================
# # =================================================================================


# ##### LOCAL_HOMOPHILY_2: #####
# # =========================================================
# heterogeneity_score_under_review='local_homophily_2'
# # ----------------------------------------------------
# # ---
# df=DF_local_heterogeneity_scores.copy()
# # ---
# celltypes=list(np.unique(df.celltype))
# print(celltypes)
# no_of_celltypes=len(celltypes)
# no_of_rows=int(no_of_celltypes/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # ---

    
# # ---
# fig = plt.figure(layout="constrained", figsize=(20, 20))
# # # _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# # # fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in celltypes:
#     print(i)
#     df_=df[df['celltype']==i]
#     # ==================================
#     pvals_=[]
#     for j in condition_combinations:
#         pvals_.append(str(mwu_pvals_per_celltype[i][j][heterogeneity_score_under_review]))
#     condComb_pval_dict=dict(zip(condition_combinations, pvals_))
#     # ==================================
#     pairs=condition_combinations
#     r1, r2 = 0, len(conditions_)-1
#     conditions_idx=createList(r1, r2)
#     conditions_idx_dict=dict(zip(conditions_, conditions_idx))
#     condition_col=list(df_['condition'])
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     df_['_idx_']=_idx_
#     # ==================================
#     # ---
#     if((float(count+1)%float(3))==0):
#         count_outer+=1
#         count_inner=-1
#     # ---
#     count+=1
#     count_inner+=1
#     print(count, count_outer, count_inner)
#     # ---
#     xlabel='\nCondition'
#     _str_=f'n({i})'
#     ylabel='Neighborhood enrichment score\n'
#     # ---
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=df_, ax=ax_array[count_outer, count_inner])
#     # ax=sns.boxplot(x="Number of Research Years", y="Total Publications", hue="Number of Research Years", data=zero, ax=ax_array[count_outer, count_inner])
#     # ---
#     # from statannotations.Annotator import Annotator
#     # annot = Annotator(ax, pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
#     # annot.new_plot(ax, pairs=pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
#     # annot.configure(loc='outside', comparisons_correction=None,)
#     # annot.set_custom_annotations(pvals_)
#     # annot.annotate()
#     # ---
#     annot_y_list=[0.14, 0.12, 0.10]
#     annot_indexing=[(0,1), (0,2), (1,2)]
#     count_annot=0
#     idx_combinations=list(itertools.combinations(conditions_idx, 2))
#     for k in idx_combinations:
#         count_annot+=1
#         if count_annot==4:
#             count_annot=1
#         star_=str(round(float(condComb_pval_dict[(conditions_[k[0]], conditions_[k[1]])]),5))
#         annot_stat(star_, annot_indexing[count_annot-1][0], annot_indexing[count_annot-1][1], annot_y_list[count_annot-1], 0.001, ax=ax)
#         print(conditions_[k[0]], conditions_[k[1]])
#     # # ---
#     # annot_stat('*', 0, 1, 1.40, 0.05, ax=ax)
#     # annot_stat('**', 0, 2, 1.20, 0.05, ax=ax)
#     # annot_stat('***', 1, 2, 1, 0.05, ax=ax)
#     # # ---
    
# # =============================================================================
# #     for i in range(len(pairs)):
# #         if i==2:
# #             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 2, 1, ax=ax)
# #             # annot_stat(pvals_[i], 0, 1, 1, 1, ax=ax)
# #             annot_stat(pvals_[i], i, i+1, 1, 0.5, ax=ax)
# # 
# #         else:
# #             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 1, 1, ax=ax)
# #             annot_stat(pvals_[i], i, i+1, 0, 0.5, ax=ax)
# # =============================================================================
#     # ---
#     # # _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # # title=f'{i}('+'${{p}_{MWU}=}$'+f'{_x_}'+')'
#     # # ax_array[count_outer, count_inner].set_title(title)
#     # ---
#     _title_=f'\n\n\n\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels(list(conditions_))
#     # ---
#     ax.legend(title='Condition', loc='upper right')
#     # ---
#     handles, labels = ax_array[count_outer, count_inner].get_legend_handles_labels()
#     # # ax_array[count_outer, count_inner].legend(handles=handles[::-1], labels=[cond1, cond2])
#     # ---
#     ax_arrays.append(ax_array[count_outer, count_inner])
#     # ---
#     # sns.move_legend(ax.figure.axes[0], "upper left", bbox_to_anchor=(1.25, 1))
#     # ================================== 
# fig.tight_layout()
# plt.suptitle(f'{heterogeneity_score_under_review}\n\n\n')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_allCellTypes_allConditions.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()    
# # ---
# # # plot_local_heterogeneity_based_differential_analyses_cumulative('local_entropy_2')

# # =================================================================================
# # =================================================================================
# # =================================================================================



# ##### EGOPHILY_2: #####
# # =========================================================
# heterogeneity_score_under_review='egophily_2'
# # ----------------------------------------------------
# # ---
# df=DF_local_heterogeneity_scores.copy()
# # ---
# celltypes=list(np.unique(df.celltype))
# print(celltypes)
# no_of_celltypes=len(celltypes)
# no_of_rows=int(no_of_celltypes/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # ---

    
# # ---
# fig = plt.figure(layout="constrained", figsize=(20, 20))
# # # _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# # # fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in celltypes:
#     print(i)
#     df_=df[df['celltype']==i]
#     # ==================================
#     pvals_=[]
#     for j in condition_combinations:
#         pvals_.append(str(mwu_pvals_per_celltype[i][j][heterogeneity_score_under_review]))
#     condComb_pval_dict=dict(zip(condition_combinations, pvals_))
#     # ==================================
#     pairs=condition_combinations
#     r1, r2 = 0, len(conditions_)-1
#     conditions_idx=createList(r1, r2)
#     conditions_idx_dict=dict(zip(conditions_, conditions_idx))
#     condition_col=list(df_['condition'])
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     df_['_idx_']=_idx_
#     # ==================================
#     # ---
#     if((float(count+1)%float(3))==0):
#         count_outer+=1
#         count_inner=-1
#     # ---
#     count+=1
#     count_inner+=1
#     print(count, count_outer, count_inner)
#     # ---
#     xlabel='\nCondition'
#     _str_=f'n({i})'
#     ylabel='Neighborhood enrichment score\n'
#     # ---
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=df_, ax=ax_array[count_outer, count_inner])
#     # ax=sns.boxplot(x="Number of Research Years", y="Total Publications", hue="Number of Research Years", data=zero, ax=ax_array[count_outer, count_inner])
#     # ---
#     # from statannotations.Annotator import Annotator
#     # annot = Annotator(ax, pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
#     # annot.new_plot(ax, pairs=pairs, data=df_, x="condition", y=heterogeneity_score_under_review)
#     # annot.configure(loc='outside', comparisons_correction=None,)
#     # annot.set_custom_annotations(pvals_)
#     # annot.annotate()
#     # ---
#     annot_y_list=[0.90, 0.80, 0.70]
#     annot_indexing=[(0,1), (0,2), (1,2)]
#     count_annot=0
#     idx_combinations=list(itertools.combinations(conditions_idx, 2))
#     for k in idx_combinations:
#         count_annot+=1
#         if count_annot==4:
#             count_annot=1
#         star_=str(round(float(condComb_pval_dict[(conditions_[k[0]], conditions_[k[1]])]),5))
#         annot_stat(star_, annot_indexing[count_annot-1][0], annot_indexing[count_annot-1][1], annot_y_list[count_annot-1], 0.02, ax=ax)
#         print(conditions_[k[0]], conditions_[k[1]])
#     # # ---
#     # annot_stat('*', 0, 1, 1.40, 0.05, ax=ax)
#     # annot_stat('**', 0, 2, 1.20, 0.05, ax=ax)
#     # annot_stat('***', 1, 2, 1, 0.05, ax=ax)
#     # # ---
    
# # =============================================================================
# #     for i in range(len(pairs)):
# #         if i==2:
# #             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 2, 1, ax=ax)
# #             # annot_stat(pvals_[i], 0, 1, 1, 1, ax=ax)
# #             annot_stat(pvals_[i], i, i+1, 1, 0.5, ax=ax)
# # 
# #         else:
# #             # annot_stat(pvals_[i], pairs[i][0], pairs[i][1], 1, 1, ax=ax)
# #             annot_stat(pvals_[i], i, i+1, 0, 0.5, ax=ax)
# # =============================================================================
#     # ---
#     # # _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # # title=f'{i}('+'${{p}_{MWU}=}$'+f'{_x_}'+')'
#     # # ax_array[count_outer, count_inner].set_title(title)
#     # ---
#     _title_=f'\n\n\n\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels(list(conditions_))
#     # ---
#     ax.legend(title='Condition', loc='upper right')
#     # ---
#     handles, labels = ax_array[count_outer, count_inner].get_legend_handles_labels()
#     # # ax_array[count_outer, count_inner].legend(handles=handles[::-1], labels=[cond1, cond2])
#     # ---
#     ax_arrays.append(ax_array[count_outer, count_inner])
#     # ---
#     # sns.move_legend(ax.figure.axes[0], "upper left", bbox_to_anchor=(1.25, 1))
#     # ================================== 
# fig.tight_layout()
# plt.suptitle(f'{heterogeneity_score_under_review}\n\n\n')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_allCellTypes_allConditions.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()    
# # ---
# # # plot_local_heterogeneity_based_differential_analyses_cumulative('local_entropy_2')

# # =================================================================================
































