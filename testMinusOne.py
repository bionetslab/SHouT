#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 02:55:23 2024

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






# LOCAL_ENTROPY_2:
# ---
condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Eczema', 'local_entropy_2']
# ========================================================
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
    condition_col=list(DF['condition'])
    conditions_idx_dict={cond1: 0, cond2: 1}
    _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
    DF['_idx_']=_idx_
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
    ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
    # ---
    _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
    # ---
    # =====
    annot_stat(str(_x_), 0, 1, 1.00, 0.05, ax=ax)
    # =====
    # ---
    _title_=f'\n{i}'
    ax_array[count_outer, count_inner].set_title(_title_)
    # ---
    ax_array[count_outer, count_inner].set_xlabel('')
    ax_array[count_outer, count_inner].set_ylabel('')
    # ---
    ax.set_xticklabels([cond1, cond2])
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
plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
plt.show()
plt.close()








    

# # =================================================================================
# # =================================================================================
# # =================================================================================













# # LOCAL_HOMOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Eczema', 'local_homophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.1, 0.001, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()














# # EGOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Eczema', 'egophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.50, 0.05, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()


# # LOCAL_ENTROPY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Psoriasis', 'local_entropy_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 1.00, 0.05, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()



# # LOCAL_HOMOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Psoriasis', 'local_homophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.05, 0.005, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()







# # EGOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['T-Cell Lymphoma', 'Psoriasis', 'egophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.50, 0.05, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()



# # LOCAL_ENTROPY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['Eczema', 'Psoriasis', 'local_entropy_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 1.00, 0.05, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()





# # LOCAL_HOMOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['Eczema', 'Psoriasis', 'local_homophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.1, 0.005, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()





# # LOCAL_EGOPHILY_2:
# # ---
# condition1, condition2, heterogeneity_score_under_review=['Eczema', 'Psoriasis', 'egophily_2']
# # ========================================================
# conditions_under_review=[condition1, condition2]
# cond1, cond2=sorted(conditions_under_review)
# # =====================
# count = 0
# _list_=[]
# _p_=[]
# not_present_list=[]
# for j in mwu_pvals_per_celltype_significant:
#     if ((cond1, cond2) in mwu_pvals_per_celltype_significant[j]) and (heterogeneity_score_under_review in mwu_pvals_per_celltype_significant[j][(cond1, cond2)]):
#         count+=1
#         _list_.append(j)
#         _p_.append(mwu_pvals_per_celltype_significant[j][(cond1, cond2)][heterogeneity_score_under_review])
#     else:
#         not_present_list.append(j)
# # ---
# _dict_=dict(zip(_list_, _p_))
# _dict_=dict(sorted(_dict_.items(), key=lambda item: item[1]))
# _list_=_dict_.keys()
# _len_=len(_list_)
# no_of_rows=int(count/2)-1
# if no_of_rows==0:
#     no_of_rows=1
# # =====================
# if len(not_present_list)>0:
#     str_='['
#     for i in not_present_list:
#         str_=str_+str(i)+', '
#     str_ = str_[:-2 or None]
#     str_+='] cells do not show significant ' + f'{heterogeneity_score_under_review} scores between conditions {cond1} and {cond2}.'
# # =====================
# fig = plt.figure(layout="constrained", figsize=(15, 15))
# _title_=f'{heterogeneity_score_under_review} scores of significantly differentiated cell types (p<0.05)\n'+f'({cond1} vs {cond2})\n'
# fig.suptitle(_title_, fontsize=12)
# ax_array = fig.subplots(no_of_rows,3, squeeze=False)
# ax_arrays=[]
# # ---
# count=-1
# count_outer=-1
# count_inner=-1
# for i in _list_:
#     DF=[]
#     _DF1_=pd.DataFrame()
#     _DF1_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond1][i][heterogeneity_score_under_review].values()
#     _DF1_['celltype']=i
#     _DF1_['condition']=cond1
#     DF.append(_DF1_)
#     _DF2_=pd.DataFrame()
#     _DF2_[heterogeneity_score_under_review]=condition_celltype_dict_dict[cond2][i][heterogeneity_score_under_review].values()
#     _DF2_['celltype']=i
#     _DF2_['condition']=cond2
#     DF.append(_DF2_)
#     DF=pd.concat(DF, axis=0)
#     # ---
#     condition_col=list(DF['condition'])
#     conditions_idx_dict={cond1: 0, cond2: 1}
#     _idx_=[conditions_idx_dict.get(item,item)  for item in condition_col]
#     DF['_idx_']=_idx_
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
#     ax=sns.boxplot(x="_idx_", y=heterogeneity_score_under_review, hue="condition", data=DF, ax=ax_array[count_outer, count_inner])
#     # ---
#     _x_=round(mwu_pvals_per_celltype_significant[i][(cond1, cond2)][heterogeneity_score_under_review],5)
#     # ---
#     # =====
#     annot_stat(str(_x_), 0, 1, 0.7, 0.005, ax=ax)
#     # =====
#     # ---
#     _title_=f'\n{i}'
#     ax_array[count_outer, count_inner].set_title(_title_)
#     # ---
#     ax_array[count_outer, count_inner].set_xlabel('')
#     ax_array[count_outer, count_inner].set_ylabel('')
#     # ---
#     ax.set_xticklabels([cond1, cond2])
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
# plt.suptitle(f'\n{heterogeneity_score_under_review} scores - {cond1} vs {cond2}')
# plt.savefig(f'boxplot_significant_{heterogeneity_score_under_review}_scores_{cond1}vs{cond2}.jpg', format='jpg', bbox_inches='tight')
# plt.show()
# plt.close()









