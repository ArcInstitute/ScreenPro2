## Copyright (c) 2022-2024 ScreenPro2 Development Team.
## All rights reserved.
## Gilbart Lab, UCSF / Arc Institute.
## Multi-Omics Tech Center, Arc Insititue.
##
## courtesy Tyler Fair (@tdfair), M. Horlbeck, (@mhorlbeck)

'''
evaluate module: evaluate essentiality detection performance
'''

import numpy as np
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
import matplotlib.lines as mlines


def calcROC(df_in, essential, nonessential, score_col, target_col='target', verbose=False):
    df = df_in.copy()
    df[target_col] = df[target_col].str.split('-').str[0]
    
    # AUC-ROC
    df['DepMap'] = np.nan
    df.loc[df[target_col].isin(essential),'DepMap'] = 'essential'
    df.loc[df[target_col].isin(nonessential),'DepMap'] = 'non_essential'
    
    y_true = df[df['DepMap'].notna()]['DepMap']
    y_scores = df[df['DepMap'].notna()][score_col]
    
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label='non_essential')
    
    if verbose: print('ROC-AUC score:', metrics.roc_auc_score(y_true, y_scores))
    
    return fpr, tpr


def calcPR(df_in, truePos, trueNeg, score_col, target_col='target', ascending=True, verbose=False):

    df = df_in.copy()

    scoreList = df.set_index(target_col)[score_col]

    truePos = truePos.intersection(scoreList.index)
    trueNeg = trueNeg.intersection(scoreList.index)
    
    cumulativeTup = [(0,1,np.nan)]
    cumulativeTP = 0.0
    cumulativeFP = 0.0
    
    tup_cross95 = None
    
    for gene, fold_change in scoreList.sort_values(inplace=False, ascending=ascending).items():
        if gene in truePos or gene in trueNeg:
            if gene in truePos:
                cumulativeTP += 1

            elif gene in trueNeg:
                cumulativeFP += 1

            cumulativeTup.append((cumulativeTP / len(truePos), cumulativeTP / (cumulativeTP + cumulativeFP), fold_change))
    
    tup_cross95 = [cross95 for cross95 in cumulativeTup[::-1] if cross95[1] >= 0.95][0]
    
    if verbose: print('Precision-Recall AUC:', tup_cross95)
    
    return cumulativeTup, tup_cross95
