"""
delta module
"""
import numpy as np
import anndata as ad
import pandas as pd


### Key functions for calculating delta phenotype score

def buildPhenotypeData(adata, score_tag, cond_ref, cond_test, growth_rate_reps=None, ctrl_label='negative_control'):
    """Calculate phenotype score for each pair of replicates.
    """
    score_name = f'{score_tag}:{cond_test}_vs_{cond_ref}'

    adat = adata.copy()
    
    #TODO: fix `_calculateGrowthFactor` and `_getTreatmentDoublingRate` to maintain same format
    # growth_factor_table = self._calculateGrowthFactor(
    #     untreated = untreated, treated = x_label, db_rate_col = db_rate_col
    # )

    adat_ctrl = adat[:, adat.var.targetType.eq(ctrl_label)].copy()
    
    results = {}
    
    if growth_rate_reps is None:
        growth_rate_reps = dict([(replicate, 1) for replicate in adat.obs.replicate.unique()])
    
    for replicate in adat.obs.replicate.unique():
        x=adat[adat.obs.query(f'condition == "{cond_ref}" & replicate == {str(replicate)}').index].X.T
        y=adat[adat.obs.query(f'condition == "{cond_test}" & replicate == {str(replicate)}').index].X.T
    
        x_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{cond_ref}" & replicate == {str(replicate)}').index].X.T
        y_ctrl=adat_ctrl[adat_ctrl.obs.query(f'condition == "{cond_test}" & replicate == {str(replicate)}').index].X.T
    
        res = calculateDelta(
            x=x,y=y,x_ctrl=x_ctrl,y_ctrl=y_ctrl,
            growth_rate=growth_rate_reps[replicate],
        )
    
        results.update({f'{score_name}::replicate_{replicate}': res.reshape(-1)})
    
    out = pd.DataFrame(
        results,
        # obs = growth_factor_table.loc[pdata_df.index,:],
        index=adat.var.index
    )

    return out


def calculateDelta(x, y, x_ctrl, y_ctrl, growth_rate):
    """Calculate phenotype score normalized by negative control and growth rate.
    
    Args:
        x (np.array): array of values
        y (np.array): array of values
        x_ctrl (np.array): array of values
        y_ctrl (np.array): array of values
        growth_rate (int): growth rate
    
    Returns:
        np.array: array of scores
    """
    # calculate control median and std 
    ctrl_median = np.median(
        calculateLog2e(x=x_ctrl, y=y_ctrl), 
        axis=0 # for each individual sample (i.e. replicate)
    )

    # calculate log2e (i.e. log2 fold change enrichment y / x)
    log2e = calculateLog2e(x=x, y=y)

    # calculate delta score normalized by control median and growth rate
    delta = (log2e - ctrl_median) / growth_rate
    
    return delta


### Utility functions

def calculateLog2e(x, y):
    return np.log2(y) - np.log2(x)


def averageBestN(scores, numToAverage):
    # Sort and find top n guide per target, see #18
    return np.mean(sorted(scores, key=abs, reverse=True)[:numToAverage])


def aggregatePhenotypeScore(scores, level=None):
    # average delta score by given level (replicate or target)
    if level == None:
        pass
    elif level == 'replicate' or level == 'col':
        scores = np.mean(scores, axis=1)
    elif level == 'target' or level == 'row':
        scores = np.mean(scores, axis=0)
    elif level == 'all':
        scores = np.mean(scores)
    else:
        raise ValueError(f'Invalid level: {level}')

    return scores


def generatePseudoGeneAnnData(adata, num_pseudogenes='auto', pseudogene_size='auto', ctrl_label='negative_control'):
    """Generate pseudogenes from negative control elements in the library.

    Args:
        adata (AnnData): AnnData object
        num_pseudogenes (int): number of pseudogenes to generate
        pseudogene_size (int): number of sgRNA elements in each pseudogene
        ctrl_label (str): control label, default is 'negative_control'
    
    Returns:
        AnnData: AnnData object with pseudogenes
    """
    if pseudogene_size == 'auto':
        # sgRNA elements / target in the library
        pseudogene_size = int(adata.var[~adata.var.targetType.eq(ctrl_label)].groupby('target').size().mean())

    if num_pseudogenes == 'auto':
        # approx number of target in the library
        num_pseudogenes = len(adata.var.loc[~adata.var.targetType.eq(ctrl_label),'target'].unique()) * pseudogene_size
    
    adata_ctrl = adata[:,adata.var.targetType.eq(ctrl_label)].copy()
    ctrl_elements = adata_ctrl.var.index.to_list()
    
    adata_pseudo_list = []
    pseudo_source_sgrna = []

    for pseudo_num in range(0, num_pseudogenes, pseudogene_size):
        pseudo_elements = np.random.choice(ctrl_elements, pseudogene_size, replace=False)
        pseudo_labels = [f'pseudo_{pseudo_num}_{i}' for i in range(1,pseudogene_size+1)]
        
        adata_pseudo = ad.AnnData(
            X = adata_ctrl.X[:,adata_ctrl.var.index.isin(pseudo_elements)],
            obs = adata_ctrl.obs   
        )
        adata_pseudo.var_names = pseudo_labels
        adata_pseudo_list.append(adata_pseudo)
        
        for element in pseudo_elements: 
            pseudo_source_sgrna.append(element)
        
    out = ad.concat(adata_pseudo_list, axis=1)
    out.var['target'] = out.var.index.str.split('_').str[:-1].str.join('_')
    out.var['targetType'] = ctrl_label
    out.var['source'] = pseudo_source_sgrna
    out.obs = adata_ctrl.obs.copy()
    
    return out
