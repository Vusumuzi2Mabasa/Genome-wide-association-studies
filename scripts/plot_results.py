#!/usr/bin/env python3
"""Reproduce PLINK figures and result tables from exported synthetic results."""
import argparse,csv,json,math,platform
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import beta,chi2

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,default=Path('results'));ap.add_argument('--data',type=Path,default=Path('data'));a=ap.parse_args()
    g=a.results/'genotype';out=a.results/'figures';out.mkdir(exist_ok=True)
    rows=read(g/'assoc.PHENO.glm.logistic.hybrid');attempted=[r for r in rows if r['TEST']=='ADD']
    valid=[r for r in attempted if r['ERRCODE']=='.' and math.isfinite(float(r['P'])) and 0<float(r['P'])<=1]
    expected=json.loads((a.results/'validation.json').read_text())['genotype']
    assert len(attempted)==expected['retained_variants'] and len(valid)==expected['valid_tests']
    pv=np.array([float(r['P']) for r in valid]);threshold=.05/len(attempted)
    ranked=[]
    for r in sorted(valid,key=lambda r:float(r['P'])):
        odds=float(r['OR']);se=float(r['LOG(OR)_SE'])
        ranked.append(dict(variant=r['ID'],chromosome=r['#CHROM'],position=r['POS'],effect_allele=r['A1'],other_allele=r['OMITTED'],effect_allele_frequency=r['A1_FREQ'],n=r['OBS_CT'],odds_ratio=odds,ci95_lower=math.exp(math.log(odds)-1.96*se),ci95_upper=math.exp(math.log(odds)+1.96*se),p=float(r['P']),bonferroni_significant=float(r['P'])<threshold,firth=r['FIRTH?']))
    write(g/'plink_ranked_results.tsv',ranked)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,ax=plt.subplots(figsize=(11,4.8),layout='constrained');offset=0;ticks=[];labels=[]
    for index,c in enumerate(sorted({int(r['#CHROM']) for r in valid})):
        rr=sorted([r for r in valid if int(r['#CHROM'])==c],key=lambda r:int(r['POS']))
        xx=np.array([int(r['POS']) for r in rr]);coords=xx-xx.min()+offset
        ax.scatter(coords,[-math.log10(float(r['P'])) for r in rr],s=22,color=['#146C94','#53A88D'][index%2])
        ticks.append((coords.min()+coords.max())/2);labels.append(str(c));offset=coords.max()+150000
    ax.axhline(-math.log10(threshold),color='#C76A15',ls='--',label=f'Bonferroni: 0.05/{len(attempted)}')
    ax.axhline(-math.log10(5e-8),color='#9B2948',ls=':',label='Genome-wide: 5 x 10^-8')
    ax.set(xticks=ticks,xticklabels=labels,xlabel='Synthetic chromosome (positions concatenated)',ylabel='-log10(P)',title=f'PLINK association results | {len(valid)} valid tests',ylim=(0,8))
    ax.legend(loc='upper right',fontsize=9);fig.savefig(out/'plink_manhattan.png');plt.close(fig)
    obs=-np.log10(np.sort(pv));n=len(pv);ranks=np.arange(1,n+1);exp=-np.log10((ranks-.5)/n)
    lo=-np.log10(beta.ppf(.975,ranks,n+1-ranks));hi=-np.log10(beta.ppf(.025,ranks,n+1-ranks))
    inflation=float(np.median(chi2.isf(pv,1))/chi2.ppf(.5,1))
    fig,ax=plt.subplots(figsize=(6,5.4),layout='constrained');ax.fill_between(exp,lo,hi,color='#D7E7EF',label='95% pointwise null envelope')
    ax.scatter(exp,obs,s=22,color='#146C94');limit=max(hi.max(),obs.max())+.1
    ax.plot([0,limit],[0,limit],color='#555555',ls='--');ax.set(xlim=(0,limit),ylim=(0,limit),xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title=f'PLINK QQ plot | lambda GC = {inflation:.3f}')
    ax.legend(fontsize=9,loc='upper left');fig.savefig(out/'plink_qq.png');plt.close(fig)
    pcs=read(g/'pcs.eigenvec');cov={r['IID']:r for r in read(a.data/'genotype/covariates.tsv')};ph={r['IID']:r for r in read(a.data/'genotype/phenotype.tsv')}
    assert len({r['IID'] for r in pcs})==len(pcs)==expected['retained_samples']
    joined=[dict(IID=r['IID'],POP=cov[r['IID']]['POP'],PHENO=ph[r['IID']]['PHENO'],**{k:float(r[k]) for k in ['PC1','PC2','PC3']}) for r in pcs]
    assert all(math.isfinite(r[k]) for r in joined for k in ['PC1','PC2','PC3']);write(g/'pca_metadata.tsv',joined)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    populations=sorted({r['POP'] for r in joined})
    for ax,y in zip(axes,['PC2','PC3']):
        for pop,color in zip(populations,['#146C94','#D97732']):
            for phenotype,marker,label in [('1','o','control'),('2','^','case')]:
                rr=[r for r in joined if r['POP']==pop and r['PHENO']==phenotype]
                ax.scatter([r['PC1'] for r in rr],[r[y] for r in rr],c=color,marker=marker,s=42,alpha=.85,label=f'{pop}: {label}')
        ax.set(xlabel='PC1 score',ylabel=y+' score');ax.axhline(0,c='#DDDDDD',lw=.7,zorder=0);ax.axvline(0,c='#DDDDDD',lw=.7,zorder=0)
    axes[1].legend(fontsize=8);fig.suptitle(f'PLINK PCA | {len(pcs)} samples, {expected["pca_markers"]} LD-pruned variants')
    fig.savefig(out/'plink_pca.png');plt.close(fig)
    metrics=dict(attempted_tests=len(attempted),plotted_tests=n,excluded_tests=len(attempted)-n,bonferroni_threshold=threshold,lambda_gc=inflation,minimum_p=float(pv.min()),top_variant=ranked[0]['variant'],pca_samples=len(pcs),pca_markers=expected['pca_markers'],python=platform.python_version(),matplotlib=matplotlib.__version__,numpy=np.__version__,scipy=scipy.__version__)
    (a.results/'plot_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
