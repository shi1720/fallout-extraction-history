"""Reproduce a bounded historical reanalysis; no inferred environmental inventories."""
from __future__ import annotations
import csv, json, math, platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parents[1]

def read_soils():
    with (ROOT/'data/soil_1959.csv').open() as f:
        return list(csv.DictReader(f))

def proportional_interval(x,y,eps):
    """True values lie in reported_value*(1 +/- eps); h=r*f, 0<=r<=1.
    Returns the intersection; lower>upper denotes infeasibility.
    eps is a sensitivity assumption, NOT a confidence level.
    """
    if not 0<=eps<1: raise ValueError('eps must lie in [0,1)')
    x,y=np.asarray(x,float),np.asarray(y,float)
    if np.any(x<=0) or np.any(y<=0):raise ValueError('positive measurements required')
    return max(0.,float(np.max(y/x*(1-eps)/(1+eps)))), min(1.,float(np.min(y/x*(1+eps)/(1-eps))))

def additive_interval(x,y,eps):
    """h=f-d, d>=0, under the same componentwise relative tolerance."""
    if not 0<=eps<1: raise ValueError('eps must lie in [0,1)')
    x,y=np.asarray(x,float),np.asarray(y,float)
    return max(0.,float(np.max((1-eps)*x-(1+eps)*y))),float(np.min((1+eps)*x-(1-eps)*y))

def minimum_tolerance(x,y,model):
    """Closed forms, exact up to floating arithmetic; no statistical independence assumed."""
    x,y=np.asarray(x,float),np.asarray(y,float)
    if model=='proportional':
        r=y/x
        # Includes the physical r<=1 requirement.
        v=math.sqrt(float(np.max(r)/np.min(r)))
        return max((v-1)/(v+1),float(np.max((y-x)/(y+x))),0.)
    if model=='additive':
        d=x-y;s=x+y
        # Intersection [d_i-e*s_i,d_i+e*s_i] with [0,infinity).
        return max(0.,float(np.max((d[:,None]-d[None,:])/(s[:,None]+s[None,:]))),float(np.max(-d/s)))
    raise ValueError(model)

def fit_models(x,y):
    X=np.column_stack([np.ones(len(x)),x]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
    r=np.dot(x,y)/np.dot(x,x);d=float(np.mean(x-y))
    return {'affine':{'intercept':float(beta[0]),'slope':float(beta[1]),'x_intercept':float(-beta[0]/beta[1]),'sse':float(np.sum((y-X@beta)**2))},
            'proportional':{'fraction':float(r),'sse':float(np.sum((y-r*x)**2))},
            'additive':{'loss':d,'sse':float(np.sum((y-(x-d))**2))}}

def main():
    rows=read_soils();x=np.array([float(r['fusion_mCi_mi2']) for r in rows]);y=np.array([float(r['hcl_mCi_mi2']) for r in rows]);ratio=y/x
    models=fit_models(x,y);loo=[]
    for i in range(len(x)):
        keep=np.arange(len(x))!=i;f=fit_models(x[keep],y[keep]);a=f['affine'];p=f['proportional'];d=f['additive']
        loo.append({'omitted_site':i+1, 'affine_x_intercept':a['x_intercept'],
                    'affine_heldout_squared_error':float((y[i]-(a['intercept']+a['slope']*x[i]))**2),
                    'proportional_heldout_squared_error':float((y[i]-p['fraction']*x[i])**2),
                    'additive_heldout_squared_error':float((y[i]-(x[i]-d['loss']))**2),
                    'proportional_min_tolerance':minimum_tolerance(x[keep],y[keep],'proportional'),
                    'additive_min_tolerance':minimum_tolerance(x[keep],y[keep],'additive')})
    # The old SSE can be reconciled with integer printing without changing x.
    a=models['affine'];pred=a['intercept']+a['slope']*x
    factor=math.sqrt(359/a['sse']);yw=pred+factor*(y-pred)
    out={'n':len(rows),'mean_percent':float(ratio.mean()*100),'sample_sd_percent':float(ratio.std(ddof=1)*100),
         'population_sd_percent':float(ratio.std(ddof=0)*100),'min_percent':float(ratio.min()*100),'max_percent':float(ratio.max()*100),
         'ratio_of_sums_percent':float(y.sum()/x.sum()*100),'n_below_fusion':int(np.sum(y<x)),
         'models':models,'minimum_tolerance':{m:minimum_tolerance(x,y,m) for m in ['proportional','additive']},
         'interval_at_10_percent':{'proportional':proportional_interval(x,y,.1),'additive':additive_interval(x,y,.1)},
         'rounding_witness':{'target_sse':359,'max_hcl_adjustment':float(np.max(np.abs(yw-y))), 'sse_recomputed':fit_models(x,yw)['affine']['sse']},
         'site10_ratio_range_under_half_unit_rounding':[100*13.5/21.5,100*14.5/20.5],
         'leave_one_out':{'x_intercept_range':[min(r['affine_x_intercept'] for r in loo),max(r['affine_x_intercept'] for r in loo)],
                         'prediction_rmse':{m:math.sqrt(np.mean([r[m+'_heldout_squared_error'] for r in loo])) for m in models}},
         'software':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__}}
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results/summary.json').write_text(json.dumps(out,indent=2)+'\n')
    with (ROOT/'results/leave_one_out.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(loo[0]));w.writeheader();w.writerows(loo)
    with (ROOT/'results/rounding_witness.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['site_id','fusion_printed','hcl_printed','hcl_compatible_unrounded']);w.writerows(zip(range(1,17),x,y,yw))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    colors={'blue':'#22577A','rust':'#B44D30','green':'#387B65','grey':'#68717A'}
    fig,axs=plt.subplots(1,2,figsize=(10,4.1),layout='constrained')
    ax=axs[0];ax.scatter(x,y,c=colors['blue'],s=35,zorder=4)
    for i,(xx,yy) in enumerate(zip(x,y),1):ax.annotate(str(i),(xx,yy),xytext=(4,4),textcoords='offset points',fontsize=8)
    xx=np.array([0,115]);ax.plot(xx,xx,':',color=colors['grey'],label='Equality')
    ax.plot(xx,a['intercept']+a['slope']*xx,color=colors['rust'],label='Free affine fit')
    ax.plot(xx,models['proportional']['fraction']*xx,'--',color=colors['green'],label='Proportional fit')
    ax.set(xlabel='Fusion (mCi per square mile)',ylabel='HCl extraction (mCi per square mile)',xlim=(0,117),ylim=(0,117),title='A  Published method comparison')
    ax.legend(fontsize=8,frameon=False,loc='upper left')
    ax=axs[1];ax.scatter(range(1,17),ratio*100,color=colors['blue'],zorder=3)
    ax.axhspan(78,101,color=colors['grey'],alpha=.13,label='Printed summary range: 78–101%')
    ax.axhline(100,color=colors['grey'],ls=':');ax.scatter([10],[ratio[9]*100],color=colors['rust'],s=55,zorder=5,label='Site 10: printed 14 / 21')
    ax.set(xlabel='Published site identifier',ylabel='HCl / fusion (%)',xticks=range(1,17),ylim=(60,108),title='B  Preserve the table / summary discrepancy')
    ax.legend(fontsize=8,frameon=False,loc='lower left')
    for ext in ['pdf','png']:fig.savefig(ROOT/f'figures/soil_comparison.{ext}',dpi=220)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    ax=axs[0]
    thresholds=[100*out['minimum_tolerance'][m] for m in ['proportional','additive']]
    ax.barh([1,0],thresholds,color=[colors['green'],colors['rust']],height=.42)
    for iy,t in zip([1,0],thresholds):ax.text(t+.25,iy,f'{t:.2f}%',va='center',fontsize=10)
    ax.axvline(10,color=colors['grey'],lw=1,ls=':',label='Historical 10% discussion scale')
    ax.set(xlabel='Minimum relative tolerance per measurement (%)',yticks=[0,1],yticklabels=['Fixed additive loss','Common extraction fraction'],xlim=(0,13),ylim=(-.6,1.6),title='A  Compatibility thresholds')
    ax.legend(frameon=False,fontsize=8,loc='upper left')
    ax=axs[1];ax.scatter(range(1,17),[z['affine_x_intercept'] for z in loo],color=colors['blue'])
    ax.axhline(a['x_intercept'],color=colors['rust'],ls='--',label=f'All 16: {a["x_intercept"]:.2f}')
    ax.axhline(0,color=colors['grey'],lw=.7);ax.set(xlabel='Omitted site identifier',ylabel='Affine x-intercept (mCi per square mile)',xticks=range(1,17),title='B  Leave-one-site-out sensitivity')
    ax.legend(frameon=False,fontsize=8)
    for ext in ['pdf','png']:fig.savefig(ROOT/f'figures/model_sensitivity.{ext}',dpi=220)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.2));ax.set_axis_off();ax.set(xlim=(0,10),ylim=(0,4.2))
    boxes=[(.1,2.3,2.6,1.45,'COLLECTION RESIN','Target: deposited fallout\nUnwanted: resin-bound activity\nAcid leaves unwanted activity behind'),(.1,.25,2.6,1.45,'SAMPLED SOIL','Target: sampled-soil inventory\nConcern: unextracted fallout\nAcid may leave wanted activity behind'),(3.4,2.3,2.8,1.45,'CORRECTION REPORTED','Unused-column checks\n0.24 mCi / mi² per column\nAcid procedure from March 1961'),(3.4,.25,2.8,1.45,'CORRECTION DEFERRED','Contaminated fusion setting\n16-soil comparison; no replicates\nNew blind tests planned at NBS'),(7,1.2,2.7,1.7,'HISTORICAL CLAIM','Method, material, and target\nmust travel with each number.\nMore complete dissolution is\nnot always more valid measurement.')]
    for j,(bx,by,bw,bh,head,body) in enumerate(boxes):
        ax.add_patch(FancyBboxPatch((bx,by),bw,bh,boxstyle='round,pad=0.08',fc='#F3F6F8',ec=colors['blue'] if j!=3 else colors['rust'],lw=1.1))
        ax.text(bx+.12,by+bh-.20,head,fontsize=9,fontweight='bold',va='top',color=colors['blue'])
        ax.text(bx+.12,by+bh-.58,body,fontsize=9,va='top',linespacing=1.6)
    for p,q in [((2.8,3),(3.25,3)),((2.8,.95),(3.25,.95)),((6.35,3),(6.85,2.45)),((6.35,.95),(6.85,1.75))]:ax.annotate('',xy=q,xytext=p,arrowprops={'arrowstyle':'->','color':colors['grey'],'lw':1.3})
    fig.tight_layout()
    for ext in ['pdf','png']:fig.savefig(ROOT/f'figures/two_extractions.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
