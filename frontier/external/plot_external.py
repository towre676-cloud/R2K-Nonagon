#!/usr/bin/env python3
import json,numpy as np,matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'confirmation_results.json').read_text())['results']
methods=['loewner','vector_fitting','cubic_log_frequency','positive_real_relaxation']
names=['Loewner','Vector\nfitting','Cubic\n(log frequency)','Positive RC\nrelaxation']
colors=['#2e5d74','#b16b45','#6f8b64','#8b648d']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(11,7.6),sharey=True)
for ax,(scenario,budget) in zip(axs.flat,[(s,b) for s in ['interleaved','missing_band'] for b in [12,20]]):
    vals=[np.array([100*x['test_error'] for x in r if x['scenario']==scenario and x['budget']==budget and x['method']==m]) for m in methods]
    bp=ax.boxplot(vals,tick_labels=names,patch_artist=True,widths=.6,showfliers=True,medianprops={'color':'black','linewidth':1.4})
    for patch,c in zip(bp['boxes'],colors):patch.set_facecolor(c);patch.set_alpha(.65)
    ax.set_yscale('log');ax.set_ylim(.025,3);ax.set_title(('Interleaved holdout' if scenario=='interleaved' else 'Contiguous missing band')+f' · {budget} training points')
    ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True)
    for i,v in enumerate(vals,1):ax.text(i,2.35,f'{np.median(v):.3f}%',ha='center',fontsize=9)
axs[0,0].set_ylabel('Heldout complex relative RMS error (%)')
axs[1,0].set_ylabel('Heldout complex relative RMS error (%)')
fig.suptitle('48 protected measured battery cells: structure must earn its accuracy',fontweight='bold',fontsize=14,y=.99)
fig.text(.5,.015,'Same training data per method · validation selects candidates · each box contains 48 physical cells · medians labeled · lower is better',ha='center',fontsize=9,color='#444444')
fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig(ROOT/'external_comparison.png',dpi=180)
