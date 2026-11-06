from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flint import arb
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'summary_results.json').read_text());rows=d['cases'];ids=sorted({r['n0'] for r in rows});Ls=[.5,1,2,4,8,16]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
fig,(ax,bx)=plt.subplots(1,2,figsize=(13,5.6),gridspec_kw={'width_ratios':[1.7,1]})
colors=plt.get_cmap('tab10').colors
for j,n in enumerate(ids):
 rr=[r for r in rows if r['n0']==n];offset=(j-(len(ids)-1)/2)*.11
 ax.errorbar(np.arange(6)+offset,[r['residual_midpoint']*1000 for r in rr],yerr=[r['nominal_block_se']*1000 for r in rr],fmt='o',capsize=2,ms=4,color=colors[j],label=f"{n:,}"+(' held out' if rr[0]['kind']=='held_out' else ''))
ax.axhline(0,color='#637186',lw=1);ax.set_xticks(range(6),[str(x) for x in Ls]);ax.set_xlabel('Interval length L');ax.set_ylabel('Empirical minus model variance (×10⁻³)');ax.set_title('Fresh windows; bars are nominal ten-block SE');ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=3);ax.grid(axis='y',alpha=.15)
x=np.arange(len(ids));model=[];emp=[];sampling=[]
for n in ids:
 rr=[r for r in rows if r['n0']==n]
 model.append(max(float(arb(r['model_infinite_ball']).rad().upper()) for r in rr))
 emp.append(max(float(arb(r['empirical_ordinate_unfolding_error_upper']).upper()) for r in rr))
 sampling.append(max(r['nominal_block_se'] for r in rr))
bx.plot(x,model,'o-',label='Certified model radius',color='#0b7285');bx.plot(x,emp,'s-',label='Certified empirical error',color='#7950a6');bx.plot(x,sampling,'^-',label='Nominal sampling SE',color='#c67618');bx.axhline(1e-6,color='#9c2542',ls='--',lw=1,label='10⁻⁶ numerical target');bx.set_yscale('log');bx.set_xticks(x,[f'{n:,}' for n in ids],rotation=30,ha='right');bx.set_ylabel('Variance scale');bx.set_title('Numerical error versus sampling uncertainty');bx.grid(axis='y',alpha=.15);bx.legend(fontsize=8,loc='best')
fig.tight_layout();fig.subplots_adjust(bottom=.25);fig.savefig(HERE/'residuals.png',dpi=160);fig.savefig(HERE/'residuals.svg');print('Plot generated',len(ids),'windows')
