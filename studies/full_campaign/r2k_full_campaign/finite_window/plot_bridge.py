from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from flint import arb
HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'finite_window_summary.json').read_text());b=json.loads((HERE/'sine_bias_results.json').read_text())
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
fig,(ax,bx)=plt.subplots(1,2,figsize=(12.5,4.8))
for kind,color,label in [('point_selected','#167d9a','Point-selected'),('fixed_span','#b45725','Fixed coordinate span')]:
 w=next(w for w in s['windows'] if w['geometry']==kind and w['n0']==20000000);tiers=w['dependence_tiers']
 ax.plot([t['batches'] for t in tiers],[t['observed_max_absolute_z'] for t in tiers],'o-',color=color,label=label+' observed max')
 ax.plot([t['batches'] for t in tiers],[t['joint_reference_95'] for t in tiers],'--',color=color,alpha=.65,label=label+' reference 95%')
ax.set_xscale('log',base=2);ax.set_xticks([5,10,20,40,80],['5','10','20','40','80']);ax.set_xlabel('Number of anchor batches');ax.set_ylabel('Maximum absolute standardized residual');ax.set_title('Index 20,000,000: sensitivity to batch choice');ax.grid(alpha=.15);ax.legend(fontsize=8)
rr=next(w for w in s['windows'] if w['geometry']=='fixed_span' and w['n0']==20000000)['rows'];xs=[r['L'] for r in rr]
bx.plot(xs,[float(arb(r['residual_ball']).abs_upper()) for r in rr],'o-',color='#b45725',label='Observed absolute residual')
bx.plot(xs,[float(arb(r['combined_deterministic_radius']).upper()) for r in rr],'s-',color='#167d9a',label='Certified numerical radius')
bx.plot(xs,[float(arb(r['expected_downward_bias_ball']).mid()) for r in b['production_fixed_span_benchmarks']],'^-',color='#7950a6',label='Conditional sine-model centering bias')
bx.set_xscale('log',base=2);bx.set_xticks(xs,[str(x) for x in xs]);bx.set_yscale('log');bx.set_xlabel('Interval length L');bx.set_ylabel('Variance scale');bx.set_title('Finite-window mean bias is a separate quantity');bx.grid(alpha=.15);bx.legend(fontsize=8)
fig.text(.5,.015,'Multiplier references are diagnostics; independence and sign symmetry are not established for the zeta data.',ha='center',fontsize=9,color='#573f40')
fig.tight_layout(rect=[0,.045,1,1]);fig.savefig(HERE/'finite_window_diagnostics.png',dpi=160);fig.savefig(HERE/'finite_window_diagnostics.svg')
