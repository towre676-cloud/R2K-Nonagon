from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=Path(__file__).resolve().parents[1];OUT=Path(__file__).parent/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#281e19','text.color':'#281e19','xtick.color':'#281e19','ytick.color':'#281e19','svg.fonttype':'none','axes.edgecolor':'#b9aa92','axes.facecolor':'none','figure.facecolor':'none','grid.color':'#b9aa92'})
ox='#7c2435';teal='#35606a';grey='#8b7966'
def save(fig,name):
 fig.tight_layout(pad=1.3);fig.savefig(OUT/(name+'.svg'),transparent=True);plt.close(fig)
 p=OUT/(name+'.svg');s=p.read_text().replace('#281e19','var(--ink)').replace('#b9aa92','var(--rule)').replace('#7c2435','var(--ox)').replace('#35606a','var(--teal)').replace('#8b7966','var(--soft)');p.write_text(s)
D=json.loads((Path(__file__).resolve().parent/'dynamic_results.json').read_text());rows=D['counts'];n=np.array([x['n'] for x in rows]);fig,ax=plt.subplots(figsize=(7.2,3.3));ax.plot(n,[x['old_real_readings'] for x in rows],'o-',color=grey,label='Three complete real-frequency inverses');ax.plot(n,[x['dimension'] for x in rows],'s-',color=ox,label='Elimination design: n(n+2) real scalars');ax.set(xlabel='Number of labeled nodes n',ylabel='Real scalar readings');ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.3);save(fig,'measurement-count')
F=json.loads((Path(__file__).resolve().parent/'forests_results.json').read_text());fig,axs=plt.subplots(1,2,figsize=(8.6,3.3));a=F['hard_families']['weak_link'];x=np.array([z['w'] for z in a]);y=np.array([z['resolution_scale'] for z in a]);axs[0].loglog(x,y,'o-',color=ox);axs[0].set(xlabel='Bridge weight w',ylabel='Edge resolution factor Γ⁻²');axs[0].grid(alpha=.25);a=F['hard_families']['distance'];axs[1].semilogy([z['depth'] for z in a],[z['resolution_scale'] for z in a],'s-',color=teal);axs[1].set(xlabel='Stem depth d',ylabel='Edge resolution factor Γ⁻²');axs[1].grid(alpha=.25);save(fig,'edge-resolution')
I=json.loads((Path(__file__).resolve().parent/'design_results.json').read_text());fig,ax=plt.subplots(figsize=(7.2,3.2));x=np.arange(3);gain=[z['greedy_logdet_gain'] for z in I['cases']];mean=[z['uniform_subset_mean_logdet_gain'] for z in I['cases']];ax.bar(x-.17,mean,.34,color=grey,label='Uniform subset: mean logdet');ax.bar(x+.17,gain,.34,color=ox,label='Greedy = exhaustive optimum');ax.set_xticks(x,['Grounded path','2 × 2 grid','Asymmetric diamond']);ax.set(ylabel='Information logdet gain');ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.25);save(fig,'information-design')
O=json.loads((Path(__file__).resolve().parent/'origin_cauchy_results.json').read_text());a=O['infinite_product_rows'];import re
lo=lambda s: float(re.search(r'[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?',s).group())
fig,ax=plt.subplots(figsize=(7.2,3.1));T=[r['T'] for r in a];yl=[lo(r['uniform_R2_over_x2_lower']) for r in a];yu=[lo(r['uniform_R2_over_x2_upper']) for r in a];m=(np.array(yl)+np.array(yu))/2;ax.errorbar(T,m,yerr=[m-np.array(yl),np.array(yu)-m],fmt='o',color=ox,capsize=4);ax.axhline(3,color=teal,ls='--',lw=1);ax.text(T[-1],3.015,'c = 3',ha='right',color=teal);ax.set(xscale='log',xlabel='Height T',ylabel='Uniform bounds on R₂(x) / x²');ax.set_ylim(2.4,3.08);ax.grid(axis='y',alpha=.25);save(fig,'origin-enclosure')
print('four scientific SVG figures generated')
