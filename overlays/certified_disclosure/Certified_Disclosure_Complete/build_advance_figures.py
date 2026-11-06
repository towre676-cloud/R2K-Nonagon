#!/usr/bin/env python3
"""Regenerate the numerical advance figures from retained certificate ledgers.

Plots display binary64 projections of exact rational/dyadic endpoints. The
certificate ledgers remain authoritative. No values are fitted or invented.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import argparse
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MultipleLocator, FuncFormatter

BASE=Path(__file__).resolve().parent
POS=BASE/'possq-0.4.0/results'
NON=BASE/'r2k_nonagon_v14.1/part_XIV'
OUT=BASE/'figures'
NAVY='#18354a';TEAL='#137e83';COPPER='#b4653f';GRAY='#7d8e98';PALE='#e6edef'

plt.rcParams.update({
 'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
 'axes.labelsize':10,'axes.labelcolor':NAVY,'text.color':NAVY,
 'xtick.color':NAVY,'ytick.color':NAVY,'axes.edgecolor':'#cad4d8',
 'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold',
 'axes.titlelocation':'left','axes.grid':False,'legend.frameon':False,
 'savefig.facecolor':'white','figure.facecolor':'white','svg.fonttype':'none',
 'pdf.fonttype':42,'lines.linewidth':1.8})

def load(path):return json.loads(path.read_text())
def dyadic(d):return float(Fraction(int(d['mantissa']))*Fraction(2)**int(d['exponent2']))
def clean(ax):
 ax.yaxis.grid(True,color=PALE,linewidth=.7);ax.set_axisbelow(True)
 ax.tick_params(length=3,color='#bac6cc',labelsize=9)
def pct(x,pos):
 return f'{100*x:,.0f}%' if abs(100*x)>=1000 else f'{100*x:.3g}%'
def save(fig,name):
 OUT.mkdir(exist_ok=True)
 for suffix in ['svg','png','pdf']:
  fig.savefig(OUT/f'{name}.{suffix}',bbox_inches='tight',dpi=220 if suffix=='png' else None)
 plt.close(fig)

def spectral():
 d=load(POS/'resolvent_campaign.json')
 count=len(d['campaigns']);cols=2 if count>=4 else count;rows=math.ceil(count/cols)
 fig,axs=plt.subplots(rows,cols,figsize=(12.8,4.15*rows),sharey=True,squeeze=False)
 axs=axs.ravel()
 truth=float(Fraction(d['truth_mass']))
 for i,(ax,row) in enumerate(zip(axs,d['campaigns'])):
  for key,color,style,label in [('baseline',GRAY,'--','Wide-band readings'),('adaptive',TEAL,'-','Nominal-model-selected readings')]:
   points=sorted(row[key],key=lambda z:z['readings'])
   x=np.array([p['readings'] for p in points]);lower=np.array([float(Fraction(p['lower']['value'])) for p in points]);upper=np.array([float(Fraction(p['upper']['value'])) for p in points])
   ax.plot(x,upper,style,color=color,marker='o' if key=='adaptive' else 's',ms=4,label=label)
   ax.plot(x,lower,style,color=color,marker='o' if key=='adaptive' else 's',ms=4)
   if key=='adaptive':ax.fill_between(x,lower,upper,color=color,alpha=.075)
  ax.axhline(truth,color=COPPER,lw=1.4,label='Reference mass = 1/16')
  exponent=round(math.log10(float(Fraction(row['noise_relative_h0']))))
  ax.set_title(f'{chr(65+i)}  Reading budget 10$^{{{exponent}}}$ $h(0)$',pad=13)
  max_readings=max(p['readings'] for key in ['baseline','adaptive'] for p in row[key])
  ax.set_xlabel('Scalar response readings');ax.set_xticks(np.arange(2,max_readings+1,2));ax.set_xlim(1.7,max_readings+.3);ax.set_ylim(-.005,.435)
  ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));clean(ax)
 for i,ax in enumerate(axs):
  if i>=count:ax.set_visible(False)
  elif i%cols==0:ax.set_ylabel('Mass in the slow-decay band')
 handles,labels=axs[0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.035),ncol=3,fontsize=9)
 fig.subplots_adjust(wspace=.14,hspace=.38,bottom=.13 if rows>1 else .22,top=.92 if rows>1 else .86)
 save(fig,'spectral_disclosure')
 return 'Exact rational upper and lower spectral-mass certificates on the 32-node, 80-edge grounded hypercube. The target is [1/64, 1/8], whose reference port-visible spectral mass is 1/16. Teal envelopes use nominal-model-selected shifts; grey envelopes use wide-band shifts at the same reading count. Selection uses the nominal model to score proposed future readings, while the certificates use only the selected observations and declared budgets. Each panel uses a declared absolute reading-error budget relative to h(0). Every endpoint is independently Sturm-verified; acquisition is a heuristic, not a proved optimum.'

def spectral_witnesses():
 d=load(POS/'resolvent_campaign.json')
 row=min(d['campaigns'],key=lambda r:float(Fraction(r['noise_relative_h0'])))
 r=max(row['adaptive'],key=lambda r:r['readings'])
 if 'upper_witness' not in r:return None
 lo,hi=map(Fraction,d['I']);a,b=map(Fraction,d['target'][0])
 shifts=list(map(Fraction,r['shifts']))
 xx=np.geomspace(float(lo),float(hi),1000)
 def evaluate(coeffs):
  cs=list(map(Fraction,coeffs))
  return np.array([float(cs[0]+sum((c/(Fraction(float(x))+s) for c,s in zip(cs[1:],shifts)),Fraction(0))) for x in xx])
 fig,axs=plt.subplots(1,2,figsize=(12.8,4.55))
 ax=axs[0]
 ax.axvspan(float(a),float(b),color=TEAL,alpha=.055)
 ax.plot(xx,evaluate(r['upper']['coeffs']),color=TEAL,label='Upper sandwich q⁺')
 ax.plot(xx,evaluate(r['lower']['coeffs']),color=COPPER,label='Lower sandwich q⁻')
 ax.plot([float(lo),float(b),float(b),float(hi)],[1,1,0,0],color=NAVY,ls='--',lw=1,label='Target indicator')
 ax.set_xscale('log');ax.set_yscale('symlog',linthresh=1,linscale=.8)
 ax.set_xlim(float(lo),float(hi));ax.set_xlabel('Decay rate λ');ax.set_ylabel('Rational dual test function')
 ax.set_title('A  The verified rational sandwich',pad=13);ax.legend(loc='lower right',fontsize=8.5)
 ax.text(.035,.93,f'{r["readings"]} raw response readings',transform=ax.transAxes,fontsize=9,color=GRAY)
 ax=axs[1]
 lambdas=list(map(Fraction,d['spectrum']));weights=list(map(Fraction,d['spectral_weights']))
 wx,ww=zip(*[(Fraction(x),Fraction(w)) for x,w in r['upper_witness']['atoms']])
 marker,stems,base=ax.stem(list(map(float,lambdas)),list(map(float,weights)),linefmt='-',markerfmt='o',basefmt=' ')
 plt.setp(marker,color=GRAY,markersize=3);plt.setp(stems,color=GRAY,linewidth=1.1,alpha=.7)
 marker,stems,base=ax.stem(list(map(float,wx)),list(map(float,ww)),linefmt='-',markerfmt='D',basefmt=' ')
 plt.setp(marker,color=TEAL,markersize=4);plt.setp(stems,color=TEAL,linewidth=1.4)
 ax.plot([],[],color=GRAY,marker='o',ms=4,label='Reference spectral measure')
 ax.plot([],[],color=TEAL,marker='D',ms=4,label='Exactly compatible atomic witness')
 ax.axvspan(float(a),float(b),color=TEAL,alpha=.055)
 ax.set_xscale('log');ax.set_xlim(float(lo)*.93,float(hi)*1.08)
 ax.set_xlabel('Decay rate λ');ax.set_ylabel('Atomic spectral mass');ax.yaxis.set_major_formatter(PercentFormatter(1))
 ax.set_title('B  Distinct measures fit the same reading boxes',pad=13);ax.legend(loc='upper left',fontsize=8.5)
 ref=float(Fraction(d['truth_mass']));wm=float(Fraction(r['upper_witness']['value']))
 ax.text(.035,.67,f'Slow-band reference mass  {ref:.2%}\nCompatible witness mass    {wm:.2%}',transform=ax.transAxes,fontsize=9,color=NAVY,linespacing=1.8)
 for ax in axs:clean(ax)
 fig.subplots_adjust(wspace=.25,bottom=.2,top=.86)
 save(fig,'spectral_witnesses')
 return f'The rational certificate and its witness at {r["readings"]} nominally selected response shifts and declared error budget {r["noise_relative_h0"]} h(0). Left shows the exact rational dual functions evaluated at plotting points, with a symmetric logarithmic ordinate outside [-1,1] to retain the complete oscillation range. Sturm verification establishes the sandwich on the full continuous support, not merely at plotted points. Right contrasts the reference spectrum with a distinct exactly compatible rational atomic witness. The witness puts {wm:.8%} mass in the slow band, versus reference mass {ref:.2%}, demonstrating remaining distributional ambiguity under the same reading boxes.'

def frontier():
 cycle,star=[load(POS/f'frontier_{g}.json') for g in ['cycle','star']]
 fig,axs=plt.subplots(1,2,figsize=(12.8,4.65))
 for graph,d,color in [('Cycle',cycle,TEAL),('Star',star,NAVY)]:
  pairs=sorted((int(k),v) for k,v in d['copies'].items())
  x=np.array([k for k,v in pairs]);y=np.array([math.log10(dyadic(v['upper'])) for k,v in pairs])
  axs[0].plot(x,y,color=color,label=f'{graph}: certified law')
  v=d['copies']['70'];axs[0].scatter([70],[math.log10(dyadic(v['upper']))],s=25,color=color,zorder=4)
 old=load(POS/'cyclecraft_ledger_cycle.json')['cases'][0]['reported_chernoff']
 x=np.arange(1,max(map(int,cycle['copies']))+1)
 axs[0].plot(x,x*np.log10(old),color=GRAY,linestyle='--',lw=1.4,label='Cycle: earlier single-tilt upper bound')
 axs[0].axhline(-100,color=COPPER,lw=.9,linestyle=':')
 axs[0].text(3,-98.2,'Risk specification 10⁻¹⁰⁰',color=COPPER,fontsize=9)
 axs[0].annotate('70 independent copies',xy=(70,math.log10(dyadic(cycle['copies']['70']['upper']))),
  xytext=(34,-113),fontsize=9,color=TEAL,arrowprops={'arrowstyle':'-','color':TEAL,'lw':.8})
 axs[0].set_title('A  Independent ten-period experiments',pad=13)
 axs[0].set_xlabel('Number of independent copies, K');axs[0].set_ylabel('log₁₀ wrong-direction probability')
 axs[0].set_xlim(0,81);axs[0].set_ylim(-119,1);axs[0].set_xticks([0,20,40,60,80]);axs[0].legend(loc='upper right',fontsize=8.5)
 horizon_max=0
 for graph,d,color in [('Cycle',cycle,TEAL),('Star',star,NAVY)]:
  path=POS/f'frontier_horizons_{graph.lower()}.json'
  if not path.exists():continue
  h=load(path);pairs=sorted((int(k),v) for k,v in h['horizons'].items())
  x=np.array([k for k,v in pairs]);y=np.array([math.log10(dyadic(v['probability']['upper'])) for k,v in pairs])
  axs[1].plot(x,y,color=color,marker='o',ms=4,label=f'{graph}: correlated trajectory')
  cw=d['optimized_all_horizon_cw'];xx=np.linspace(10,max(x),300)
  yy=math.log10(dyadic(cw['C_upper']))+xx*math.log10(dyadic(cw['rho_upper']))
  axs[1].plot(xx,yy,color=color,ls='--',lw=1.1,alpha=.75,label=f'{graph}: all-horizon upper bound')
  horizon_max=max(horizon_max,max(x))
 axs[1].axvspan(10,80,color=PALE,alpha=.5,zorder=-2)
 axs[1].text(87,-7,'Earlier horizon: 80 periods',fontsize=8.5,color=GRAY)
 axs[1].set_title('B  One uninterrupted trajectory',pad=13)
 axs[1].set_xlabel('Consecutive protocol periods, n');axs[1].set_ylabel('log₁₀ wrong-direction probability')
 axs[1].set_xlim(0,max(660,horizon_max*1.025));axs[1].set_ylim(-105,1);axs[1].legend(loc='lower left',fontsize=8.2)
 for ax in axs:clean(ax)
 fig.subplots_adjust(wspace=.25,bottom=.18,top=.86)
 save(fig,'counting_frontier')
 return 'Two distinct certified counting experiments. Left: sums over independent copies of a ten-period experiment, compared with the inherited single-tilt Chernoff route. Right: one uninterrupted trajectory with inter-period correlations retained by the joint count-state kernel. Solid curves show rigorous finite-event enclosures whose widths are below plotting resolution; dashed right-panel curves are certified Collatz–Wielandt upper bounds. Both graphs first reach the 10^-100 independent-copy specification at K=70, with every smaller copy count excluded by a certified lower bound.'

def access():
 old=load(NON/'portloop14_case118.json')
 new=load(NON/'rational_access/rational_case118_pocket6.json')
 hw=load(NON/'rational_access/rational_case118_hw20.json')
 fig,axs=plt.subplots(1,2,figsize=(12.8,4.45),gridspec_kw={'width_ratios':[1.1,1]})
 hist=old['history'];x=[r['ports'] for r in hist];y=[r['rank'] for r in hist]
 axs[0].plot(x,y,color=GRAY,marker='o',label='Inherited placement loop')
 axs[0].scatter([len(new['ports'])],[new['rational_rank']],marker='D',s=55,color=TEAL,zorder=4,label='Pocket-directed placement')
 axs[0].axhline(new['unknowns'],color=COPPER,linestyle=':',lw=1.2,label='All 179 parameter directions')
 axs[0].annotate('Bus 55 closes the pocket',xy=(len(new['ports']),new['rational_rank']),xytext=(5.12,179.35),fontsize=9,color=TEAL,
  arrowprops={'arrowstyle':'-','lw':.8,'color':TEAL})
 axs[0].set_xlim(4.8,7.25);axs[0].set_ylim(177.7,179.7);axs[0].set_xticks([5,6,7]);axs[0].set_yticks([178,179])
 axs[0].set_xlabel('Measuring ports');axs[0].set_ylabel('Certified rank lower bound')
 axs[0].set_title('A  Placement completes the rank certificate',pad=13);axs[0].legend(loc='lower right',fontsize=8.5)
 values=[hist[-1]['readings'],new['readings']]
 colors=[GRAY,TEAL];labels=['Inherited seven-port design','Pocket-directed six-port design']
 bars=axs[1].barh([1,0],values,height=.45,color=colors)
 for bar,v in zip(bars,values):axs[1].text(v+12,bar.get_y()+bar.get_height()/2,f'{v:,} readings',va='center',fontsize=10,color=NAVY)
 axs[1].set_yticks([1,0],labels,fontsize=9);axs[1].set_xlim(0,max(values)*1.29);axs[1].set_ylim(-.7,1.7)
 axs[1].set_xlabel('Scalar readings at full rank');axs[1].set_title('B  Same complete local identification',pad=13)
 saving=1-values[1]/values[0]
 axs[1].text(.03,.92,f'{saving:.1%} fewer readings · one fewer port',transform=axs[1].transAxes,color=TEAL,fontsize=11,fontweight='bold')
 for ax in axs:clean(ax)
 axs[1].yaxis.grid(False);axs[1].xaxis.grid(True,color=PALE,linewidth=.7)
 fig.subplots_adjust(wspace=.62,bottom=.2,top=.85)
 save(fig,'access_design')
 return f'IEEE 118 access design. The inherited loop used {len(old["port_buses"])} ports and {values[0]} readings to certify full rank 179. Directing the sixth port to bus 55 certifies rank 179 with {values[1]} readings, a {saving:.2%} acquisition reduction. The inherited intermediate ranks are modular rank lower bounds; the new endpoint is an exact rational full-rank certificate. Full rank certifies local identification under the declared intervention model; this is an improved explicit construction, not a proof of minimum port count. Separately, the restricted shunt/generator scheme has exact rational rank {hw["rational_rank"]}/{hw["unknowns"]}, with {hw["first_order_identifiable_branches"]}/{hw["branches"]} branch coordinates certified identifiable to first order.'

def strength_scan():
 path=NON/'two_strength_cert_results.json'
 if not path.exists():return None
 d=load(path)
 cases=d.get('cases',[])
 if isinstance(cases,dict):cases=[dict(case=k,**v) for k,v in cases.items()]
 if not cases:return None
 fig,axs=plt.subplots(1,len(cases),figsize=(12.8,4.45),squeeze=False)
 for i,(ax,case) in enumerate(zip(axs[0],cases)):
  rows=case['rows'];name=case.get('case',case.get('name',f'case {i+1}'))
  for sigma,color in [('1/100000000',TEAL),('1/1000000',COPPER)]:
   rs=sorted([r for r in rows if r['sigma']==sigma and r.get('median_relative_radius') is not None],key=lambda r:float(Fraction(r['strength_scale'])))
   if not rs:continue
   x=[float(Fraction(r['strength_scale'])) for r in rs];y=[r['median_relative_radius'] for r in rs]
   exponent=round(math.log10(float(Fraction(sigma))))
   ax.loglog(x,y,color=color,marker='o',ms=4,label=f'Calibrated inputs, σ = 10$^{{{exponent}}}$')
   rr=[r for r in rs if r.get('raw_single_median_relative_radius') is not None]
   if rr:ax.loglog([float(Fraction(r['strength_scale'])) for r in rr],[r['raw_single_median_relative_radius'] for r in rr],color=color,marker='s',ls='--',ms=4,label=f'Raw single-node calibration, 10$^{{{exponent}}}$')
  ax.set_title(f'{chr(65+i)}  {name.replace("case","IEEE ")}',pad=13);ax.set_xlabel('Pair strength × median Gⱼⱼ');ax.set_ylabel('Median certified relative radius')
  ax.yaxis.set_major_formatter(FuncFormatter(pct));clean(ax)
 handles,labels=axs[0,0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.025),ncol=2,fontsize=8.5)
 fig.subplots_adjust(wspace=.3,bottom=.28,top=.85)
 save(fig,'two_strength_strengths')
 return 'Deterministic interval certification from two pair strengths. Solid curves condition on calibrated diagonal and port-column inputs; dashed curves propagate the raw single-node calibration uncertainties as well. Relative radii are scored against the reference network after inversion, while root selection itself uses only the observations. Increasing intervention strength can improve signal but can also reduce separation between the two response differences.'

def roots():
 path=NON/'two_strength_design_results.json'
 if not path.exists():return None
 d=load(path);rows=d['rows']
 if len(rows)<6:return None
 fig,axs=plt.subplots(1,2,figsize=(12.8,4.5),sharey=True)
 for i,(ax,kind,title) in enumerate(zip(axs,['calibrated','raw-single'],['Calibrated diagonal and port-column inputs','Calibration from raw single-node readings'])):
  rs=sorted([r for r in rows if r['calibration']==kind],key=lambda r:float(Fraction(r['sigma'])))
  for n,color,style,label in [('two',GRAY,'--','Two pair strengths'),('four',TEAL,'-','Four-strength intersection')]:
   xx=[float(Fraction(r['sigma'])) for r in rs]
   yy=[r['median_radius_'+n] for r in rs]
   ax.loglog(xx,yy,color=color,ls=style,marker='o' if n=='four' else 's',ms=5,label=label)
   # Thin matching curves retain the 90th-percentile branch radius without
   # presenting them as confidence or probability intervals.
   p90=[r['p90_radius_'+n] for r in rs]
   ax.loglog(xx,p90,color=color,ls=style,alpha=.4,lw=1)
  ax.set_title(f'{chr(65+i)}  {title}',pad=13)
  ax.set_xlabel('Declared reading-error budget / G₀₀')
  ax.set_xticks([1e-8,1e-6,1e-4]);ax.set_xticklabels(['10⁻⁸','10⁻⁶','10⁻⁴'])
  ax.set_xlim(4e-9,2.5e-4);ax.legend(loc='upper left',fontsize=9)
  ax.yaxis.set_major_formatter(FuncFormatter(pct));clean(ax)
 axs[0].set_ylabel('Certified relative interval radius')
 fig.text(.5,.025,'Heavy curves: median radius · thin curves: 90th-percentile radius',ha='center',fontsize=8.5,color=GRAY)
 fig.subplots_adjust(wspace=.18,bottom=.21,top=.84)
 save(fig,'two_strength_intervals')
 return f'IEEE 118 deterministic interval certification evaluated on all {d["pairs"]} lean pairs. Four declared pair strengths ({", ".join(d["scales"])}, normalized by median Gⱼⱼ) produce six pairwise elimination enclosures whose intersection can never widen a valid baseline enclosure. Heavy curves show the median scored relative radius among accepted pairs; thin curves show the 90th percentile. Left conditions on calibrated diagonal and port-column inputs. Right propagates uncertainty from raw single-node calibration. Reference values enter only the final inclusion and radius scores, not root selection. The four-strength protocol uses four pair readings, versus two for the baseline.'

def raw_cancellation():
 path=NON/'two_strength_raw_design_results.json'
 if not path.exists():return None
 d=load(path);rows=d['rows']
 if len(rows)<9:return None
 scales=sorted(set(r['strength_scale'] for r in rows),key=lambda x:float(Fraction(x)))
 fig,axs=plt.subplots(1,len(scales),figsize=(12.8,4.4),sharey=True)
 for i,(ax,scale) in enumerate(zip(axs,scales)):
  rs=sorted([r for r in rows if r['strength_scale']==scale],key=lambda r:float(Fraction(r['sigma'])))
  for key,color,style,label in [('before',GRAY,'--','Diagonal substitution'),('after',TEAL,'-','Raw algebraic cancellation')]:
   xx=[float(Fraction(r['sigma'])) for r in rs]
   yy=[r['median_radius_'+key] for r in rs]
   p90=[r['p90_radius_'+key] for r in rs]
   ax.loglog(xx,yy,color=color,ls=style,marker='o' if key=='after' else 's',ms=4,label=label)
   ax.loglog(xx,p90,color=color,ls=style,alpha=.35,lw=1)
  ax.set_title(f'{chr(65+i)}  Pair-strength scale {scale}',pad=13)
  ax.set_xlabel('Reading-error budget / G₀₀')
  ax.set_xticks([1e-8,1e-6,1e-4]);ax.set_xticklabels(['10⁻⁸','10⁻⁶','10⁻⁴']);ax.set_xlim(4e-9,2.5e-4)
  ax.yaxis.set_major_formatter(FuncFormatter(pct));clean(ax)
 axs[0].set_ylabel('Certified relative interval radius')
 handles,labels=axs[0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.025),ncol=2,fontsize=9)
 fig.text(.5,.04,'Heavy curves: median · thin curves: 90th percentile · identical raw reading budget',ha='center',fontsize=8.5,color=GRAY)
 fig.subplots_adjust(wspace=.15,bottom=.24,top=.85)
 save(fig,'raw_cancellation')
 return 'Exact interval inversion with identical raw data. Cancelling the inferred diagonal variables algebraically before interval evaluation reduces dependency inflation. Each panel fixes a pair-strength scale; each point is scored over accepted IEEE 118 lean pairs at the declared deterministic raw-reading budget. Heavy curves show median relative radius and thin curves the 90th percentile. The comparison uses the same four single-node and two pair differences, with no additional observations; truth is used only for inclusion and radius scoring.'

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--roots-only',action='store_true');a=parser.parse_args()
 captions={}
 if not a.roots_only:
  for name,fun in [('spectral_disclosure',spectral),('spectral_witnesses',spectral_witnesses),('counting_frontier',frontier),('access_design',access)]:
   captions[name]=fun()
 rootcap=roots()
 if rootcap:captions['two_strength_intervals']=rootcap
 scancap=strength_scan()
 if scancap:captions['two_strength_strengths']=scancap
 rawcap=raw_cancellation()
 if rawcap:captions['raw_cancellation']=rawcap
 if a.roots_only and (OUT/'captions.json').exists():captions={**load(OUT/'captions.json'),**captions}
 (OUT/'captions.json').write_text(json.dumps(captions,indent=2))
 (OUT/'CAPTIONS.md').write_text('\n\n'.join(f'**{k}.** {v}' for k,v in captions.items())+'\n')
 sources=[POS/'resolvent_campaign.json',POS/'frontier_cycle.json',POS/'frontier_star.json',
          POS/'frontier_horizons_cycle.json',POS/'frontier_horizons_star.json',
          POS/'cyclecraft_ledger_cycle.json',NON/'portloop14_case118.json',
          NON/'rational_access/rational_case118_pocket6.json',NON/'rational_access/rational_case118_hw20.json',
          NON/'two_strength_cert_results.json',NON/'two_strength_design_results.json']
 sources.append(NON/'two_strength_raw_design_results.json')
 provenance={'schema':'nonagon.possq.figures.v1','generator':str(Path(__file__).relative_to(BASE)),
             'sources':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources if p.exists()},
             'display':'Binary64 projections of retained rational/dyadic endpoints; scores use their retained display values.'}
 (OUT/'figure_provenance.json').write_text(json.dumps(provenance,indent=2))
 print(json.dumps({'figures':list(captions),'output':str(OUT)}))

if __name__=='__main__':main()
