"""Observed-pilot gate; no scenario labels, clean data or true parameters."""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'adaptive'))
import numpy as np
import acquisition as adaptive
import reference_core as core

def pilot_features(omega,noisy_values,order,band=(.03,8.)):
    s=1j*np.asarray(omega);y=np.asarray(noisy_values)
    try:
        model=core.fit('loewner_plus_vf',s,y,order,band)
        other=core.fit('fixed_seed_vf',s,y,order,band)
        p=model['poles'];positive=p[p.imag>1e-8]
        gap=1.
        for i in range(len(positive)):
            for j in range(i):gap=min(gap,abs(positive[i]-positive[j])/max(abs(positive[i]),abs(positive[j]),1e-20))
        r=core.residue_fit(s,y,p);rr=float(np.min(abs(r))/max(np.max(abs(r)),1e-30))
        pool=1j*np.geomspace(*band,1025);a=model['predict'](pool);b=other['predict'](pool)
        disagreement=float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-30))
        damping=float(np.min(-p.real/np.maximum(abs(p),1e-30)))
        residual=float(np.linalg.norm(model['predict'](s)-y)/max(np.linalg.norm(y),1e-30))
        result=dict(valid=True,min_relative_pole_gap=float(gap),residue_ratio=rr,
                    model_disagreement=disagreement,min_damping_ratio=damping,training_nrmse=residual)
        if not all(np.isfinite(v) for k,v in result.items() if k!='valid'):raise ValueError('nonfinite features')
        return result
    except Exception as exc:
        return dict(valid=False,error=f'{type(exc).__name__}: {exc}')

def matches(features,rule):
    if rule['kind']=='always_log':return False
    if rule['kind']=='always_adaptive':return bool(features.get('valid'))
    if not features.get('valid'):return False
    if rule['kind']=='and':return all(matches(features,r) for r in rule['conditions'])
    value=features[rule['feature']]
    return value<=rule['threshold'] if rule['op']=='le' else value>=rule['threshold']

def choose_strategy(omega,noisy_values,order,band=(.03,8.),rule=None):
    if rule is None:rule=json.loads((Path(__file__).resolve().parent/'selected_rule.json').read_text())['rule']
    features=pilot_features(omega,noisy_values,order,band)
    return ('adaptive_disagreement' if matches(features,rule) else 'log_uniform'),features

def candidate_rules():
    rules=[dict(id='always_log',kind='always_log'),dict(id='always_adaptive',kind='always_adaptive')]
    specs=[('min_relative_pole_gap','le',[.005,.01,.02,.05,.1,.2,.4]),
           ('residue_ratio','le',[1e-6,1e-5,1e-4,.001,.01,.1]),
           ('model_disagreement','ge',[.001,.01,.1,1.]),
           ('training_nrmse','ge',[.0001,.001,.01,.05]),
           ('min_damping_ratio','le',[.005,.01,.025,.05])]
    for feature,op,thresholds in specs:
        for t in thresholds:rules.append(dict(id=f'{feature}_{op}_{t:g}',kind='threshold',feature=feature,op=op,threshold=t))
    for gap in [.01,.05,.1]:
        for disagreement in [.01,.1]:
            rules.append(dict(id=f'gap{gap:g}_and_disagreement{disagreement:g}',kind='and',conditions=[
                dict(kind='threshold',feature='min_relative_pole_gap',op='le',threshold=gap),
                dict(kind='threshold',feature='model_disagreement',op='ge',threshold=disagreement)]))
    return rules
