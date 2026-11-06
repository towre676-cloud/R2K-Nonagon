#!/usr/bin/env python3
"""Fixed-split measured multi-cell impedance benchmark, with all candidates retained."""
from __future__ import annotations
import hashlib,json,time,warnings,platform,argparse
from pathlib import Path
import numpy as np
from scipy import linalg,interpolate,optimize
import scipy,skrf
from skrf.vectorFitting import VectorFitting
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def score(y,p):return float(linalg.norm(y-p)/linalg.norm(y))

def loewner(freq,z,order):
    scale=float(np.exp(np.mean(np.log(freq))))
    s=1j*freq/scale;ia=np.arange(0,len(s),2);ib=np.arange(1,len(s),2)
    mu=np.r_[s[ia],s[ia].conj()];nu=np.r_[s[ib],s[ib].conj()]
    v=np.r_[z[ia],z[ia].conj()];w=np.r_[z[ib],z[ib].conj()]
    L=(v[:,None]-w[None,:])/(mu[:,None]-nu[None,:])
    Ls=(mu[:,None]*v[:,None]-nu[None,:]*w[None,:])/(mu[:,None]-nu[None,:])
    if order>min(L.shape):raise ValueError('order exceeds pencil dimension')
    Y=linalg.svd(np.c_[L,Ls],full_matrices=False)[0][:,:order]
    X=linalg.svd(np.r_[L,Ls],full_matrices=False)[2].conj().T[:,:order]
    A=Y.conj().T@Ls@X;E=Y.conj().T@L@X;B=Y.conj().T@v;C=w@X
    poles=linalg.eigvals(A,E)*scale*2*np.pi
    def predict(f):return np.array([C@linalg.solve(A-(1j*x/scale)*E,B) for x in f])
    return predict,{'descriptor_order':order,'finite_rhp_poles':int(np.sum(np.isfinite(poles)&(poles.real>1e-8*scale)))},None

def vector_fitting(freq,z,real_poles,complex_pairs,max_iterations):
    ntwk=skrf.Network(frequency=skrf.Frequency.from_f(freq,unit='hz'),z=z[:,None,None],z0=50)
    vf=VectorFitting(ntwk);vf.max_iterations=max_iterations
    vf.vector_fit(n_poles_real=real_poles,n_poles_cmplx=complex_pairs,init_pole_spacing='log',parameter_type='z',fit_constant=True,fit_proportional=True,enforce_dc=False)
    return lambda f:vf.get_model_response(0,0,f),{'initial_real_poles':real_poles,'initial_complex_pairs':complex_pairs,'fitted_order':int(vf.get_model_order(vf.poles)),'finite_rhp_poles':int(np.sum(vf.poles.real>0))},None

def relaxation(freq,z,grid_size,ridge):
    tau=np.geomspace(.01/(2*np.pi*freq.max()),100/(2*np.pi*freq.min()),grid_size)
    def basis(f):
        s=2j*np.pi*np.asarray(f)
        return np.c_[np.ones(len(s)),s,1/s,1/(1+s[:,None]*tau)]
    B=basis(freq);weights=1/np.maximum(abs(z),np.finfo(float).tiny)
    A=np.r_[(B*weights[:,None]).real,(B*weights[:,None]).imag]
    y=np.r_[(z*weights).real,(z*weights).imag]
    norms=linalg.norm(A,axis=0)
    if np.any(norms==0):raise ValueError('zero design column')
    A=A/norms
    if ridge:
        penalty=np.diag(np.r_[np.zeros(3),np.ones(grid_size)])*np.sqrt(ridge)
        A=np.r_[A,penalty];y=np.r_[y,np.zeros(3+grid_size)]
    coef,residual=optimize.nnls(A,y,maxiter=10000)
    coef=coef/norms
    if not np.all(np.isfinite(coef)) or np.any(coef<0):raise ValueError('invalid NNLS coefficients')
    model={'representation':'positive_real_debye_v1','frequency_convention':'s=+i*2*pi*f_hz','units':{'s':'1/second','Z':'ohm','R0':'ohm','L':'henry','Cinv':'1/farad','a':'ohm','tau':'second'},'R0':repr(float(coef[0])),'L':repr(float(coef[1])),'Cinv':repr(float(coef[2])),'terms':[{'a':repr(float(a)),'tau':repr(float(t))} for a,t in zip(coef[3:],tau)]}
    return lambda f:basis(f)@coef,{'grid_size':grid_size,'ridge':ridge,'active_debye_terms':int(np.sum(coef[3:]>0)),'nonnegative_model_by_construction':True},model

def data_for(device,protocol,manifest):
    path=ROOT/'data'/f'{device}.csv'
    expected=next(x for x in manifest['files'] if x['name']==path.name)
    if sha(path)!=expected['sha256']:raise ValueError('dataset hash mismatch')
    d=np.genfromtxt(path,delimiter=',',skip_header=1)
    d=d[np.argsort(d[:,0])]
    if d.shape!=(protocol['expected_rows'],3) or not np.all(np.isfinite(d)) or np.any(d[:,0]<=0) or len(np.unique(d[:,0]))!=len(d):raise ValueError('invalid cell data')
    return d[:,0],d[:,1]+1j*d[:,2]

def run(phase):
    protocol=json.loads((ROOT/'protocol.json').read_text());manifest=json.loads((ROOT/'data/download_manifest.json').read_text())
    devices=protocol[phase+'_devices'];records=[];candidates=[];models=[];predictions=[]
    for device in devices:
        freq,z=data_for(device,protocol,manifest)
        for case in protocol['cases']:
            ti=np.array(case['train']);vi=np.array(case['validation']);qi=np.array(case['test'])
            if set(ti)&set(vi) or set(ti)&set(qi) or set(vi)&set(qi):raise ValueError('split leakage')
            for method in protocol['methods']:
                if method=='loewner':configs=[{'order':r} for r in protocol['loewner_orders']]
                elif method=='vector_fitting':configs=[{'real_poles':r,'complex_pairs':c,'max_iterations':protocol['vf_max_iterations']} for r in protocol['vf_real_poles'] for c in protocol['vf_complex_pairs']]
                elif method=='positive_real_relaxation':configs=[{'grid_size':g,'ridge':r} for g in protocol['relaxation_grid_sizes'] for r in protocol['relaxation_ridge']]
                else:configs=[{}]
                fitted=[]
                for cfg in configs:
                    start=time.perf_counter();rec={'phase':phase,'device':device,'scenario':case['scenario'],'budget':case['budget'],'method':method,'config':cfg}
                    try:
                        with warnings.catch_warnings(record=True) as wr:
                            warnings.simplefilter('always')
                            if method=='loewner':predict,info,model=loewner(freq[ti],z[ti],**cfg);complexity=cfg['order']
                            elif method=='vector_fitting':predict,info,model=vector_fitting(freq[ti],z[ti],**cfg);complexity=cfg['real_poles']+2*cfg['complex_pairs']
                            elif method=='positive_real_relaxation':predict,info,model=relaxation(freq[ti],z[ti],**cfg);complexity=cfg['grid_size']
                            else:
                                spline=interpolate.CubicSpline(np.log10(freq[ti]),z[ti]);predict=lambda f,spline=spline:spline(np.log10(f));info={'coordinate':'log10_frequency','boundary_condition':'not-a-knot'};model=None;complexity=0
                            vp=predict(freq[vi])
                            if not np.all(np.isfinite(vp)):raise ValueError('nonfinite validation predictions')
                            rec.update(status='ok',validation_error=score(z[vi],vp),model_info=info,warnings=[str(w.message) for w in wr])
                        fitted.append((rec['validation_error'],complexity,predict,rec,model))
                    except Exception as exc:rec.update(status='failed',error=type(exc).__name__+': '+str(exc))
                    rec['seconds']=time.perf_counter()-start;candidates.append(rec)
                base={'phase':phase,'device':device,'scenario':case['scenario'],'budget':case['budget'],'method':method,'n_train':len(ti),'n_validation':len(vi),'n_test':len(qi)}
                if not fitted:records.append(dict(base,status='all_candidates_failed'));continue
                _,_,predict,best,model=min(fitted,key=lambda x:(x[0],x[1]))
                try:
                    qp=predict(freq[qi])
                    if not np.all(np.isfinite(qp)):raise ValueError('nonfinite selected test predictions')
                except Exception as exc:records.append(dict(base,status='test_prediction_failed',error=str(exc)));continue
                out=dict(base,status='ok',selected_config=best['config'],model_info=best['model_info'],validation_error=best['validation_error'],test_error=score(z[qi],qp),test_max_absolute_error_ohm=float(np.max(abs(z[qi]-qp))),warnings=best['warnings'],selected_fit_validation_seconds=best['seconds'])
                if model is not None:
                    model_id=f'{device}_{case["scenario"]}_{case["budget"]}'
                    model.update(model_id=model_id,phase=phase,device=device,scenario=case['scenario'],budget=case['budget'])
                    models.append(model);out['passive_model_id']=model_id
                records.append(out)
                predictions.append(dict(base,frequency_hz=freq[qi].tolist(),measured_re=z[qi].real.tolist(),measured_im=z[qi].imag.tolist(),predicted_re=qp.real.tolist(),predicted_im=qp.imag.tolist()))
        print(phase,device,'complete',flush=True)
        # checkpoint completed devices without conflating partial status with final run.
        (ROOT/f'{phase}_checkpoint.json').write_text(json.dumps({'partial':True,'completed_device':device,'results':records},indent=2)+'\n')
    receipt={'phase':phase,'protocol_sha256':sha(ROOT/'protocol.json'),'source_sha256':sha(__file__),'download_manifest_sha256':sha(ROOT/'data/download_manifest.json'),'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'scikit_rf':skrf.__version__},'device_count':len(devices),'case_count':len(devices)*len(protocol['cases']),'selected_result_count':len(records),'candidate_count':len(candidates),'candidate_failure_count':sum(r['status']!='ok' for r in candidates),'selected_failure_count':sum(r['status']!='ok' for r in records),'results':records}
    for name,obj in [(f'{phase}_results.json',receipt),(f'{phase}_candidates.json',candidates),(f'{phase}_predictions.json',predictions),(f'{phase}_passive_models.json',{'models':models})]:
        (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    (ROOT/f'{phase}_checkpoint.json').unlink(missing_ok=True)
    allmodels=[]
    for p in [ROOT/'development_passive_models.json',ROOT/'confirmation_passive_models.json']:
        if p.exists():allmodels.extend(json.loads(p.read_text())['models'])
    (ROOT/'selected_passive_models.json').write_text(json.dumps({'models':allmodels},indent=2)+'\n')
    print('DONE',phase,'candidates',len(candidates),'selectedfailures',receipt['selected_failure_count'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['development','confirmation'],required=True);run(ap.parse_args().phase)
