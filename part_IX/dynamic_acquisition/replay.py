#!/usr/bin/env python3
"""Finite unit-step records -> bounded real-Laplace observations -> K,C,M.

Synthetic acquisition demonstration for C=2*beta*M, diagonal M>0.  Tail and
quadrature envelopes are analytic; sensor, force normalization, and numerical
rounding envelopes are declared assumptions, not measurements of an instrument.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
from random import Random
import sys

import numpy as np
from scipy.integrate import simpson
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
CORE = HERE.parent
sys.path.insert(0, str(CORE))
from bounded_intervention import Interval, canonical_plan, infer_bounds, load_key


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solve(matrix, rhs):
    """Independent Fraction Gaussian forward solve; no inference calls."""
    n = len(matrix)
    w = [[F(x) for x in row] + [F(rhs[i])] for i, row in enumerate(matrix)]
    for j in range(n):
        p = next(i for i in range(j, n) if w[i][j])
        w[j], w[p] = w[p], w[j]
        d = w[j][j]
        w[j] = [x/d for x in w[j]]
        for i in range(n):
            if i != j:
                d = w[i][j]
                w[i] = [a-d*b for a,b in zip(w[i],w[j])]
    return [row[-1] for row in w]


def strings(matrix):
    return [[str(x) for x in row] for row in matrix]


def fixture(n, seed):
    rng = Random(seed)
    grounds = [F(rng.randint(8,12),4) for _ in range(n)]
    k = [[F(i == j)*grounds[i] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            w = F(rng.randint(4,8),4)
            k[i][i] += w; k[j][j] += w
            k[i][j] = k[j][i] = -w
    masses = [F(4+i,4) for i in range(n)]
    m = [[F(i == j)*masses[i] for j in range(n)] for i in range(n)]
    beta = F(1,2)
    c = [[2*beta*x for x in row] for row in m]
    return dict(id=f'dynamic_n{n}_{seed}', n=n, seed=seed, K=k, C=c, M=m,
                beta=beta, grounds=grounds, masses=masses)


def validate(f):
    n=f['n']; k=f['K']; c=f['C']; m=f['M']; beta=f['beta']
    if beta <= 0: raise ValueError('beta must be positive')
    if any(m[i][j] != 0 for i in range(n) for j in range(n) if i != j):
        raise ValueError('Tail proof requires diagonal mass')
    if any(m[i][i] <= 0 for i in range(n)):
        raise ValueError('Mass must be strictly positive')
    if any(c[i][j] != 2*beta*m[i][j] for i in range(n) for j in range(n)):
        raise ValueError('This tail proof requires C=2*beta*M')
    if any(k[i][j] != k[j][i] or (i != j and k[i][j] > 0)
           for i in range(n) for j in range(n)):
        raise ValueError('Symmetric grounded spring matrix required')
    g = [sum(row) for row in k]
    if min(g) <= 0:
        raise ValueError('This tail proof requires strictly positive grounding at every site')
    ell = min(g)/max(m[i][i] for i in range(n))
    if ell <= beta*beta:
        raise ValueError('This tail proof requires underdamped spectral lower margin')
    return ell


def sqrt_lower(q, bits=64):
    scale=1 << bits
    return F(isqrt(q.numerator*scale*scale//q.denominator),scale)


def exp_negative_upper(q):
    # e > 2, so exp(-q) <= 2^(-floor(q)) for every rational q >= 0.
    require(q >= 0, 'Nonnegative exponential argument required')
    return F(1, 1 << (q.numerator//q.denominator))


def loaded(f, loads):
    k=[row[:] for row in f['K']]
    for j,t in loads.items(): k[int(j)][int(j)] += F(t)
    return k


def operator(f,s):
    n=f['n']
    return [[f['K'][i][j]+s*f['C'][i][j]+s*s*f['M'][i][j]
             for j in range(n)] for i in range(n)]


def modal_samples(f, loads, time):
    """Finite step displacement at the same fixed port, using modal synthesis."""
    k=loaded(f,loads); n=f['n']; beta=float(f['beta'])
    masses=np.array([float(f['M'][i][i]) for i in range(n)])
    invroot=1/np.sqrt(masses)
    a=np.array(k,dtype=float)*invroot[:,None]*invroot[None,:]
    lam,v=np.linalg.eigh(a)
    omega=np.sqrt(lam-beta*beta)
    weights=v[0,:]**2/masses[0]
    response=np.zeros_like(time)
    for q,w,z in zip(lam,weights,omega):
        response += w/q*(1-np.exp(-beta*time)*(np.cos(z*time)+beta/z*np.sin(z*time)))
    return response, (lam,v,omega,weights)


def state_sample(f,loads,time):
    """Separate first-order matrix-exponential step solution for spot checks."""
    n=f['n']; k=np.array(loaded(f,loads),dtype=float)
    m=np.array(f['M'],dtype=float); c=np.array(f['C'],dtype=float)
    a=np.block([[np.zeros((n,n)),np.eye(n)],[-np.linalg.solve(m,k),-np.linalg.solve(m,c)]])
    z=np.zeros(2*n); z[:n]=np.linalg.solve(k,np.eye(n)[0])
    # x(t)=x_equilibrium - exp(A t) x_equilibrium, with zero initial state.
    return (z-expm(a*time)@z)[0]


def validate_prior(f, prior):
    """Offline fixture qualification against externally supplied class bounds."""
    validate(f)
    if f['beta'] != F(prior['beta']):
        raise ValueError('Fixture decay does not satisfy declared beta')
    masses=[f['M'][i][i] for i in range(f['n'])]
    if min(masses)<F(prior['mass_min']) or max(masses)>F(prior['mass_max']):
        raise ValueError('Fixture mass outside declared class')
    if min(sum(row) for row in f['K'])<F(prior['ground_min']):
        raise ValueError('Fixture grounding outside declared class')
    if max(sum(abs(x) for x in row) for row in f['K'])>F(prior['stiffness_norm_inf_max']):
        raise ValueError('Fixture stiffness outside declared class')


def envelopes(prior, loads, s, horizon, step, sample_error, rounding_error):
    """Bounds depend only on declared class metadata and commanded loads."""
    ell=F(prior['ground_min'])/F(prior['mass_max'])
    beta=F(prior['beta'])
    require(ell>beta*beta,'Declared class needs a positive underdamped margin')
    omega_min=sqrt_lower(ell-beta*beta)
    # Weighted Rayleigh quotient: A<=||K||_inf/m_min*I.
    upper=(F(prior['stiffness_norm_inf_max'])+max((F(x) for x in loads.values()),default=F(0)))/F(prior['mass_min'])
    omega_sq_max=upper-beta*beta
    amplitude=1/(F(prior['mass_min'])*ell)
    if s == 0:
        tail=amplitude*(1+beta/omega_min)*exp_negative_upper(beta*horizon)
        quadrature=F(0); observation=sample_error
    else:
        a=s+beta
        tail=amplitude*exp_negative_upper(s*horizon)*(1+(1+beta/omega_min)*exp_negative_upper(beta*horizon))
        fourth=amplitude*(s**4+(a*a+omega_sq_max)**2*(1+beta/omega_min))
        quadrature=s*horizon*step**4*fourth/180
        # Simpson has positive weights summing to T; exp(-s*t)<=1.
        observation=s*horizon*sample_error
    return dict(tail=tail, quadrature=quadrature, sample=observation,
                numerical=rounding_error,
                total=tail+quadrature+observation+rounding_error,
                spectral_lower=ell, spectral_upper=upper,
                omega_lower=omega_min)


def box_matrix(result):
    return [[Interval(F(x[0]),F(x[1])) for x in row] for row in result['K_intervals']]


def combine(a,b,c,ca,cb,cc):
    return [[ca*a[i][j]+cb*b[i][j]+cc*c[i][j] for j in range(len(a))] for i in range(len(a))]


def infer_dynamic(n, packets):
    # Observer receives only scalar intervals, known labels and loads.
    recovered=[]
    for s in (F(0),F(1),F(2)):
        observations={load_key(x['loads']):Interval(F(x['interval'][0]),F(x['interval'][1]))
                      for x in packets if F(x['s'])==s}
        result=infer_bounds(n,observations)
        if result['status'] != 'enclosed':
            return {'status':'unresolved','s':str(s),'reason':result.get('reason')}
        recovered.append(result)
    p0,p1,p2=map(box_matrix,recovered)
    masses=combine(p0,p1,p2,F(1,2),F(-1),F(1,2))
    dampers=combine(p0,p1,p2,F(-3,2),F(2),F(-1,2))
    return dict(status='enclosed', operators=recovered,
                K_intervals=[[x.json() for x in row] for row in p0],
                C_intervals=[[x.json() for x in row] for row in dampers],
                M_intervals=[[x.json() for x in row] for row in masses])


def replay(output_dir):
    output_dir.mkdir(parents=True,exist_ok=True)
    protocol={
        'schema':'nonagon.dynamic-step.protocol.v1',
        'excitation':'calibrated unit step force at port 0; zero initial state',
        'readout':'displacement at port 0',
        's':['0','1','2'], 'horizon':'80', 'sample_interval':'1/2048',
        'sample_absolute_error':'1/1000000000000',
        'transformed_numerical_error':'1/1000000000000',
        'load_relative_error':'0',
        'certification_prior':{
            'beta':'1/2','mass_min':'1','mass_max':'7/4',
            'ground_min':'2','stiffness_norm_inf_max':'16',
            'source':'fixed class bounds, supplied independently of scalar records and hidden matrices',
        },
        'inverse':'existing finite t=1,2t=2,pair t=1 interval protocol',
        'proof_scope':'C=2*beta*M, diagonal M in [1,7/4], ground>=2, norm_inf(K)<=16, beta=1/2',
        'calibration_scope':'declared simultaneous sensor error; exact step force and loads; no drift/model mismatch beyond this envelope',
        'fixture_seeds':[910301,910302,910401,910402],
        'no_physical_hardware':True,
    }
    horizon=F(protocol['horizon']); step=F(protocol['sample_interval'])
    sample_error=F(protocol['sample_absolute_error']); rounding_error=F(protocol['transformed_numerical_error'])
    prior=protocol['certification_prior']
    count=int(horizon/step); require(count%2==0,'Even Simpson subinterval count required')
    times=np.arange(count+1,dtype=float)*float(step)
    fixtures=[fixture(n,seed) for n,seed in [(3,910301),(3,910302),(4,910401),(4,910402)]]
    records=[]; summaries=[]; checks=0; max_error=0.; max_ratio=0.; max_state_error=0.; max_radius=0.
    for f in fixtures:
        validate_prior(f,prior); packets=[]; modal_errors=[]
        for loads in canonical_plan(f['n']):
            values, modal=modal_samples(f,loads,times)
            for time in (F(1,8),F(1),F(7,2)):
                index=int(time/step)
                error=abs(float(values[index])-state_sample(f,loads,float(time)))
                max_state_error=max(max_state_error,error); checks+=1
                require(error<1e-12,'Modal and state-space transient disagree')
            # Deterministic bounded sample errors vary in sign; same trace reused for s=1,2.
            phase=np.arange(count+1,dtype=np.int64)
            noise=np.where((phase%7)<3,1.,-1.)*float(sample_error)
            observed=values+noise
            for s in (F(0),F(1),F(2)):
                center=float(observed[-1]) if s==0 else float(s)*float(simpson(np.exp(-float(s)*times)*observed,dx=float(step)))
                budget=envelopes(prior,loads,s,horizon,step,sample_error,rounding_error)
                p=operator(f,s)
                for j,t in loads.items(): p[int(j)][int(j)]+=F(t)
                exact=solve(p,[F(i==0) for i in range(f['n'])])[0]
                # str(float) is an exact rational center of the submitted packet.
                rational_center=F(str(center)); error=abs(rational_center-exact)
                radius=budget['total']; checks+=1
                require(error<=radius,'Analytic observation envelope missed exact Laplace response')
                max_error=max(max_error,float(error)); max_ratio=max(max_ratio,float(error/radius))
                packet=dict(fixture=f['id'],s=str(s),loads={str(j):str(t) for j,t in loads.items()},
                            center=str(rational_center),interval=[str(rational_center-radius),str(rational_center+radius)],
                            error_budget={k:str(v) for k,v in budget.items()},exact_reference=str(exact),
                            observed_absolute_error=str(error))
                packets.append(packet)
        inferred=infer_dynamic(f['n'],packets)
        require(inferred['status']=='enclosed','Dynamic interval reconstruction unresolved')
        for name in ('K','C','M'):
            for i in range(f['n']):
                for j in range(f['n']):
                    lo,hi=map(F,inferred[name+'_intervals'][i][j]); truth=f[name][i][j]
                    checks+=1; require(lo<=truth<=hi,'Dynamic operator enclosure missed truth')
                    max_radius=max(max_radius,float((hi-lo)/2))
        # Independent exact static forward references establish nominal KCM identities.
        p0,p1,p2=[operator(f,s) for s in (F(0),F(1),F(2))]
        for i in range(f['n']):
            for j in range(f['n']):
                require((p2[i][j]-2*p1[i][j]+p0[i][j])/2==f['M'][i][j],'Mass identity failed')
                require((4*p1[i][j]-p2[i][j]-3*p0[i][j])/2==f['C'][i][j],'Damping identity failed')
                checks+=2
        records.append(dict(fixture=dict(id=f['id'],seed=f['seed'],K=strings(f['K']),C=strings(f['C']),M=strings(f['M']),beta=str(f['beta'])),
                            packets=packets,reconstruction=inferred))
        summaries.append(dict(id=f['id'],order=f['n'],traces=len(canonical_plan(f['n'])),
                              scalar_settings=len(packets),status=inferred['status']))
    controls=[]
    from copy import deepcopy
    short_time=np.arange(int(F(1)/step)+1,dtype=float)*float(step)
    short_values,_=modal_samples(fixtures[0],{},short_time)
    for s in (F(0),F(1),F(2)):
        short_center=(float(short_values[-1]) if s==0 else
                      float(s)*float(simpson(np.exp(-float(s)*short_time)*short_values,dx=float(step))))
        short_budget=envelopes(prior,{},s,F(1),step,sample_error,rounding_error)
        exact=solve(operator(fixtures[0],s),[F(i==0) for i in range(fixtures[0]['n'])])[0]
        error=abs(F(str(short_center))-exact)
        require(error<=short_budget['total'],'Short-record tail envelope missed true transform')
        require(error>short_budget['total']-short_budget['tail'],'Short-record test did not detect omitted tail')
        controls.append(dict(id='short_record_tail_required_s'+str(s),status='contained_with_tail',
                             observed_error=str(error),tail_bound=str(short_budget['tail'])))
        checks+=2
    def excessive_damping(f):
        f['beta']=F(10)
        f['C']=[[20*x for x in row] for row in f['M']]
    mutations=[('nonproportional_damping',lambda f:f['C'][0].__setitem__(0,f['C'][0][0]+1)),
               ('nondiagonal_mass',lambda f:f['M'][0].__setitem__(1,F(1,8))),
               ('negative_mass',lambda f:f['M'][0].__setitem__(0,F(-1))),
               ('nonpositive_decay',lambda f:f.__setitem__('beta',F(0))),
               ('missing_grounding',lambda f:f['K'][0].__setitem__(0,-sum(f['K'][0][1:]))),
               ('no_underdamped_margin',excessive_damping)]
    for name,mutate in mutations:
        bad=deepcopy(fixtures[0]); mutate(bad)
        try: validate(bad)
        except ValueError as exc: controls.append(dict(id=name,status='rejected',reason=str(exc))); checks+=1
        else: raise RuntimeError('Unsupported tail model accepted: '+name)
    # Increasing response uncertainty must not silently retain the qualified result.
    bad_packets=deepcopy(records[0]['packets'])
    for row in bad_packets:
        mid=F(row['center']); row['interval']=[str(mid-F(1,100)),str(mid+F(1,100))]
    wide=infer_dynamic(fixtures[0]['n'],bad_packets)
    require(wide['status']=='unresolved','Wide-uncertainty control should remain unresolved'); checks+=1
    controls.append(dict(id='wide_observation_error',status=wide['status'],reason=wide['reason']))
    (output_dir/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    (output_dir/'acquisition_records.json').write_text(json.dumps(records,indent=2)+'\n')
    receipt=dict(schema='nonagon.dynamic-step.receipt.v1',status='PASSED',
                 fixtures=len(fixtures),step_transient_records=sum(x['traces'] for x in summaries),
                 real_laplace_scalar_settings=sum(x['scalar_settings'] for x in summaries),
                 sampled_points_per_trace=count+1,operator_enclosures=3*len(fixtures),
                 checked_predicates=checks,controls=controls,summaries=summaries,
                 max_observation_error=float(max_error),max_error_fraction_of_envelope=float(max_ratio),
                 max_modal_vs_state_space_error=float(max_state_error),max_operator_entry_radius=float(max_radius),
                 scope=protocol['proof_scope'],no_physical_hardware=True,
                 envelope_metadata_source='fixed external class prior in protocol; no hidden matrix enters envelopes()',
                 conditional_assumptions=['exact calibrated unit step and local loads','declared simultaneous per-sample error','declared transformed numerical rounding envelope','linear time-invariant model and zero initial state'],
                 dependency_versions=dict(numpy=np.__version__,python=sys.version.split()[0]),
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 interval_source_sha256=hashlib.sha256((CORE/'bounded_intervention.py').read_bytes()).hexdigest(),
                 protocol_sha256=hashlib.sha256((output_dir/'protocol.json').read_bytes()).hexdigest(),
                 records_sha256=hashlib.sha256((output_dir/'acquisition_records.json').read_bytes()).hexdigest())
    (output_dir/'dynamic_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','fixtures','step_transient_records','real_laplace_scalar_settings','operator_enclosures','checked_predicates','max_observation_error','max_operator_entry_radius')},indent=2))
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    args=parser.parse_args(); replay(args.output_dir)
