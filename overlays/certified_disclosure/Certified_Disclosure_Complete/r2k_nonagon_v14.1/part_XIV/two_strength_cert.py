"""Exact interval elimination for two-strength pair tomography.

The certificate encloses every x compatible with the declared input intervals
and both pair equations. It does not assert that arbitrary inconsistent data
have a physical realization. No root is chosen by proximity or by a known x.

Run ``python two_strength_cert.py --campaign`` from any directory. The campaign
uses unmodified MATPOWER rational data, exact inverses, and deterministic bounded
reading perturbations. Output separates calibrated-input and raw-reading paths.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import argparse
import hashlib
import json
import sys
import time

# Exact inverses of public grid matrices legitimately have long numerators.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        object.__setattr__(self, 'lo', F(self.lo))
        object.__setattr__(self, 'hi', F(self.hi))
        if self.lo > self.hi:
            raise ValueError('empty interval')

    @classmethod
    def point(cls, x):
        return cls(F(x), F(x))

    @classmethod
    def ball(cls, x, eps):
        if F(eps) < 0:
            raise ValueError('negative error budget')
        return cls(F(x) - F(eps), F(x) + F(eps))

    def __add__(self, other):
        other = as_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -as_interval(other)

    def __rsub__(self, other):
        return as_interval(other) + -self

    def __mul__(self, other):
        other = as_interval(other)
        vals = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(vals), max(vals))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.contains(0):
            raise ZeroDivisionError('interval denominator contains zero')
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self * as_interval(other).reciprocal()

    def __rtruediv__(self, other):
        return as_interval(other) * self.reciprocal()

    def square(self):
        if self.contains(0):
            return Interval(0, max(self.lo*self.lo, self.hi*self.hi))
        return Interval(min(self.lo*self.lo, self.hi*self.hi), max(self.lo*self.lo, self.hi*self.hi))

    def sqrt(self, bits=160):
        if self.lo < 0:
            raise ValueError('square-root input not certified nonnegative')
        def floor_sqrt(x):
            q = (x.numerator << (2*bits)) // x.denominator
            return F(isqrt(q), 1 << bits)
        l, h = floor_sqrt(self.lo), floor_sqrt(self.hi)
        if h*h != self.hi:
            h += F(1, 1 << bits)
        return Interval(l, h)

    def contains(self, x):
        return self.lo <= F(x) <= self.hi

    def intersect(self, other):
        return Interval(max(self.lo, other.lo), min(self.hi, other.hi))

    def rounded(self, bits=192):
        """Outward dyadic rounding; keep singleton exact identities unchanged."""
        if self.lo == self.hi:
            return self
        scale = 1 << bits
        return Interval(F((self.lo*scale).__floor__(), scale),
                        F((self.hi*scale).__ceil__(), scale))

    @property
    def mid(self):
        return (self.lo + self.hi)/2

    @property
    def width(self):
        return self.hi-self.lo

    @property
    def distance_from_zero(self):
        return F(0) if self.contains(0) else min(abs(self.lo), abs(self.hi))

    def as_list(self):
        return [str(self.lo), str(self.hi)]


def as_interval(x):
    return x if isinstance(x, Interval) else Interval.point(x)


def reconstruct_single(theta, delta1, delta2, sign=1):
    """Enclose G_jj and G_0j from two raw single-node response differences.

    delta1=h0-h(theta), delta2=h0-h(2theta). The caller must provide the
    calibrated sign of G_0j. A grounded connected passive M-matrix has sign +1.
    The two delta intervals may share a baseline error; ignoring correlation is
    conservative. A denominator crossing zero causes explicit rejection.
    """
    theta = F(theta)
    if theta <= 0 or sign not in (-1, 1):
        raise ValueError('positive strength and calibrated sign required')
    d1, d2 = as_interval(delta1), as_interval(delta2)
    diagonal = (2*d1-d2)/(2*theta*(d2-d1))
    column_square = (d1*(1+theta*diagonal)/theta).intersect(d1*d2/(2*theta*(d2-d1)))
    column = column_square.sqrt()
    return diagonal, column if sign == 1 else -column


def pair_delta(theta, u, v, A, B, x):
    """Exact response difference, for exact rational arguments."""
    theta, u, v, A, B, x = map(F, (theta, u, v, A, B, x))
    aa, bb = 1/theta+A, 1/theta+B
    return (u*u*bb+v*v*aa-2*u*v*x)/(aa*bb-x*x)


def certify_pair(theta1, theta2, u, v, A, B, delta1, delta2, iterations=6, raw_interval=False):
    """Return a rational interval enclosing every compatible common root.

    All arguments after theta2 may be Interval or exact rationals. The exact
    uniqueness statement is conditional on u*v*(delta2-delta1) != 0 and on both
    loaded matrices being nonsingular; physical consistency is a separate issue.
    """
    theta1, theta2 = F(theta1), F(theta2)
    if theta1 == theta2 or theta1 == 0 or theta2 == 0:
        raise ValueError('two distinct nonzero strengths required')
    u, v, A, B, d1, d2 = map(as_interval, (u, v, A, B, delta1, delta2))
    r1, r2 = 1/theta1, 1/theta2
    aa1, bb1, aa2, bb2 = r1+A, r1+B, r2+A, r2+B
    a = u*v
    C1 = u.square()*bb1+v.square()*aa1-d1*aa1*bb1
    C2 = u.square()*bb2+v.square()*aa2-d2*aa2*bb2
    denominator = 2*a*(d2-d1)
    raw = (d2*C1-d1*C2)/denominator
    # Algebraic cancellation before interval evaluation reduces dependency loss.
    numerator = ((d2-d1)*(u.square()*B+v.square()*A)
                 +(d2*r1-d1*r2)*(u.square()+v.square())
                 -d1*d2*(r1-r2)*(r1+r2+A+B))
    reduced = numerator/denominator
    enclosure = raw.intersect(reduced).rounded()
    initial_width = enclosure.width
    for _ in range(iterations):
        old = enclosure
        for dl, C in ((d1, C1), (d2, C2)):
            derivative = 2*dl*enclosure-2*a
            if not derivative.contains(0):
                m = enclosure.mid
                residual_at_m = dl*m*m-2*a*m+C
                newton = m-residual_at_m/derivative
                enclosure = enclosure.intersect(newton).rounded()
        if enclosure == old:
            break
    residuals = [dl*enclosure.square()-2*a*enclosure+C for dl, C in ((d1, C1), (d2, C2))]
    if not all(r.contains(0) for r in residuals):
        raise ValueError('input box rejected by pair-equation residual')
    if raw_interval:
        return enclosure
    return {
        'schema': 'nonagon.two_strength_interval.v1',
        'theta': [str(theta1), str(theta2)],
        'inputs': {k: val.as_list() for k, val in zip(
            ('G0i', 'G0j', 'Gii', 'Gjj', 'delta1', 'delta2'), (u, v, A, B, d1, d2))},
        'enclosure': enclosure.as_list(),
        'denominator': denominator.as_list(),
        'denominator_margin': str(denominator.distance_from_zero),
        'initial_width': str(initial_width),
        'contraction_iterations_requested': iterations,
        'residuals': [r.as_list() for r in residuals],
        'arithmetic': 'exact rational endpoints',
        'claim': 'all physically compatible common roots lie in this interval',
    }


def verify_certificate(cert):
    """Deterministic replay of the exact arithmetic certificate."""
    try:
        inp = cert['inputs']
        data = [Interval(*map(F, inp[k])) for k in ('G0i','G0j','Gii','Gjj','delta1','delta2')]
        got = certify_pair(*map(F, cert['theta']), *data,
                           iterations=cert['contraction_iterations_requested'])
        return got == cert
    except (ValueError, ZeroDivisionError, KeyError, TypeError):
        return False


def certify_many(strengths, deltas, u, v, A, B, raw_interval=False):
    """Fuse all certified pair eliminations from three or more strengths.

    Selection and intersection depend only on observation intervals. Increasing
    this acquisition set can only tighten the resulting enclosure, provided
    each declared interval contains its true response difference.
    """
    strengths=list(map(F,strengths));deltas=list(map(as_interval,deltas))
    if len(strengths)!=len(deltas) or len(strengths)<2:
        raise ValueError('at least two matched strengths and observations required')
    if len(set(strengths))!=len(strengths):
        raise ValueError('strengths must be distinct')
    inputs=list(map(as_interval,(u,v,A,B)))
    enclosure=None;used=[];skipped=[]
    for i in range(len(strengths)):
        for j in range(i+1,len(strengths)):
            try:
                X=certify_pair(strengths[i],strengths[j],*inputs,deltas[i],deltas[j],raw_interval=True)
            except (ZeroDivisionError,ValueError):
                skipped.append([i,j]);continue
            enclosure=X if enclosure is None else enclosure.intersect(X)
            used.append([i,j])
    if enclosure is None:
        raise ValueError('no strength pair has a certified elimination denominator')
    if raw_interval:
        return enclosure
    return dict(schema='nonagon.multi_strength_interval.v1',strengths=list(map(str,strengths)),
                deltas=[d.as_list() for d in deltas],
                inputs={k:x.as_list() for k,x in zip(('G0i','G0j','Gii','Gjj'),inputs)},
                enclosure=enclosure.as_list(),used_pairs=used,skipped_pairs=skipped,
                claim='intersection enclosing every physically compatible common root')


def verify_many_certificate(cert):
    try:
        inp=cert['inputs']
        got=certify_many(list(map(F,cert['strengths'])),[Interval(*map(F,d)) for d in cert['deltas']],
                         *[Interval(*map(F,inp[k])) for k in ('G0i','G0j','Gii','Gjj')])
        return got==cert
    except (ValueError,ZeroDivisionError,KeyError,TypeError):
        return False


def certify_raw_pair(theta, single_i_1, single_i_2, single_j_1, single_j_2,
                     delta1, delta2, signs=(1,1), raw_interval=False):
    """Eliminate directly in six raw response differences, preserving algebra.

    This intersects the propagated-calibration certificate with a second exact
    expression that cancels the inferred diagonal entries symbolically. Signs
    refer to the calibrated port-column signs; passive M-matrices give (+,+).
    """
    theta=F(theta)
    si1,si2,sj1,sj2,d1,d2=map(as_interval,
        (single_i_1,single_i_2,single_j_1,single_j_2,delta1,delta2))
    Ai,ui=reconstruct_single(theta,si1,si2,signs[0])
    Bj,vj=reconstruct_single(theta,sj1,sj2,signs[1])
    propagated=certify_pair(theta,2*theta,ui,vj,Ai,Bj,d1,d2,raw_interval=True)
    uv_square=si1*si2*sj1*sj2/(4*theta*theta*(si2-si1)*(sj2-sj1))
    uv=uv_square.sqrt()*(signs[0]*signs[1])
    ratio=(d2*(si1+sj1-d1)/(si1*sj1)-d1*(si2+sj2-d2)/(si2*sj2))/(2*(d2-d1))
    direct=uv*ratio
    X=propagated.intersect(direct).rounded()
    if raw_interval:
        return X
    return dict(schema='nonagon.raw_two_strength_interval.v1',theta=str(theta),
                inputs={k:d.as_list() for k,d in zip(('si1','si2','sj1','sj2','delta1','delta2'),
                                                    (si1,si2,sj1,sj2,d1,d2))},
                signs=list(signs),enclosure=X.as_list(),propagated=propagated.as_list(),
                direct=direct.as_list(),claim='all compatible raw-reading models enclosed')


def verify_raw_certificate(cert):
    try:
        inp=cert['inputs']
        got=certify_raw_pair(F(cert['theta']),
                            *[Interval(*map(F,inp[k])) for k in ('si1','si2','sj1','sj2','delta1','delta2')],
                            signs=cert['signs'])
        return got==cert
    except (ValueError,ZeroDivisionError,KeyError,TypeError):
        return False


def load_grid(case):
    """Read exact model; return inverse and lean-pair indices without cwd aliases."""
    import flint
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root/'part_XI'))
    import gridlib
    path = root/'part_XI'/'data'/f'{case}.m'
    nodes, edges, ground, meta = gridlib.network(path)
    rec = json.loads((root/'part_XIII'/f'record_{case}.json').read_text())
    port = nodes.index(rec['port_bus'])
    n = len(nodes)
    K = flint.fmpq_mat(n, n)
    def fq(q): return flint.fmpq(q.numerator, q.denominator)
    for (i,j), w in edges.items():
        w=fq(w);K[i,i]+=w;K[j,j]+=w;K[i,j]-=w;K[j,i]-=w
    for i,g in enumerate(ground): K[i,i]+=fq(g)
    G=K.inv()
    def fr(q): return F(int(q.p), int(q.q))
    diagonal = [fr(G[i,i]) for i in range(n)]
    column = [fr(G[port,i]) for i in range(n)]
    adj = {i:set() for i in range(n)}
    for i,j in edges: adj[i].add(j);adj[j].add(i)
    pairs = {(min(i,j), max(i,j)) for i in range(n) for j in adj[i] if port not in (i,j)}
    idx = {b:i for i,b in enumerate(nodes)}
    widened = {idx[int(b)]:v for b,v in rec.get('widened_sources',{}).items()}
    for bus in rec['extra_source_buses']:
        s=idx[bus]; seen={s};front=[s]
        for _ in range(widened.get(s,1)+1):
            front=[o for u in front for o in adj[u] if o not in seen and not seen.add(o)]
        pairs.update((min(s,j),max(s,j)) for j in seen if j not in (port,s))
    samples=[(i,j,column[i],column[j],diagonal[i],diagonal[j],fr(G[i,j])) for i,j in sorted(pairs)]
    return samples, sorted(diagonal)[n//2], column[port], meta, hashlib.sha256(path.read_bytes()).hexdigest()


def campaign(cases=('case118','case300')):
    import numpy as np
    out={'schema':'nonagon.two_strength_campaign.v1','cases':[],
         'noise_model':'deterministic bounded response differences, radius 2*sigma*G00',
         'scope':'grid primary path conditions on exact calibrated G0i,G0j,Gii,Gjj; raw path propagates single-node uncertainty'}
    for case in cases:
        started=time.perf_counter()
        samples, med, G00, meta, sha=load_grid(case)
        rows=[]
        for c in (F(1,10), F(1), F(10), F(100)):
            theta=c/med
            for sigma in (F(0),F(1,10**10),F(1,10**8),F(1,10**6),F(1,10**4)):
                eps=2*sigma*abs(G00)
                widths=[]; successful=0;rejected=0;truth_contained=0;minus_fail=0; nearest_fail=0
                raw_success=0;raw_contained=0;raw_widths=[]
                for k,(i,j,u,v,A,B,x) in enumerate(samples):
                    d1=pair_delta(theta,u,v,A,B,x); d2=pair_delta(2*theta,u,v,A,B,x)
                    # Both signs occur and centres are perturbed within declared budgets.
                    e1=eps*F((k*17%9)-4,4);e2=eps*F((k*23%9)-4,4)
                    box1,box2=Interval.ball(d1+e1,eps),Interval.ball(d2+e2,eps)
                    try:
                        X=certify_pair(theta,2*theta,u,v,A,B,box1,box2,raw_interval=True)
                        successful+=1;truth_contained+=int(X.contains(x));widths.append(float(X.width/(2*abs(x))) if x else float('inf'))
                    except (ZeroDivisionError,ValueError): rejected+=1
                    # Point-rule comparison uses the same noisy centres, never the true x for selection.
                    def roots(t,d):
                        aa,bb=1/t+float(A),1/t+float(B);uf,vf=float(u),float(v)
                        disc=(d*aa-uf*uf)*(d*bb-vf*vf)
                        if disc<0 or d==0:return None
                        return np.array([(uf*vf-np.sqrt(disc))/d,(uf*vf+np.sqrt(disc))/d])
                    r1,r2=roots(float(theta),float(box1.mid)),roots(float(2*theta),float(box2.mid))
                    if r1 is None or r2 is None:
                        minus_fail+=1;nearest_fail+=1
                    else:
                        tol=1e-6*max(abs(float(x)),1e-300)
                        minus_fail+=int(abs(r1[0]-float(x))>tol)
                        pick=r1[np.argmin(np.min(abs(r1[:,None]-r2[None,:]),axis=1))]
                        nearest_fail+=int(abs(pick-float(x))>tol)
                    if case=='case118' and sigma in (F(1,10**8), F(1,10**6)):
                        # Passive M-matrix sign is known positive; raw single-node differences carry independent budgets.
                        try:
                            def single(t,gcol,gdiag):return t*gcol*gcol/(1+t*gdiag)
                            Ai,ui=reconstruct_single(theta,Interval.ball(single(theta,u,A),eps),Interval.ball(single(2*theta,u,A),eps))
                            Bj,vj=reconstruct_single(theta,Interval.ball(single(theta,v,B),eps),Interval.ball(single(2*theta,v,B),eps))
                            X=certify_pair(theta,2*theta,ui,vj,Ai,Bj,box1,box2,raw_interval=True)
                            raw_success+=1;raw_contained+=int(X.contains(x));raw_widths.append(float(X.width/(2*abs(x))))
                        except (ZeroDivisionError,ValueError):pass
                def quantile(vals,q):return float(np.quantile(vals,q)) if vals else None
                row=dict(strength_scale=str(c),sigma=str(sigma),pairs=len(samples),certified=successful,rejected=rejected,
                         truth_contained=truth_contained,median_relative_radius=quantile(widths,.5),p90_relative_radius=quantile(widths,.9),
                         minus_point_miss_at_1e6_relative=minus_fail,nearest_point_miss_at_1e6_relative=nearest_fail,
                         raw_single_calibration_certified=raw_success,raw_single_truth_contained=raw_contained,
                         raw_single_median_relative_radius=quantile(raw_widths,.5))
                assert truth_contained==successful and raw_contained==raw_success
                rows.append(row)
                print(case, str(c), str(sigma), 'cert', successful,'raw',raw_success,'radius',row['median_relative_radius'],flush=True)
        out['cases'].append(dict(case=case,pairs=len(samples),input_sha256=sha,negative_reactance_branches=meta['negative_reactance_branches'],seconds=time.perf_counter()-started,rows=rows))
    return out


def summarize_campaign(result):
    rows=[r for c in result['cases'] for r in c['rows']]
    return dict(settings=len(rows),pair_boxes=sum(r['pairs'] for r in rows),
                certified=sum(r['certified'] for r in rows),rejected=sum(r['rejected'] for r in rows),
                truth_covered=sum(r['truth_contained'] for r in rows),
                raw_calibration_certified=sum(r['raw_single_calibration_certified'] for r in rows),
                raw_calibration_truth_covered=sum(r['raw_single_truth_contained'] for r in rows),
                all_certified_cover_truth=all(r['certified']==r['truth_contained'] and
                    r['raw_single_calibration_certified']==r['raw_single_truth_contained'] for r in rows))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--campaign',action='store_true')
    parser.add_argument('--case',action='append')
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('two_strength_cert_results.json'))
    args=parser.parse_args()
    if args.campaign:
        result=campaign(tuple(args.case or ('case118','case300')))
        result['summary']=summarize_campaign(result)
        result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(args.output)
