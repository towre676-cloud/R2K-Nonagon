"""Rigorous rational interval inference from one-port local-load observations.

The physical class is a connected grounded spring/conductance network in known
coordinates, with externally enclosed normalized responses and independently
calibrated positive diagonal-load intervals. Bounds must be supplied externally.
The program cannot establish instrument calibration from the observations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
import argparse
import json
from pathlib import Path


class Unresolved(ValueError):
    """A sufficient arithmetic margin could not be established."""


def rational(x):
    if isinstance(x, (float, bool)):
        raise TypeError("Use rational strings or integers, never binary floats")
    return x if isinstance(x, F) else F(x)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        object.__setattr__(self, 'lo', rational(self.lo))
        object.__setattr__(self, 'hi', rational(self.hi))
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    @staticmethod
    def point(x):
        x = rational(x)
        return Interval(x, x)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Interval) else Interval.point(x)

    def __add__(self, other):
        b = self.coerce(other)
        return Interval(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        b = self.coerce(other)
        p = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return Interval(min(p), max(p))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise Unresolved('A denominator interval contains zero')
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def square(self):
        if self.lo <= 0 <= self.hi:
            return Interval(F(0), max(self.lo*self.lo, self.hi*self.hi))
        return Interval(min(self.lo*self.lo, self.hi*self.hi),
                        max(self.lo*self.lo, self.hi*self.hi))

    def sqrt(self, bits=96):
        if self.lo < 0:
            raise Unresolved('The discriminant interval is not nonnegative')
        def bounds(q):
            a, b = isqrt(q.numerator), isqrt(q.denominator)
            if a*a == q.numerator and b*b == q.denominator:
                v = F(a, b)
                return v, v
            scale = 1 << bits
            lower = F(isqrt((q.numerator*scale*scale)//q.denominator), scale)
            return lower, lower + F(1, scale)
        return Interval(bounds(self.lo)[0], bounds(self.hi)[1])

    def strictly_positive(self, name):
        if self.lo <= 0:
            raise Unresolved(name + ' is not separated from zero')
        return self

    def positive_class_intersection(self, name):
        if self.hi <= 0:
            raise Unresolved(name + ' is incompatible with the positive-network class')
        return Interval(max(F(0), self.lo), self.hi)

    def dyadic(self,bits=96):
        # Outward compression avoids denominator explosion across matrix entries.
        # Exact point intervals are retained without widening.
        if self.lo==self.hi:return self
        scale=1<<bits
        lower=F((self.lo.numerator*scale)//self.lo.denominator,scale)
        upper=F(-((-self.hi.numerator*scale)//self.hi.denominator),scale)
        return Interval(lower,upper)

    def contains(self, x):
        return self.lo <= rational(x) <= self.hi

    def json(self):
        return [str(self.lo), str(self.hi)]

    @property
    def midpoint(self):
        return (self.lo+self.hi)/2

    @property
    def radius(self):
        return (self.hi-self.lo)/2


def inverse(matrix):
    """Exact Gauss-Jordan implementation local to this bounded module."""
    n = len(matrix)
    a = [[rational(v) for v in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    if any(len(row) != 2*n for row in a):
        raise ValueError('Matrix must be square')
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            raise Unresolved('The midpoint compliance is singular')
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][k]
        a[k] = [x/p for x in a[k]]
        for i in range(n):
            if i != k:
                q = a[i][k]
                a[i] = [x-q*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def load_key(loads):
    return tuple(sorted((int(k), rational(v)) for k, v in loads.items() if rational(v)))


def _infer_bounds(n, observations, t1=F(1), t2=F(2), gap=None, bits=96, relative_load_error=F(0)):
    """Return conditional enclosures, or an explicit unresolved verdict.

    observations maps canonical tuples ((zero_based_site, rational_load), ...)
    to Interval scalar compliance readings. No hidden stiffness is accepted.
    """
    t1, t2 = rational(t1), rational(t2)
    if type(bits) is not int or not 1 <= bits <= 512:
        raise ValueError('Square-root enclosure bits must lie in 1..512')
    if type(n) is not int or n < 1 or not (0 < t1 < t2):
        raise ValueError('Require known n >= 1 and 0 < t1 < t2')
    if gap is not None:
        gap = rational(gap)
        if gap <= 0:
            raise ValueError('The optional minimum positive spring gap must be positive')
    relative_load_error=rational(relative_load_error)
    if not 0<=relative_load_error<1:
        raise ValueError('Relative load error must satisfy 0 <= e < 1')
    def actual_load(nominal):
        return Interval(nominal*(1-relative_load_error),nominal*(1+relative_load_error))
    contract = {
        'known_labeled_sites': n,
        'force_and_readout_site_zero_based': 0,
        'response_normalization': 'external intervals enclose true displacement/force or voltage/current ratio',
        'local_grounding_loads': 'independent relative calibration intervals',
        'relative_load_error': str(relative_load_error),
        'measurement_intervals': 'externally justified simultaneous bounds',
        'physical_class': 'connected grounded nonnegative spring/conductance network',
        'minimum_positive_edge_weight': None if gap is None else str(gap),
        'sqrt_enclosure_bits': bits,
    }
    try:
        def read(loads):
            key = load_key(loads)
            if key not in observations:
                raise ValueError('Missing load configuration: ' + str(key))
            return Interval.coerce(observations[key])
        h = read({}).strictly_positive('Baseline compliance').dyadic(bits)
        g = [[None for _ in range(n)] for _ in range(n)]
        g[0][0] = h
        a2 = {}
        for j in range(1, n):
            d1 = (h-read({j:t1})).strictly_positive('First load difference')
            d2 = (h-read({j:t2})).strictly_positive('Second load difference')
            u1,u2=actual_load(t1),actual_load(t2)
            separation=(u2-u1).strictly_positive('Single-load calibration separation')
            y1, y2 = u1/d1, u2/d2
            slope = (y2-y1)/separation
            intercept = (y1-u1*slope).strictly_positive('Single-load intercept')
            a2[j] = 1/intercept
            g[0][j] = g[j][0] = a2[j].sqrt(bits).dyadic(bits)
            g[j][j] = (slope/intercept).positive_class_intersection('Compliance diagonal').dyadic(bits)
        for j in range(1, n):
            for k in range(j+1, n):
                delta = (h-read({j:t1,k:t1})).strictly_positive('Pair-load difference')
                product = g[0][j]*g[0][k]
                # Each paired load has its own calibration interval; equality
                # of the two physical shunts is not assumed.
                u,v=actual_load(t1),actual_load(t1)
                aa = a2[j]*(1+v*g[k][k])/v + a2[k]*(1+u*g[j][j])/u
                bb = (1+u*g[j][j])*(1+v*g[k][k])
                disc = product.square()-delta*aa+delta.square()*bb/(u*v)
                c = (product-disc.sqrt(bits))/delta
                g[j][k] = g[k][j] = c.positive_class_intersection('Off-diagonal compliance').dyadic(bits)
        midpoint = [[v.midpoint for v in row] for row in g]
        radii = [[v.radius for v in row] for row in g]
        r = inverse(midpoint)
        rnorm = max(sum(abs(v) for v in row) for row in r)
        enorm = max(sum(row) for row in radii)
        rho = rnorm*enorm
        if rho >= 1:
            return {'format':'nonagon-bounded-intervention/1','status':'unresolved',
                    'reason':'Inverse perturbation margin is not below one',
                    'contract':contract,'G_intervals':[[v.json() for v in row] for row in g],
                    'inverse_margin':str(rho)}
        error = rnorm*rho/(1-rho)
        kbox = [[Interval(v-error,v+error).dyadic(bits) for v in row] for row in r]
        edges=[]
        for i in range(n):
            for j in range(i+1,n):
                raw = -kbox[i][j]
                if raw.hi < 0:
                    raise Unresolved('A spring interval contradicts the physical class')
                # The physical class supplies nonnegative spring weights.
                weight = Interval(max(F(0),raw.lo),max(F(0),raw.hi))
                status='unresolved'
                if weight.lo > 0:
                    status='present'
                elif gap is not None and weight.hi < gap:
                    status='absent_under_gap'
                if gap is not None and 0 < weight.lo and weight.hi < gap:
                    status='inconsistent_with_gap'
                edges.append({'sites_zero_based':[i,j],'weight_interval':weight.json(),'status':status})
        return {'format':'nonagon-bounded-intervention/1','status':'enclosed',
                'contract':contract,'G_intervals':[[v.json() for v in row] for row in g],
                'K_intervals':[[v.json() for v in row] for row in kbox],
                'inverse_witness':{'midpoint_inverse':[[str(v) for v in row] for row in r],
                                   'midpoint_inverse_norm_inf':str(rnorm),
                                   'compliance_radius_norm_inf':str(enorm),
                                   'rho':str(rho),'stiffness_entry_error':str(error)},
                'edges':edges}
    except Unresolved as exc:
        return {'format':'nonagon-bounded-intervention/1','status':'unresolved',
                'reason':str(exc),'contract':contract}


def infer_bounds(n, observations, t1=F(1), t2=F(2), gap=None, bits=96, relative_load_error=F(0)):
    """Return a self-contained interval record, with exact input enclosures."""
    result=_infer_bounds(n,observations,t1,t2,gap,bits,relative_load_error)
    result['nominal_load_amplitudes']=[str(rational(t1)),str(rational(t2))]
    result['acquisition_intervals']=[
        {'nominal_loads':{str(j):str(v) for j,v in key},
         'normalized_response_interval':Interval.coerce(value).json()}
        for key,value in sorted(observations.items())]
    return result


def canonical_plan(n,t1=F(1),t2=F(2)):
    result=[{}]
    for j in range(1,n):
        result.extend([{j:t1},{j:t2}])
    result.extend({j:t1,k:t1} for j in range(1,n) for k in range(j+1,n))
    return result


def observation_intervals(n,oracle,eta,t1=F(1),t2=F(2)):
    eta=rational(eta)
    if eta<0:raise ValueError('eta must be nonnegative')
    return {load_key(x):Interval(oracle(x)-eta,oracle(x)+eta)
            for x in canonical_plan(n,t1,t2)}


def bound_document(document, eta, gap=None, bits=96, relative_load_error=F(0)):
    """Adapt a canonical exact-center acquisition document to interval inputs.

    eta is an external simultaneous absolute error bound. The incoming scalar
    values are interpreted as centers; their exact-center metadata is not used
    to conclude that the underlying physical readings are exact.
    """
    if document.get('schema') != 'nonagon.intervention.observations.v1':
        raise ValueError('Unsupported observation document schema')
    meta=document.get('contract',{})
    n=meta.get('network_order')
    if type(n) is not int or not 2 <= n <= 100:
        raise ValueError('Known network order 2..100 required')
    expected={'input_coordinate':0,'readout_coordinate':0,
              'stiffness_class':'connected_grounded_spring','known_site_set':True,
              'force_amplitude':'1','sensor_gain':'1','calibration':'exact'}
    for key,value in expected.items():
        if type(meta.get(key)) is not type(value) or meta.get(key)!=value:
            raise ValueError('Unsupported calibration/access contract: '+key)
    labels=meta.get('coordinate_labels',[])
    if len(labels)!=n or len(set(labels))!=n:
        raise ValueError('Known unique site labels required')
    t=rational(document.get('plan',{}).get('base_load'))
    eta=rational(eta)
    if eta<0:raise ValueError('eta must be nonnegative')
    obs={}
    rows=document.get('observations',[])
    plan=canonical_plan(n,t,2*t)
    if len(rows)!=len(plan):raise ValueError('Complete finite load plan required')
    for row,loads in zip(rows,plan):
        key=load_key(row['loads'])
        if key!=load_key(loads):raise ValueError('Incorrect or reordered load configuration')
        value=rational(row['response'])
        obs[key]=Interval(value-eta,value+eta)
    result=infer_bounds(n,obs,t,2*t,gap,bits,relative_load_error)
    result['measurement_error_bound']=str(eta)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',help='JSON with n, eta and scalar observations [{loads:{site:load},response:...}]')
    parser.add_argument('--output')
    parser.add_argument('--eta', help='Externally justified simultaneous absolute scalar error bound')
    parser.add_argument('--gap', help='Declared minimum positive edge weight')
    parser.add_argument('--relative-load-error',default='0',help='Externally justified relative conductance/stiffness load error')
    args=parser.parse_args()
    doc=json.loads(Path(args.input).read_text())
    eta=rational(args.eta if args.eta is not None else doc.get('eta','0'))
    if args.gap is not None:doc['gap']=args.gap
    if doc.get('schema')=='nonagon.intervention.observations.v1':
        result=bound_document(doc,eta,doc.get('gap'),relative_load_error=rational(args.relative_load_error))
    else:
        obs={load_key(x['loads']):Interval(rational(x['response'])-eta,rational(x['response'])+eta)
             for x in doc['observations']}
        result=infer_bounds(doc['n'],obs,rational(doc.get('t1','1')),rational(doc.get('t2','2')),doc.get('gap'),relative_load_error=rational(args.relative_load_error))
    out=json.dumps(result,indent=2)
    if args.output:Path(args.output).write_text(out+'\n')
    else:print(out)

if __name__=='__main__':main()
