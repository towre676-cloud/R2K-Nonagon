"""Direct Arb argument-principle and Chebyshev-Hermite contour engine.

No archive Taylor/DFT module is imported. Arb's direct adaptive integrator and
function ball evaluations are trusted. A separate dyadic partition proves that
the complete rectangular boundary is zero-free before contour integration.
"""
import hashlib
import importlib.metadata
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPS = HERE.parents[1] / 'nonagon_deps'
if DEPS.exists():
    sys.path.insert(0, str(DEPS))
from flint import acb, arb, acb_series, arb_mat, ctx


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()


def sha(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def rational(text):
    q = Fraction(str(text))
    return arb(q.numerator) / q.denominator


def serialize(value):
    def exact(component):
        if not component.is_finite():
            return None
        m,e=component.mid().man_exp()
        r,f=component.rad().man_exp()
        return {'mid_mantissa':str(m),'mid_exponent':int(e),
                'rad_mantissa':str(r),'rad_exponent':int(f)}
    return {'real': value.real.str(35, radius=True),
            'imag': value.imag.str(35, radius=True),
            'exact_real':exact(value.real),'exact_imag':exact(value.imag)}


def deserialize_exact(value):
    def part(data):
        require(data is not None,'nonfinite exact component')
        mid=arb(int(data['mid_mantissa']))*arb(2)**data['mid_exponent']
        rad=arb(int(data['rad_mantissa']))*arb(2)**data['rad_exponent']
        return arb(mid,rad)
    return acb(part(value['exact_real']),part(value['exact_imag']))


class AnalyticFunction:
    def __init__(self, spec):
        self.spec = spec
        self.kind = spec['kind']
        self.value_calls = 0
        self.derivative_calls = 0
        if self.kind == 'polynomial':
            self.coeff = [acb(rational(x)) for x in spec['coefficients_ascending']]
            require(len(self.coeff) >= 2, 'nonconstant polynomial required')
            self.deriv = [j*self.coeff[j] for j in range(1, len(self.coeff))]
        else:
            require(self.kind == 'riemann_xi_t', 'unsupported analytic function')
            self.I = acb(0, 1)
            self.logpi = arb.pi().log()

    @staticmethod
    def horner(coeff, z):
        value = acb(0)
        for c in coeff[::-1]:
            value = value*z+c
        return value

    def boundary_value(self, z):
        self.value_calls += 1
        if self.kind == 'polynomial':
            return self.horner(self.coeff, z)
        # At Re(t)>1 the elementary and Gamma factors completing zeta are
        # nonzero and finite. Thus proving zeta nonzero proves Xi nonzero.
        return (acb('0.5')+self.I*z).zeta()

    def logarithmic_derivative(self, z):
        self.derivative_calls += 1
        if self.kind == 'polynomial':
            return self.horner(self.deriv,z)/self.horner(self.coeff,z)
        require(ctx.cap >= 2, 'series precision cap below derivative order')
        s = acb('0.5')+self.I*z
        series = acb_series([s,1],prec=2).zeta()
        require(series.prec >= 2, 'zeta derivative series was truncated')
        return self.I*(1/s+1/(s-1)-self.logpi/2+(s/2).digamma()/2+series[1]/series[0])


def boundary_certificate(fun, corners, settings):
    """Full boundary coverage by exact-dyadic subintervals excluding zero."""
    records = []
    max_depth = settings['boundary_max_depth']
    limit = settings['boundary_eval_limit']
    start = fun.value_calls
    for side in range(4):
        A, D = corners[side], corners[(side+1)%4]-corners[side]
        # Explicit point checks detect exact endpoint/midpoint boundary zeros
        # in controls without waiting for an integrator budget to be exhausted.
        for tau in [Fraction(0),Fraction(1,2),Fraction(1)]:
            v = fun.boundary_value(A+D*rational(tau))
            if not v.is_finite() or v.contains(0):
                return {'proved_zero_free':False,'reason':'boundary point not excluded from zero',
                        'side':side,'tau':str(tau),'value':serialize(v),'segments':records}
        stack = [(0,1,0)]
        while stack:
            left,right,depth=stack.pop()
            denominator=2**depth
            q0,q1=Fraction(left,denominator),Fraction(right,denominator)
            mid,rad=(q0+q1)/2,(q1-q0)/2
            tau=arb(rational(mid).mid(),rational(rad).upper())
            z=A+D*tau
            value=fun.boundary_value(z)
            if value.is_finite() and not value.contains(0):
                records.append({'side':side,'tau0':str(q0),'tau1':str(q1),
                                'function_value':serialize(value),
                                'abs_lower':value.abs_lower().str(30,radius=True)})
                continue
            if depth>=max_depth or fun.value_calls-start>=limit:
                return {'proved_zero_free':False,'reason':'boundary exclusion budget exhausted',
                        'side':side,'tau0':str(q0),'tau1':str(q1),'depth':depth,
                        'segments':records}
            # Right is pushed first so accepted records have increasing tau.
            stack.append((2*left+1,2*right,depth+1))
            stack.append((2*left,2*left+1,depth+1))
    # Check exact-dyadic coverage, independent of traversal assumptions.
    for side in range(4):
        selected=[r for r in records if r['side']==side]
        end=Fraction(0)
        for rec in selected:
            require(Fraction(rec['tau0'])==end,'boundary partition gap')
            end=Fraction(rec['tau1'])
        require(end==1,'boundary partition incomplete')
    return {'proved_zero_free':True,'segments':records,
            'value_evaluations':fun.value_calls-start}


def tile(spec, geometry, settings):
    """Return an interval certificate or an explicit unresolved result."""
    ctx.prec=settings['precision_bits']
    ctx.cap=8
    require(ctx.cap>=2,'derivative cap')
    require(settings['precision_bits']>=80,'precision too small for this protocol')
    started=time.perf_counter()
    loq,hiq,bq=[Fraction(str(geometry[k])) for k in ['lo','hi','b']]
    require(loq<hiq and bq>0,'invalid rectangle')
    if spec['kind']=='riemann_xi_t':
        require(loq>1,'Xi prefactor nonzero argument requires Re(t)>1')
    lo,hi,b=map(rational,[loq,hiq,bq])
    center,a=(lo+hi)/2,(hi-lo)/2
    corners=[acb(lo,-b),acb(hi,-b),acb(hi,b),acb(lo,b)]
    fun=AnalyticFunction(spec)
    out={'function':spec,'geometry':{k:str(v) for k,v in zip(['lo','hi','b'],[loq,hiq,bq])},
         'input_sha256':sha({'function':spec,'geometry':geometry,'settings':settings}),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'precision_bits':ctx.prec,'series_cap':ctx.cap,'method':'direct-acb-integral-full-contour',
         'status':'unresolved','moments':[],'side_integrals':[]}
    try:
        boundary=boundary_certificate(fun,corners,settings)
        out['boundary']=boundary
        if not boundary['proved_zero_free']:
            out['reason']=boundary['reason'];return out
        mu=[]
        for order in range(2*settings['max_count']-1):
            total=acb(0);sides=[]
            for side in range(4):
                A,D=corners[side],corners[(side+1)%4]-corners[side]
                def callback(t,analytic,A=A,D=D):
                    # These logarithmic derivatives are meromorphic; interval
                    # pole evaluations are nonfinite, as required by Arb's
                    # analytic=True contract. There are no branch-cut methods.
                    z=A+D*t
                    return ((z-center)/a).chebyshev_t(order)*fun.logarithmic_derivative(z)*D
                value=acb.integral(callback,0,1,rel_tol=arb(2)**-settings['tolerance_bits'],
                    abs_tol=arb(2)**-settings['tolerance_bits'],eval_limit=settings['integral_eval_limit'],
                    depth_limit=settings['integral_depth_limit'])
                if not value.is_finite():
                    out['reason']='direct integrator returned nonfinite ball';return out
                sides.append(serialize(value));total+=value
            value=total/(2*acb.pi()*acb(0,1));mu.append(value)
            out['moments'].append(serialize(value));out['side_integrals'].append(sides)
            if order==0:
                n=int(round(float(value.real.mid())))
                if not (0<=n<=settings['max_count'] and value.real.contains(n)
                        and (value.real-n).abs_upper()<arb('0.25') and value.imag.contains(0)
                        and value.imag.abs_upper()<arb('0.25')):
                    out['reason']='count does not isolate and contain a nonnegative integer';return out
                out['count']=n
                if n==0:
                    out['status']='certified_empty';out['all_real_simple']=True;return out
            if order>=2*n-2:break
        if not all(v.imag.contains(0) and v.imag.abs_upper()<arb('0.25') for v in mu):
            out['reason']='moments not enclosed as real';return out
        m=[v.real for v in mu]
        M=arb_mat(n,n,[(m[i+j]+m[abs(i-j)])/2 for i in range(n) for j in range(n)])
        minors=[arb_mat(k,k,[M[i,j] for i in range(k) for j in range(k)]).det() for k in range(1,n+1)]
        out['leading_minors']=[v.str(35,radius=True) for v in minors]
        if all(v>0 for v in minors):
            out['status']='certified_real_simple';out['all_real_simple']=True
        elif any(v<0 for v in minors):
            out['status']='certified_not_all_real_simple';out['all_real_simple']=False
        else:
            out['reason']='positive definiteness undecided';out['all_real_simple']=None
        if n==1 and out['status']=='certified_real_simple':
            # Count-one tiles require an extra moment to locate their zero;
            # this campaign only asks for counts and line/simple certification.
            out['zero_location']='not requested'
        if n==2:
            determinant=M.det()
            out['gram_determinant']=determinant.str(35,radius=True)
            if determinant>0:out['gap_from_determinant']=(a*determinant.sqrt()).str(35,radius=True)
        return out
    except (ValueError,RuntimeError,ZeroDivisionError,OverflowError) as error:
        out['reason']=str(error);return out
    finally:
        out['seconds']=time.perf_counter()-started
        out['boundary_value_calls']=fun.value_calls
        out['logarithmic_derivative_calls']=fun.derivative_calls
