"""Exact rational positive-real certificate, including lossless boundary poles.

This is an independent corrected implementation of the v10 Cayley/Brune branch.
It leaves the Claude source untouched. Coefficients are ascending: [a0,a1,...].
Use integers, Fraction, rational strings, or sympy Rational; ambiguous float
coefficients are deliberately rejected. SymPy is used only for exact QQ algebra,
square-free factorization, Sturm root counts, and exact determinants.
"""
from pathlib import Path
from fractions import Fraction
import json
import sympy as sp

s,x=sp.symbols('s x',real=True)
HERE=Path(__file__).resolve().parent

def rational(v):
    if isinstance(v,float):
        raise TypeError('Use an exact rational coefficient, not binary64 float')
    if isinstance(v,Fraction): return sp.Rational(v.numerator,v.denominator)
    q=sp.Rational(v)
    if not q.is_Rational: raise TypeError('Coefficient must be rational')
    return q

def polynomial(coefficients):
    return sp.Poly(sum(rational(c)*s**i for i,c in enumerate(coefficients)),s,domain=sp.QQ)

def strict_hurwitz(poly):
    """Routh-Hurwitz via all exact leading Hurwitz determinants.

    Any nonzero constant is Hurwitz. A negative overall polynomial scale has
    no effect; monic normalization makes the determinant test scale invariant.
    """
    if poly.is_zero: return dict(strictly_hurwitz=False,degree=None,reason='zero polynomial',determinants=[])
    n=poly.degree()
    if n==0: return dict(strictly_hurwitz=True,degree=0,determinants=[])
    coeff=poly.monic().all_coeffs()
    H=sp.zeros(n)
    for i in range(n):
        for j in range(n):
            k=2*j-i+1
            if 0<=k<=n: H[i,j]=coeff[k]
    determinants=[H[:k,:k].det() for k in range(1,n+1)]
    return dict(strictly_hurwitz=all(d>0 for d in determinants),degree=int(n),
                determinants=[str(d) for d in determinants])

def axis_polynomial(N,D):
    E=(N.as_expr()*D.as_expr().subs(s,-s)+N.as_expr().subs(s,-s)*D.as_expr())/2
    E=sp.Poly(sp.expand(E),s,domain=sp.QQ)
    return sp.Poly(sum(E.nth(2*j)*(-1)**j*x**j for j in range(E.degree()//2+1)),x,domain=sp.QQ) if not E.is_zero else sp.Poly(0,x,domain=sp.QQ)

def nonnegative_halfline(P):
    """Complete P>=0 on x>=0, allowing even positive roots.

    Remove the zero root. A positive initial sign and the absence of positive
    roots of odd multiplicity are necessary and sufficient. Square-free factors
    determine multiplicity; exact Sturm counts determine positive roots.
    """
    if P.is_zero: return dict(nonnegative=True,identically_zero=True,zero_multiplicity=None,factors=[])
    zero=0
    while P.nth(zero)==0: zero+=1
    reduced=sp.Poly(P.as_expr()/x**zero,x,domain=sp.QQ)
    initial_sign=sp.sign(reduced.nth(0))
    content,factors=sp.sqf_list(reduced)
    details=[]; odd_positive_roots=0
    for factor,multiplicity in factors:
        count=int(factor.count_roots(0,sp.oo))
        details.append(dict(factor=str(factor.as_expr()),multiplicity=int(multiplicity),positive_roots=count))
        if multiplicity%2: odd_positive_roots+=count
    return dict(nonnegative=bool(initial_sign>0 and odd_positive_roots==0),identically_zero=False,
                zero_multiplicity=zero,right_of_zero_sign=int(initial_sign),
                odd_multiplicity_positive_roots=odd_positive_roots,factors=details)

def certificate(Nasc,Dasc,reduce_common=False):
    """Return an exact positive-real decision for a real rational function.

    Default contract requires coprime polynomials. Set reduce_common=True to
    cancel the exact gcd first; cancellations are recorded in the certificate.
    """
    N,D=polynomial(Nasc),polynomial(Dasc)
    if D.is_zero: return dict(valid_input=False,positive_real=False,reason='zero denominator')
    gcd=sp.gcd(N,D); gcd_degree=int(gcd.degree())
    if gcd_degree>0 and not reduce_common:
        return dict(valid_input=False,positive_real=None,reason='coprime input required; use exact reduction',gcd=str(gcd.as_expr()))
    if gcd_degree>0:
        N=N.exquo(gcd); D=D.exquo(gcd)
    S=N+D; P=axis_polynomial(N,D)
    hurwitz=strict_hurwitz(S); nonnegative=nonnegative_halfline(P)
    return dict(valid_input=True,positive_real=bool(hurwitz['strictly_hurwitz'] and nonnegative['nonnegative']),
                numerator=str(N.as_expr()),denominator=str(D.as_expr()),exact_cancelled_degree=gcd_degree if reduce_common else 0,
                N_plus_D=str(S.as_expr()),hurwitz=hurwitz,axis_polynomial=str(P.as_expr()),axis=nonnegative)

def main():
    eps=Fraction(1,10**50)
    cases={
        'unstable_1_over_1_minus_s':([1],[1,-1],False),
        'lossless_foster_s_over_s2_plus_1':([0,1],[1,0,1],True),
        'capacitor_1_over_s':([1],[0,1],True),
        'inductor_s':([0,1],[1],True),
        'foster_interlacing':([1,0,1],[0,4,0,1],True),
        'lossy_first_order':([2,1],[1,1],True),
        'nonPR_s_minus_1_over_s_plus_1':([-1,1],[1,1],False),
        'rational_damping_one_tenth':([0,1],[1,Fraction(1,10),1],True),
        'nonPR_negative_axis_band':([1,0,1],[2,0,1],False),
        'even_positive_root_counterexample':([1,0,1],[1,1,1],True),
        'two_even_positive_roots':([4,0,5,0,1],[4,2,5,1,1],True),
        'zero_at_origin':([0,1],[1,1],True),
        'zero_function':([0],[1],True),
        'zero_function_negative_scale':([0],[-1],True),
        'constant_one':([1],[1],True),
        'constant_one_negative_scale':([-1],[-1],True),
        'negative_constant':([-1],[1],False),
        'negative_constant_negative_scale':([1],[-1],False),
        'negative_inductor':([0,-1],[1],False),
        'second_order_pole_at_zero':([1],[0,0,1],False),
        'wrong_lossless_residue':([0,-1],[1,0,1],False),
        'quadratic_at_infinity':([0,0,1],[1],False),
        'cubic_at_infinity':([0,0,0,1],[1],False),
        'near_boundary_positive_damping':([1,0,1],[1,eps,1],True),
        'near_boundary_negative_damping':([1,0,1],[1,-eps,1],False),
        'near_tangency_bad_band_above':([1+eps,0,1],[1,1,1],False),
        'near_tangency_bad_band_below':([1-eps,0,1],[1,1,1],False),
    }
    results={}
    for name,(N,D,expected) in cases.items():
        c=certificate(N,D)
        assert c['valid_input'] and c['positive_real']==expected, (name,c)
        # Scaling both exact polynomials must never affect the decision.
        scaled=certificate([-rational(q) for q in N],[-rational(q) for q in D])
        assert scaled['positive_real']==expected
        c['expected_positive_real']=expected; c['negative_scale_agrees']=True
        results[name]=c
    # Exact input/cancellation guards are separate from positive-real decisions.
    guards={
        'common_factor_strict':certificate([-1,1],[-1,1]),
        'common_factor_exactly_reduced':certificate([-1,1],[-1,1],reduce_common=True),
        'zero_over_nonconstant_strict':certificate([0],[1,1]),
        'zero_over_nonconstant_reduced':certificate([0],[1,1],reduce_common=True),
        'zero_denominator':certificate([1],[0]),
        'N_plus_D_zero':certificate([-1],[1]),
    }
    assert guards['common_factor_strict']['valid_input'] is False
    assert guards['common_factor_exactly_reduced']['positive_real'] is True
    assert guards['zero_over_nonconstant_reduced']['positive_real'] is True
    try:
        certificate([.1],[1])
        raise AssertionError('float input was accepted')
    except TypeError:
        guards['ambiguous_float_rejected']=True
    out=dict(method='exact rational Cayley/Brune, exact Hurwitz determinants, square-free multiplicities + Sturm counts',
             total_decision_cases=len(results),polarity_checks=len(results),all_expected_decisions_agree=True,cases=results,guards=guards)
    (HERE/'lossless_passivity_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(cases=len(results),polarity_checks=len(results),all_passed=True,
                         counterexample=results['even_positive_root_counterexample']),indent=2))

if __name__=='__main__': main()
