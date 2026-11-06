from fractions import Fraction as F
import random

import pytest

from two_strength_cert import (Interval, certify_pair, pair_delta, reconstruct_single,
                               verify_certificate, certify_many, verify_many_certificate,
                               certify_raw_pair, verify_raw_certificate)


def test_interval_arithmetic_and_sqrt_enclose_rationals():
    a=Interval(F(-2,3),F(3,4));b=Interval(F(1,5),F(2,5))
    for x in (a.lo,a.mid,a.hi):
        for y in (b.lo,b.mid,b.hi):
            assert (a+b).contains(x+y)
            assert (a*b).contains(x*y)
            assert (a/b).contains(x/y)
    r=Interval(F(2),F(3)).sqrt()
    assert r.lo*r.lo<=2 and r.hi*r.hi>=3


def test_exact_common_root_replay_without_square_root_selection():
    rng=random.Random(314)
    for _ in range(24):
        # Nonsingular local response blocks, including indefinite examples.
        u,v,A,B,x=[F(rng.randint(1,25),11) for k in range(5)]
        u*=(-1 if rng.random()<.5 else 1)
        v*=(-1 if rng.random()<.5 else 1)
        x*=(-1 if rng.random()<.5 else 1)
        t1,t2=F(3,7),F(9,7)
        d1=pair_delta(t1,u,v,A,B,x);d2=pair_delta(t2,u,v,A,B,x)
        if d1==d2:continue
        cert=certify_pair(t1,t2,u,v,A,B,d1,d2)
        assert list(map(F,cert['enclosure']))==[x,x]
        assert verify_certificate(cert)


def test_nonzero_budgets_enclose_all_sampled_consistent_roots():
    rng=random.Random(42)
    for _ in range(25):
        u,v,A,B,x=F(4,5),F(3,4),F(5,2),F(7,3),F(rng.randint(-8,8),20)
        t1,t2=F(1,2),F(2)
        d1,d2=pair_delta(t1,u,v,A,B,x),pair_delta(t2,u,v,A,B,x)
        eps=F(1,10**8)
        cert=certify_pair(t1,t2,Interval.ball(u,eps),Interval.ball(v,eps),
                          Interval.ball(A,eps),Interval.ball(B,eps),
                          Interval.ball(d1+eps/3,eps),Interval.ball(d2-eps/2,eps))
        assert Interval(*map(F,cert['enclosure'])).contains(x)
        assert verify_certificate(cert)


def test_denominator_gate_and_tamper_rejection():
    with pytest.raises(ZeroDivisionError):
        certify_pair(1,2,1,1,2,2,Interval(0,2),Interval(1,3))
    with pytest.raises(ZeroDivisionError):
        certify_pair(1,2,Interval(-1,1),1,2,2,1,2)
    u,v,A,B,x=map(F,(1,2,3,4,1))
    cert=certify_pair(1,2,u,v,A,B,pair_delta(1,u,v,A,B,x),pair_delta(2,u,v,A,B,x))
    cert['enclosure']=['0','0']
    assert not verify_certificate(cert)


def test_single_node_calibration_propagates_raw_reading_error():
    theta=F(2,3);u=F(4,7);A=F(5,4);eps=F(1,10**9)
    d1=theta*u*u/(1+theta*A);d2=2*theta*u*u/(1+2*theta*A)
    AA,uu=reconstruct_single(theta,Interval.ball(d1,eps),Interval.ball(d2,eps))
    assert AA.contains(A) and uu.contains(u)
    Aneg,uneg=reconstruct_single(theta,Interval.ball(d1,eps),Interval.ball(d2,eps),sign=-1)
    assert Aneg.contains(A) and uneg.contains(-u)


def test_single_node_degeneracy_fails_explicitly():
    with pytest.raises(ZeroDivisionError):
        reconstruct_single(1,Interval(0,1),Interval(0,1))
    with pytest.raises(ValueError):
        reconstruct_single(0,1,2)


def test_one_zero_pair_difference_is_allowed():
    # Q1 is linear here; elimination still works without dividing by delta1.
    u,v,A,B,x=map(F,(1,2,1,3,3))
    d1=pair_delta(1,u,v,A,B,x);d2=pair_delta(2,u,v,A,B,x)
    assert d1==0 and d2!=0
    cert=certify_pair(1,2,u,v,A,B,d1,d2)
    assert cert['enclosure']==['3','3']
    assert verify_certificate(cert)


def test_multiple_strengths_fuse_without_truth_and_never_widen():
    u,v,A,B,x=map(F,(1,2,3,4,1));theta=[F(1,10),F(1),F(10),F(100)]
    ds=[Interval.ball(pair_delta(t,u,v,A,B,x),F(1,10**8)) for t in theta]
    cert=certify_many(theta,ds,u,v,A,B)
    X=Interval(*map(F,cert['enclosure']))
    assert X.contains(x) and verify_many_certificate(cert)
    for i in range(len(theta)):
        for j in range(i+1,len(theta)):
            Y=certify_pair(theta[i],theta[j],u,v,A,B,ds[i],ds[j],raw_interval=True)
            assert X.width<=Y.width
    cert['enclosure']=['0','0']
    assert not verify_many_certificate(cert)


def test_raw_algebra_cancellation_preserves_coverage_and_tightens():
    theta=F(2,3);u,v,A,B,x=F(4,7),F(5,8),F(5,4),F(7,5),F(3,10)
    eps=F(1,10**8)
    def s(t,col,diag):return Interval.ball(t*col*col/(1+t*diag),eps)
    si1,si2,sj1,sj2=s(theta,u,A),s(2*theta,u,A),s(theta,v,B),s(2*theta,v,B)
    d1=Interval.ball(pair_delta(theta,u,v,A,B,x),eps)
    d2=Interval.ball(pair_delta(2*theta,u,v,A,B,x),eps)
    cert=certify_raw_pair(theta,si1,si2,sj1,sj2,d1,d2)
    X=Interval(*map(F,cert['enclosure']));Y=Interval(*map(F,cert['propagated']))
    assert X.contains(x) and X.width<=Y.width
    assert verify_raw_certificate(cert)
