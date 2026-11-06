"""Meaningful algebraic regressions for exact access certificates."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import flint
from rational_identifiability import P, rational_kernel, modular_counterexample, certify, exact_columns
from rational_access_replay import verify_certificate


def test_modular_coordinate_counterexample():
    c=modular_counterexample()
    assert c['mod_p_apparently_identifiable_coordinates']==[0]
    assert c['rational_first_order_identifiable_coordinates']==[]
    R,r,piv,free,N=rational_kernel(flint.fmpq_mat([[1,P]]))
    assert r==1 and N==[[-P,1]]


def toy():
    # Two hidden buses, one branch and a ground at each bus.
    return {'case':'two_bus','n':2,'nodes':[10,20],'el':[(0,1)],'gs':[0,1],
            'port':0,'adj':{0:{1},1:{0}},'hw':[0], 'M':[[2,1],[1,3]],'d':5,'U':3}


def test_denominator_cleared_intervention_vectors():
    D=toy();K=flint.fmpq_mat([[3,-1],[-1,2]])
    for j in (0,1):
        Kj=flint.fmpq_mat(K);Kj[j,j]+=1;G=Kj.inv()
        X=exact_columns(D,[0,1],j)
        for p in (0,1):
            assert flint.fmpq(X[p][0])*G[1,p]==flint.fmpq(X[p][1])*G[0,p]


def test_exact_partial_access_certificate_and_replay():
    D=toy();c=certify(D,[0],[],include_witnesses=True);c['scheme']='one_port_no_intervention'
    assert c['rational_rank']==1 and c['rational_nullity']==2
    assert c['first_order_identifiable_unknowns']==0
    assert verify_certificate(D,c)['verified']


def test_full_access_certificate_and_replay():
    D=toy();c=certify(D,[0,1],[],include_witnesses=True);c['scheme']='two_ports'
    assert c['rational_rank']==3 and c['rational_nullity']==0
    assert c['first_order_identifiable_unknowns']==3
    assert verify_certificate(D,c)['verified']


def test_replay_rejects_false_kernel_witness():
    import pytest
    D=toy();c=certify(D,[0],[],include_witnesses=True);c['scheme']='corrupted'
    c['kernel_generators'][0][0]+=1
    with pytest.raises(AssertionError):verify_certificate(D,c)


def test_replay_rejects_unmatched_rank_rows():
    import pytest
    D=toy();c=certify(D,[0,1],[],include_witnesses=True);c['scheme']='tampered_rank'
    c['rank_witness_rows']=[c['readings']+i for i in range(c['rational_rank'])]
    with pytest.raises(AssertionError):verify_certificate(D,c)


def test_deficient_certificate_with_explicit_observable_coordinate():
    D={'case':'disconnected_two_bus','n':2,'nodes':[10,20],'el':[],'gs':[0,1],
       'port':0,'adj':{0:set(),1:set()},'hw':[],'M':[[3,0],[0,2]],'d':6,'U':2}
    c=certify(D,[0],[],include_witnesses=True);c['scheme']='one_observable_ground'
    assert c['rational_rank']==1 and c['identifiable_coordinate_indices']==[0]
    assert len(c['coordinate_row_span_witnesses'])==1
    assert verify_certificate(D,c)['verified']


def test_complete_kernel_duality_also_certifies_coordinates():
    D={'case':'disconnected_two_bus','n':2,'nodes':[10,20],'el':[],'gs':[0,1],
       'port':0,'adj':{0:set(),1:set()},'hw':[],'M':[[3,0],[0,2]],'d':6,'U':2}
    c=certify(D,[0],[],include_witnesses=False);c['scheme']='kernel_duality'
    assert c['identifiable_coordinate_indices']==[0] and c['coordinate_row_span_witnesses']==[]
    assert verify_certificate(D,c)['verified']
