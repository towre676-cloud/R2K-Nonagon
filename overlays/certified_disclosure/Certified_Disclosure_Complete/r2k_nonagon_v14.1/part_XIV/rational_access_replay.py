"""Replay exact rational access certificates without solving a rational nullspace.

Every inverse is rebuilt from the MATPOWER rational inputs. Every retained row,
rank minor, kernel generator and coordinate row-span witness is checked exactly.
"""
import argparse, hashlib, json, time
from pathlib import Path
import flint
from rational_identifiability import HERE, P, load_exact, exact_rows, hash_matrix, modular_counterexample


def verify_certificate(D, cert):
    if not __debug__:
        raise RuntimeError("Exact certificate verification requires assertions; run without python -O")
    ports=[D['nodes'].index(b) for b in cert['ports']]
    where=[D['nodes'].index(b) for b in cert['intervention_buses']]
    U=D['U'];r=cert['rational_rank'];N=cert['kernel_generators'];free=cert['kernel_free_coordinates']
    assert cert['case']==D['case']
    assert cert['rank_witness_prime']==P and cert['unknowns']==U and cert['branches']==len(D['el'])
    assert isinstance(r,int) and 0<=r<=U
    assert r+len(N)==U and cert['rational_nullity']==U-r and len(free)==len(N)
    assert len(set(free))==len(free) and all(isinstance(f,int) and 0<=f<U for f in free)
    assert all(len(v)==U and all(type(a) is int for a in v) for v in N)
    assert all(N[k][f]!=0 and all(N[l][f]==0 for l in range(len(N)) if l!=k) for k,f in enumerate(free))
    rowids=cert['rank_witness_rows'];cols=cert['rank_witness_columns']
    assert len(rowids)==r and len(set(rowids))==r
    assert all(isinstance(i,int) and 0<=i<cert['readings'] for i in rowids)
    assert len(cols)==r and len(set(cols))==r and all(isinstance(j,int) and 0<=j<U for j in cols)
    rowpos={idx:k for k,idx in enumerate(rowids)}
    actualrows=[None]*r
    sparse=[[(i,a) for i,a in enumerate(v) if a] for v in N]
    h=hashlib.sha256();m=checks=matched=0
    for idx,row in enumerate(exact_rows(D,ports,where,cert['driving_point_only_after_intervention'])):
        if idx in rowpos:
            actualrows[rowpos[idx]]=row
            matched+=1
        h.update((','.join(map(str,row))+'\n').encode())
        for v in sparse:
            assert sum(row[i]*a for i,a in v)==0
            checks+=1
        m+=1
    assert m==cert['readings'] and matched==r
    assert h.hexdigest()==cert['all_integer_jacobian_rows_sha256']
    assert checks==cert['exact_kernel_dot_products_checked']
    assert hash_matrix(actualrows)==cert['selected_integer_rows_sha256']
    if 'selected_integer_rows' in cert:assert cert['selected_integer_rows']==actualrows
    minor=flint.nmod_mat([[actualrows[i][j]%P for j in cols] for i in range(r)],P)
    assert int(minor.det())==cert['rank_witness_minor_determinant_mod_p']!=0
    ident=[i for i in range(U) if all(v[i]==0 for v in N)]
    assert ident==cert['identifiable_coordinate_indices']
    assert cert['first_order_identifiable_unknowns']==len(ident)
    assert cert['first_order_identifiable_branches']==sum(i<len(D['el']) for i in ident)
    seen=[]
    for wit in cert['coordinate_row_span_witnesses']:
        i=wit['coordinate'];lam=flint.fmpq_mat([[flint.fmpq(x) for x in wit['coefficients']]])
        assert lam.ncols()==r
        v=lam*flint.fmpq_mat(actualrows)
        assert all(v[0,j]==int(i==j) for j in range(U))
        seen.append(i)
    assert seen==ident or (r==U and cert.get('full_rank_proves_all_coordinates_without_individual_row_span_witnesses') and not seen) or (cert.get('coordinate_classification_via_complete_kernel') is True and not seen)
    return {'case':cert['case'],'scheme':cert['scheme'],'readings':m,'rational_rank':r,'rational_nullity':len(N),
            'exact_kernel_products':checks,'coordinate_row_span_witnesses':len(seen),'verified':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('certificates',nargs='*');ap.add_argument('--receipt',type=Path,default=HERE/'rational_access'/'replay_receipt.json');a=ap.parse_args()
    paths=[Path(x) for x in a.certificates] or sorted((HERE/'rational_access').glob('rational_*.json'))
    models={};out=[];t=time.time()
    for p in paths:
        c=json.loads(p.read_text());case=c['case']
        if case not in models:models[case]=load_exact(case)
        r=verify_certificate(models[case],c);r['certificate_sha256']=hashlib.sha256(p.read_bytes()).hexdigest();out.append(r);print(json.dumps(r),flush=True)
    result={'all_passed':True,'certificates':out,'counterexample':modular_counterexample(),'seconds':time.time()-t}
    a.receipt.write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()
