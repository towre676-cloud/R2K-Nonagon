"""Replay stored enclosure/tail and positivity certificates, without FFT or eigensolver."""
from pathlib import Path
import json
from flint import arb,arb_mat,ctx
HERE=Path(__file__).resolve().parent
ctx.prec=192

def main():
    p=json.loads((HERE/'cue_covariance_results.json').read_text());K=p['cutoff'];checks=0
    assert p['status']=='complete' and K>max(map(int,p['models']))
    margins=[]
    for nn,m in p['models'].items():
        N=int(nn);tail=4*arb(N)**2/(3*arb.pi()**4*K**3)*(1+4/(1-arb(N)/K)**2)
        A=[]
        for i,row in enumerate(m['covariance']):
            aa=[]
            for j,c in enumerate(row):
                finite=arb(c['finite']['ball']);whole=arb(c['enclosure']['ball'])
                assert whole.contains(finite+arb(0,tail.upper()))
                assert finite.overlaps(arb(c['Gaussian_pairing_finite']['ball'])+arb(c['fourth_cumulant_finite']['ball']))
                assert c==m['covariance'][j][i]
                aa.append(whole);checks+=1
            A.append(aa)
        cert=m['positivity_certificate']
        P=[[arb(float.fromhex(x)) for x in row] for row in cert['preconditioner_hex']]
        det=arb_mat(P).det();assert abs(det)>0
        B=[[sum((P[i][k]*A[k][l]*P[j][l] for k in range(6) for l in range(6)),arb(0)) for j in range(6)] for i in range(6)]
        gs=[B[i][i]-sum((abs(B[i][j]) for j in range(6) if i!=j),arb(0)) for i in range(6)]
        assert all(g>0 for g in gs)
        margins.append({'rank':N,'replayed_Gershgorin_margins':[str(x) for x in gs]})
    out={'status':'PASS','entry_tail_enclosure_checks':checks,'invertible_congruence_certificates':len(margins),'positive_Gershgorin_margins':margins}
    (HERE/'covariance_validation.json').write_text(json.dumps(out,indent=2));print('PASS',checks,'entry enclosures; three positivity certificates')

if __name__=='__main__':main()
