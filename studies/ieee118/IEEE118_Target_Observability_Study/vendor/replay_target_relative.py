#!/usr/bin/env python3
"""Fast replay for the target-relative IEEE-118 development receipt."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import target_relative_panel as T
import exact_target_rank as E
from fixtures import CASES
HERE=Path(__file__).resolve().parent


def main():
    suite=json.load(open(HERE/'target_relative_results.json'))
    x=suite['ieee118_adaptive']; exact=json.load(open(HERE/'target_relative_exact_rank.json'))
    buses,K,B,wtrue,keys=T.edge_parameterization(CASES['case118'])
    nominal,scale,t,h0=T.nominal_design(B,K); pats=T.patterns(len(buses),t)
    ytrue,_=T.evaluate(wtrue,B,pats)
    sel=x['local_kernel_panel']['selected_indices']; sel2=x['adaptive_panel']['selected_indices']; add=x['adaptive_branch_breaker']['index']
    if sel2 != sel+[add]: raise ValueError('adaptive panel is not minimal panel plus saved branch breaker')
    # Numerical local rank at the actual frozen nominal design.
    _,J=T.evaluate(nominal,B,pats,want_jacobian=True)
    norms=np.linalg.norm(J,axis=0); A=J/norms
    rank=int(np.linalg.matrix_rank(A[sel],tol=1e-10))
    if rank != len(keys): raise ValueError('saved local panel lost full rank')
    # Saved exact-fit branches must really agree on all 179 measurements.
    S=x['saved_branch_weights_for_replay']
    wf=np.asarray(S['minimal_false_branch']); wt=np.asarray(S['minimal_true_branch'])
    yf,_=T.evaluate(wf,B,pats,sel); yt,_=T.evaluate(wt,B,pats,sel)
    scale_y=max(abs(ytrue[sel][0]),float(np.ptp(ytrue[sel])),1e-12)
    rf=float(np.linalg.norm((yf-ytrue[sel])/scale_y)); rt=float(np.linalg.norm((yt-ytrue[sel])/scale_y))
    if rf>1e-9 or rt>1e-9: raise ValueError('saved minimal branches do not replay as exact fits')
    if np.linalg.norm(wf-wt)/np.linalg.norm(wtrue) < 1e-5: raise ValueError('saved branches are not distinct')
    # The adaptive reading must split the branches and the final saved candidates must fit it.
    yfa,_=T.evaluate(wf,B,pats,[add]); yta,_=T.evaluate(wt,B,pats,[add]); yadd=float(ytrue[add])
    false_offset=float(yfa[0]-yadd); true_offset=float(yta[0]-yadd)
    if abs(false_offset)<1e-8 or abs(true_offset)>1e-10: raise ValueError('saved branch-breaker no longer distinguishes the branches')
    finals={k:np.asarray(S[k]) for k in ('adaptive_from_false','adaptive_from_true','adaptive_independent')}
    final_errors={}
    for name,w in finals.items():
        y,_=T.evaluate(w,B,pats,sel2); residual=float(np.linalg.norm((y-ytrue[sel2])/max(abs(ytrue[sel2][0]),float(np.ptp(ytrue[sel2])),1e-12)))
        rel=float(np.linalg.norm(w-wtrue)/np.linalg.norm(wtrue))
        final_errors[name]={'residual_norm':residual,'relative_weight_error':rel}
        if residual>1e-9 or rel>1e-8: raise ValueError(f'final saved solution failed replay: {name}')
    # Recompute the finite-field determinant witnesses rather than trusting JSON.
    Bi=np.rint(B).astype(np.int64); Kunit=Bi@Bi.T; punit=T.patterns(len(buses),1.0); mod=[]
    for cert in exact['certificates']:
        p=int(cert['prime']); G=E.inv_mod_matrix(Kunit,p); d=E.certificate_for_prime(Bi,Kunit,G,punit,sel,p)
        mod.append({'prime':p,'determinant_mod_prime':int(d)})
        if d==0 or d!=int(cert['determinant_mod_prime']): raise ValueError('finite-field rank certificate mismatch')
    receipt={
      'schema':'nonagon.target-relative-replay.v1','status':'PASS','case':'case118',
      'parameter_count':len(keys),'full_protocol_count':len(pats),'local_panel_count':len(sel),'adaptive_panel_count':len(sel2),
      'local_rank':rank,'minimal_false_branch_residual':rf,'minimal_true_branch_residual':rt,
      'branch_breaker_false_offset':false_offset,'branch_breaker_true_offset':true_offset,
      'final_saved_solutions':final_errors,'exact_rank_modular_certificates':mod,
    }
    (HERE/'target_relative_replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()
