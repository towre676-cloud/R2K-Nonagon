"""JSON-friendly research adapter; fixed log acquisition remains the default."""
import numpy as np
from gate import choose_strategy,adaptive

def decide(payload):
    mode=payload.get('mode','fixed_log')
    if mode=='fixed_log':
        return {'strategy':'log_uniform','experimental':False}
    if mode!='experimental_gate':raise ValueError('unknown policy mode')
    omega=np.asarray(payload['omega'],dtype=float)
    values=np.asarray(payload['response_real'])+1j*np.asarray(payload['response_imag'])
    if omega.shape!=(33,) or values.shape!=(33,):raise ValueError('requires exactly33commonpilotvalues')
    if not np.all(np.isfinite(omega)) or not np.all(np.isfinite(values)):raise ValueError('nonfinitepilot')
    order=int(payload['order']);band=tuple(payload.get('band',[.03,8.]))
    if order not in (6,8,10) or band!=(.03,8.):raise ValueError('outside tested order/band')
    if not np.allclose(omega,adaptive.POOL_W[adaptive.PILOT],rtol=1e-12,atol=0):raise ValueError('requires the tested common pilot')
    strategy,features=choose_strategy(omega,values,order,band)
    return {'strategy':strategy,'features':features,'experimental':True,
            'limitation':'Gate did not confirm its primary mean-risk advantage; do not treat as tail protection.'}
