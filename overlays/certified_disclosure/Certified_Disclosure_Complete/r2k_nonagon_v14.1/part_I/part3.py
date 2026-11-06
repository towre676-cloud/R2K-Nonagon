import numpy as np, json
R=json.load(open('results.json'))
t=R['E']['tiles']; print('tiles max err',max(x['err'] for x in t),'max mu0dev',max(x['mu0_dev'] for x in t),'max cons',max(x['cheb_consistency'] for x in t), len(t))
zeta_d,tau=0.1,1.0
fD=lambda s,k: s**2+2*zeta_d*s+1+k*np.exp(-s*tau); dfD=lambda s,k: 2*s+2*zeta_d-k*tau*np.exp(-s*tau)
gx,gw=np.polynomial.legendre.leggauss(32)
def count(k,Y=40.0):
    cs=[0-1j*Y,6-1j*Y,6+1j*Y,0+1j*Y,0-1j*Y]; tot=0
    for side,(A,B) in enumerate(zip(cs[:-1],cs[1:])):
        pn=40000 if side==3 else 400
        a=A+(B-A)*np.arange(pn)/pn; b=A+(B-A)*np.arange(1,pn+1)/pn
        s=((a+b)/2)[:,None]+((b-a)/2)[:,None]*gx[None,:]
        tot+=np.sum(dfD(s,k)/fD(s,k)*gw[None,:]*((b-a)/2)[:,None])
    return tot/(2j*np.pi)
kc=R['D']['kappa_crit_analytic']
grid=np.linspace(0,1.2,121); cnt=np.array([count(k) for k in grid])
R['D']['kappa_grid']=grid.tolist(); R['D']['rhp_counts']=cnt.real.tolist()
mask=np.abs(grid-kc)>0.004
R['D']['count_integrality_dev']=float(np.max(np.abs(cnt[mask]-np.round(cnt[mask].real))))
print(R['D']['count_integrality_dev'], [ (round(k,2), round(c,9)) for k,c in zip(grid,cnt.real)][20:30])
json.dump(R,open('results.json','w'),indent=1)
