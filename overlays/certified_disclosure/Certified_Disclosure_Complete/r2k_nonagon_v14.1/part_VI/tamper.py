import json, copy, hashlib, r2kverify, certgen
base = json.load(open('certs/zero_1_b0p45.json'))
def rehash(c):
    c.pop('sha256', None); c['sha256'] = hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest(); return c
tests = []
def T(name, c, full=False):
    v = r2kverify.verify(c, full=full); tests.append({'test': name, 'mode': 'full' if full else 'light', 'status': v['status'], 'reason': v.get('reason', '')})
    print(tests[-1], flush=True)
T('untouched', copy.deepcopy(base)); T('untouched', copy.deepcopy(base), True)
c = copy.deepcopy(base); c['segments'][3]['contrib'][0][0] = '0.1'; T('ball altered, hash not updated', c)
c = copy.deepcopy(base); s = c['segments'][3]['contrib'][0]; s[2] = str(float(s[2]) + 3.14159); rehash(c); T('mu0 contribution shifted by pi*i (+1 zero), rehashed', c)
T('mu0 contribution shifted by pi*i (+1 zero), rehashed', copy.deepcopy(c), True)
c = copy.deepcopy(base); s = c['segments'][3]['contrib'][1]; s[0] = str(float(s[0]) + 1e-6); rehash(c); T('mu1 contribution perturbed 1e-6, rehashed', c, True)
c = copy.deepcopy(base); del c['segments'][5]; rehash(c); T('segment deleted, rehashed', c)
c = copy.deepcopy(base); sg = c['segments'][0]; sg['D'] = [-sg['D'][0], -sg['D'][1]]; rehash(c); T('side orientation reversed, rehashed', c)
c = copy.deepcopy(base); sg = c['segments'][2]; sg['x1'] = (float.fromhex(sg['x1']) + 2**-20).hex(); rehash(c); T('segments overlap, rehashed', c)
c = copy.deepcopy(base); c['rectangle']['b'] = (0.5).hex(); rehash(c); T('rectangle claim widened to b=0.5, rehashed', c)
c = copy.deepcopy(base); c['segments'][4]['R2'] = (1e-5).hex(); rehash(c); T('disk radius shrunk (segment too long for disk), rehashed', c, True)
c = copy.deepcopy(base); c['method']['K'] = 8; rehash(c); T('Taylor order lowered to K=8, rehashed', c, True)
leg = json.load(open('certs/legacy_v5/cert_tile1.json')); T('legacy v5 record (ctx.cap bug), light', leg); T('legacy v5 record (ctx.cap bug), full replay', leg, True)
try:
    certgen.generate('xi', 10.0, 14.134725141734693, 0.45, [0, 1], True, mode='series'); tests.append({'test': 'generator on contour through first zero', 'status': 'PRODUCED (BAD)'})
except Exception as e: tests.append({'test': 'generator on contour through first zero', 'status': 'REFUSED', 'reason': str(e)[:160]})
print(tests[-1])
json.dump(tests, open('tamper_tests.json', 'w'), indent=1)
