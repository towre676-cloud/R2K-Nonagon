import json, sys, glob, time, r2kverify
mode = sys.argv[1]; pat = sys.argv[2] if len(sys.argv) > 2 else 'certs/*.json'
out = {}
try: out = json.load(open('verify_%s.json' % mode))
except Exception: pass
for fn in sorted(glob.glob(pat)):
    if 'log' in fn or fn in out: continue
    c = json.load(open(fn)); t0 = time.time()
    v = r2kverify.verify(c, full=(mode == 'full')); v['time'] = time.time() - t0
    out[fn] = v; print(fn, v['status'], v.get('count'), v.get('inertia'), v.get('reason', ''), '%.1fs' % v['time'], flush=True)
    json.dump(out, open('verify_%s.json' % mode, 'w'), indent=1)
