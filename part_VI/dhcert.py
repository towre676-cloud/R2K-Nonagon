import json, certgen, r2kverify
c = certgen.generate('dh', 84.7, 86.7, 0.45, [0, 1, 2], False, mode='series', R2max=0.05)
json.dump(c, open('certs/NEG_dh_84p7_86p7.json', 'w'), indent=1)
print('segments', len(c['segments']), c['generator_moments'][:3])
for full in (False, True): print(json.dumps(r2kverify.verify(c, full=full)))
