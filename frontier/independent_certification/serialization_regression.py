"""Human-readable Arb balls need not preserve a visible nonzero margin."""
import json
from pathlib import Path
from direct_engine import arb,acb,ctx,serialize,deserialize_exact,require

ctx.prec=100
cases=[]
for mid,rad in [('1','0.9'),('-2','1.9'),('1/100','9/1000')]:
    def q(s):
        if '/' in s:
            n,d=s.split('/');return arb(int(n))/int(d)
        return arb(s)
    value=acb(arb(q(mid),q(rad)),arb(0))
    saved=serialize(value);back=deserialize_exact(saved)
    display=acb(arb(saved['real']),arb(saved['imag']))
    require(not value.contains(0) and not back.contains(0),'exact dyadic roundtrip lost zero exclusion')
    require(back.contains(value),'dyadic roundtrip does not enclose original')
    cases.append({'display':saved['real'],'display_contains_zero':bool(display.contains(0)),
                  'original_excludes_zero':True,'exact_roundtrip_excludes_zero':True})
require(any(r['display_contains_zero'] for r in cases),'fixture did not exercise display coarsening')
result={'passed':len(cases),'cases':cases,'interpretation':'display strings are diagnostics; exact dyadic fields preserve saved boundary balls'}
(Path(__file__).resolve().parent/'serialization_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
