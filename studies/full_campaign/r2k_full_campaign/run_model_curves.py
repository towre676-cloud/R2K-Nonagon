from numerical_campaign import Evaluator,curve,ROOT
import json
archive=json.loads((ROOT/'part_VI/bk_multi_100000_1000000_10000000.json').read_text())
e=Evaluator(2000000)
for h in archive['heights']:
 for i,T in enumerate(h['T_quantiles']):curve(h['n0'],T,str(i),e)
print('All nine model curves completed.',flush=True)
