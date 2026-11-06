"""Summarize final completed campaign, replay, controls and retained limitation."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from direct_engine import require,sha

HERE=Path(__file__).resolve().parent


def load(name):return json.loads((HERE/name).read_text())


def main():
    current=load('results.json');campaign=load('campaign_results.json');replay=load('replay_results.json')
    controls=load('controls_results.json');serialization=load('serialization_results.json');p=load('protocol.json')
    require(hashlib.sha256((HERE/'direct_engine.py').read_bytes()).hexdigest()==p['engine_sha256'],'source changed after campaign')
    require(current['protocol_sha256']==campaign['protocol_sha256']==sha(p),'protocol hashes disagree')
    require(replay['fresh_replay']['all_agree'],'fresh replay did not complete')
    prior=load('initial_display_run/lehmer_results.json')
    prior['engine_sha256']=load('initial_display_run/protocol.json')['engine_sha256']
    prior['source_difference_from_current']='interval serialization supplement only; distinct source hash retained'
    result={'protocol_sha256':sha(p),'engine_sha256':p['engine_sha256'],
            'campaign':{'continuous_t_range':p['range'],'s_real_range':['-1/10','11/10'],
                'total_zero_count':campaign['total_zero_count'],'all_real_simple':campaign['all_regions_certified_real_simple_or_empty'],
                'attempts':len(campaign['attempts']),'leaf_tiles':len(campaign['leaves']),
                'count_histogram':dict(Counter(r['count'] for r in campaign['leaves'])),
                'unresolved_leaves':len(campaign['unresolved_leaves']),'seconds':campaign['seconds'],
                'boundary_value_calls':sum(r['boundary_value_calls'] for r in campaign['leaves']),
                'logarithmic_derivative_calls':sum(r['logarithmic_derivative_calls'] for r in campaign['leaves'])},
            'fresh_replay':{'all_agree':True,'tiles':len(replay['fresh_replay']['tiles']),'seconds':replay['fresh_replay']['seconds']},
            'analytic_controls_passed':controls['passed'],'integrity_tamper_rejections':len(replay['tamper_checks']),
            'interval_serialization_controls_passed':serialization['passed'],
            'post_campaign_reference':load('reference_check.json'),
            'retained_measured_limitations':{'direct_lehmer_attempt':prior},
            'trust_boundary':'independent direct contour algorithm; Arb function balls and adaptive integration trusted; no archive Taylor/DFT kernel',
            'scope':'reproduction of known finite-height zero facts over exactly[10,100]; no new height record or global RH claim'}
    (HERE/'final_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
