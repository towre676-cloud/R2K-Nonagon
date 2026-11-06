#!/usr/bin/env python3
"""Replay the integrated Nonagon release from any working directory.

Python 3.11+ is required. The numerical theory and load-planning stages also use
NumPy/SciPy. No instrument, simulator binary, network request, or physical
hardware access is performed. Retained engineering simulation packets are
checked independently; this command does not rerun ngspice.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'cross_checks/integration_receipt.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'cross_checks/replay_outputs')
    parser.add_argument('--skip-numeric-audit', action='store_true',
                        help='Skip only the optional earlier-chapter endpoint audit; other stages still run')
    parser.add_argument('--timeout', type=int, default=1800,
                        help='Maximum seconds for each individual replay stage')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    logs = out / 'logs'
    logs.mkdir(exist_ok=True)
    started = time.monotonic()
    stages = []
    receipt = {'schema': 'nonagon.integration.release-replay.v1', 'status': 'running',
               'release': '9.0', 'python': platform.python_version(),
               'physical_hardware_measurements': 0,
               'scope': 'Exact certificates, independent cross-implementation recovery, retained engineering simulation packets, and fresh integration extensions',
               'stages': stages}

    def run(name, relative_script, arguments=(), optimized=False, expected=None):
        script = ROOT / relative_script
        require(script.is_file(), 'Required replay stage is missing: ' + relative_script)
        if expected:
            promised = out / expected
            if promised.exists():
                promised.unlink()
        argv = [sys.executable] + (['-O'] if optimized else []) + [str(script)] + [str(v) for v in arguments]
        begin = time.monotonic()
        stdout, stderr = logs / (name + '.stdout.txt'), logs / (name + '.stderr.txt')
        print(json.dumps({'stage': name, 'status': 'running'}), flush=True)
        with stdout.open('w', encoding='utf-8') as so, stderr.open('w', encoding='utf-8') as se:
            process = subprocess.run(argv, cwd=ROOT, stdout=so, stderr=se,
                                     text=True, timeout=args.timeout)
        row = {'stage': name, 'script': relative_script, 'optimized': optimized,
               'source_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
               'returncode': process.returncode,
               'elapsed_seconds': round(time.monotonic() - begin, 3),
               'stdout_sha256': hashlib.sha256(stdout.read_bytes()).hexdigest(),
               'stderr_sha256': hashlib.sha256(stderr.read_bytes()).hexdigest()}
        stages.append(row)
        require(process.returncode == 0, name + ' failed; inspect ' + str(stderr))
        if expected:
            target = out / expected
            require(target.is_file(), name + ' did not produce its promised receipt')
            record = read(target)
            require(str(record.get('status', '')).lower() in ('passed', 'certified'),
                    name + ' receipt does not report success')
            row['receipt_sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
        row['status'] = 'passed'
        print(json.dumps({'stage': name, 'status': 'passed', 'elapsed_seconds': row['elapsed_seconds']}), flush=True)

    try:
        run('provenance', 'cross_checks/audit_provenance.py',
            ['--output', out / 'provenance.json'], expected='provenance.json')
        common = ['--output', out / 'production.json']
        if args.skip_numeric_audit:
            common += ['--skip-numeric-audit']
        run('production', 'part_IX/replay.py', common, expected='production.json')
        run('production_optimized', 'part_IX/replay.py',
            ['--skip-numeric-audit', '--output', out / 'production_optimized.json'],
            optimized=True, expected='production_optimized.json')
        run('engineering_independent', 'part_IX/engineering_tests/verify_engineering.py',
            ['--output', out / 'engineering.json'], expected='engineering.json')
        for optimized in (False, True):
            suffix = '_optimized' if optimized else ''
            run('cross_schema' + suffix, 'cross_checks/bridge_claude_observations.py',
                ['--output', out / ('cross_schema' + suffix + '.json'),
                 '--records', out / ('bridge_records' + suffix + '.jsonl')],
                optimized=optimized, expected='cross_schema' + suffix + '.json')
        normal, optimized = read(out / 'cross_schema.json'), read(out / 'cross_schema_optimized.json')
        require(normal['records_sha256'] == optimized['records_sha256'],
                'Normal and optimized cross-schema scalar records differ')
        run('reference_input_controls', 'cross_checks/verifier_controls/claude_verifier_controls.py')
        controls = read(ROOT / 'cross_checks/verifier_controls/claude_verifier_controls.json')
        require(controls['cases'] == 9 and controls['runs'] == 18,
                'Expected all 18 original-oracle interface controls')
        run('theory', 'part_VIII/replay_theory.py',
            ['--output', out / 'theory.json'], expected='theory.json')
        for optimized in (False, True):
            suffix = '_optimized' if optimized else ''
            run('finite_objective_controls' + suffix, 'part_IX/load_planning/independent_checks.py',
                ['--out', out / ('finite_objective_controls' + suffix + '.json')],
                optimized=optimized, expected='finite_objective_controls' + suffix + '.json')
            run('finite_planner_input_controls' + suffix, 'part_IX/load_planning/planner_input_controls.py',
                ['--out', out / ('finite_planner_input_controls' + suffix + '.json')],
                optimized=optimized, expected='finite_planner_input_controls' + suffix + '.json')
        run('finite_load_campaign', 'part_IX/load_planning/benchmark.py',
            ['--output-dir', out / 'load_planning'],
            expected='load_planning/benchmark_receipt.json')
        run('finite_load_independent', 'part_IX/load_planning/verify_benchmark.py',
            ['--records', out / 'load_planning/benchmark_cases.jsonl',
             '--receipt', out / 'load_planning/benchmark_receipt.json',
             '--output', out / 'load_planning/independent_verification.json'],
            expected='load_planning/independent_verification.json')
        run('dynamic_acquisition', 'part_IX/dynamic_acquisition/replay.py',
            ['--output-dir', out / 'dynamic_acquisition'],
            expected='dynamic_acquisition/dynamic_receipt.json')
        run('dynamic_independent', 'part_IX/dynamic_acquisition/verify_records.py',
            [out / 'dynamic_acquisition', '--output', out / 'dynamic_acquisition/independent_receipt.json'],
            expected='dynamic_acquisition/independent_receipt.json')
        receipt['summary'] = {
            'exact_certificates': normal['certificates'],
            'unchanged_scalar_observations': normal['scalar_observations'],
            'heldout_observations': normal['heldout_observations'],
            'strict_certificate_checks': normal['strict_certificate_checks'],
            'independent_cross_schema_checks': normal['independent_operator_and_response_checks'],
            'normal_and_optimized_bridge_records_identical': True,
            'independent_engineering_predicates': read(out / 'engineering.json')['checks'],
            'preserved_input_files': read(out / 'provenance.json')['preserved_files'],
            'reference_interface_controls': controls['runs'],
            'finite_objective_exact_predicates': read(out / 'finite_objective_controls.json')['exact_predicates'],
            'finite_planner_rejected_inputs': read(out / 'finite_planner_input_controls.json')['rejected'],
            'fresh_load_fixture_count': read(out / 'load_planning/benchmark_receipt.json')['fixture_count'],
            'fresh_load_prospective_cases': read(out / 'load_planning/benchmark_receipt.json')['prospective_cases'],
            'fresh_load_arm_cases': read(out / 'load_planning/benchmark_receipt.json')['arm_cases'],
            'fresh_load_independent_predicates': read(out / 'load_planning/independent_verification.json')['exact_predicates'],
            'dynamic_fixtures': read(out / 'dynamic_acquisition/dynamic_receipt.json')['fixtures'],
            'dynamic_step_traces': read(out / 'dynamic_acquisition/dynamic_receipt.json')['step_transient_records'],
            'dynamic_laplace_scalar_settings': read(out / 'dynamic_acquisition/dynamic_receipt.json')['real_laplace_scalar_settings'],
            'dynamic_operator_enclosures': read(out / 'dynamic_acquisition/dynamic_receipt.json')['operator_enclosures'],
            'dynamic_independent_predicates': read(out / 'dynamic_acquisition/independent_receipt.json')['checks']}
        receipt['status'] = 'passed'
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
        receipt['status'] = 'failed'
        receipt['reason'] = str(exc)
    receipt['elapsed_seconds'] = round(time.monotonic() - started, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'stages'}, indent=2))
    return 0 if receipt['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
