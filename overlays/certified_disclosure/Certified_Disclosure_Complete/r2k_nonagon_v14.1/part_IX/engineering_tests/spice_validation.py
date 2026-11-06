"""Independent ngspice validation of the Nonagon electrical acquisition.

All currents and voltages here are external circuit simulations, not laboratory
measurements. The source fixture is retained as exact rational conductances.
Run with --ngspice PATH or with ngspice on PATH. No automatic installation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
CORE = next((p for p in [HERE.parent, HERE.parent / "part_VIII", HERE.parent / "part8_work" / "part_VIII"]
             if (p / "intervention.py").is_file()), HERE.parent / "part_VIII")
sys.path.insert(0, str(CORE))
from intervention import acquisition_plan, response_vector, strings
from campaign import retained_twins, random_network

SCALE = F(1, 1000)  # physical siemens per normalized conductance
CURRENT = F(1, 10000)  # physical ampere
CAPACITANCE = 1e-6  # one microfarad at every node for the transient experiment


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+"\n")


def fmt(x):
    return format(float(x), ".17g")


def resistor_lines(K, loads):
    lines = []
    n = len(K)
    for i in range(n):
        for j in range(i+1, n):
            g = -K[i][j]*SCALE
            if g > 0:
                lines.append(f"Redge_{i+1}_{j+1} n{i+1} n{j+1} {fmt(1/g)}")
        g = sum(K[i])*SCALE
        if g > 0:
            lines.append(f"Rground_{i+1} n{i+1} 0 {fmt(1/g)}")
    for j, load in loads.items():
        g = F(load)*SCALE
        lines.append(f"Rload_{int(j)+1} n{int(j)+1} 0 {fmt(1/g)}")
    return lines


def run_spice(binary, circuit, log, env):
    command = [str(binary), "-n", "-b", "-o", str(log), str(circuit)]
    result = subprocess.run(command, text=True, capture_output=True, env=env, timeout=45)
    # ngspice control-mode .op may return nonzero after quit with valid output;
    # parse explicit results, but never accept a recorded fatal error.
    output = log.read_text() if log.exists() else ""
    if result.returncode or re.search(r"(?im)^fatal error|^error:", output):
        raise RuntimeError(f"ngspice failed: {circuit.name}, rc={result.returncode}: {output[-1500:]}")
    return command, output, result.stdout, result.stderr


def numeric_reconstruct(n, readings):
    """Float inverse consumes only normalized scalar readouts and public labels."""
    h = readings["baseline"]
    G = np.zeros((n, n)); G[0, 0] = h
    for j in range(1, n):
        d1 = h-readings[f"single_{j}_a"]
        d2 = h-readings[f"single_{j}_b"]
        if min(d1, d2) <= 0:
            raise ValueError("unobservable or inconsistent single load")
        y1, y2 = 1/d1, 2/d2
        slope, intercept = y2-y1, 2*y1-y2
        if min(slope, intercept) <= 0:
            raise ValueError("nonpositive single-load extraction")
        G[0, j] = G[j, 0] = math.sqrt(1/intercept)
        G[j, j] = slope/intercept
    for j in range(1, n):
        for k in range(j+1, n):
            d = h-readings[f"pair_{j}_{k}"]
            a, b, p, q = G[0,j], G[0,k], G[j,j], G[k,k]
            A = a*a*(1+q)+b*b*(1+p)
            B = (1+p)*(1+q)
            rad = (a*b)**2-d*A+d*d*B
            if rad < -1e-14:
                raise ValueError("negative physical-root discriminant")
            # Rationalized smaller root.
            c = (A-d*B)/(a*b+math.sqrt(max(0,rad)))
            G[j,k] = G[k,j] = c
    return np.linalg.inv(G)


def fixtures():
    K, Kt, *_ = retained_twins()
    values = [("v7_path", K), ("v7_triangle", Kt)]
    rng = random.Random(2026103009)
    for case in range(12):
        n = 3+case%6
        values.append((f"random_{case:02d}_n{n}", random_network(n, rng)))
    return values


def dc_campaign(binary, env, output):
    cases = []
    total_observations = 0
    for name, K in fixtures():
        n = len(K); records=[]; readings={}
        plan = acquisition_plan(n)
        heldout = [{"id":"heldout_0","loads":{"0":"1/3",str(n-1):"2/5"}},
                   {"id":"heldout_1","loads":{"1":"3/7",str(n-1):"5/11"}},
                   {"id":"heldout_2","loads":{"0":"1/13",str(n-1):"7/11"}}]
        folder = output/"netlists"/name
        folder.mkdir(parents=True, exist_ok=True)
        for row in plan+heldout:
            cid, loads = row["id"], row["loads"]
            circuit, log = folder/f"{cid}.cir", folder/f"{cid}.log"
            lines = [f"Nonagon external DC circuit {name} {cid}",
                     ".options reltol=1e-12 abstol=1e-15 vntol=1e-12 gmin=1e-18",
                     f"Iinput 0 n1 DC {fmt(CURRENT)}"]+resistor_lines(K,loads)
            lines += [".control", "set numdgt=17", "op",
                      "print "+" ".join(f"v(n{j+1})" for j in range(n)),
                      "quit", ".endc", ".end"]
            circuit.write_text("\n".join(lines)+"\n")
            command, text, stdout, stderr = run_spice(binary,circuit,log,env)
            parsed = {int(i)-1:float(v) for i,v in re.findall(r"v\(n(\d+)\)\s*=\s*([-+0-9.eE]+)",text,re.I)}
            if len(parsed) != n:
                raise RuntimeError(f"Incomplete ngspice voltage output in {log}")
            exact_u = response_vector(K,loads)
            exact_v = [u*CURRENT/SCALE for u in exact_u]
            volts = [parsed[j] for j in range(n)]
            h = volts[0]*float(SCALE/CURRENT)
            readings[cid] = h
            records.append({"id":cid,"loads":loads,"kind":"acquisition" if cid in {r['id'] for r in plan} else "heldout",
                "ngspice_voltage_V":volts,"exact_voltage_V":[str(v) for v in exact_v],
                "normalized_response":h,"exact_normalized_response":str(exact_u[0]),
                "max_voltage_error_V":max(abs(v-float(e)) for v,e in zip(volts,exact_v)),
                "command":command,"netlist":str(circuit.relative_to(output)),
                "raw_log":str(log.relative_to(output)),"netlist_sha256":hashlib.sha256(circuit.read_bytes()).hexdigest()})
            total_observations += 1
        recovered = numeric_reconstruct(n,readings)
        true = np.array(K,dtype=float)
        error = float(np.max(np.abs(recovered-true)))
        heldout_error = 0.0
        for row in heldout:
            A = recovered.copy()
            for j,t in row['loads'].items(): A[int(j),int(j)]+=float(F(t))
            heldout_error=max(heldout_error,abs(float(np.linalg.inv(A)[0,0])-readings[row['id']]))
        cases.append({"name":name,"n":n,"exact_K_normalized":strings(K),
                      "acquisition_configurations":len(plan),"heldout_configurations":len(heldout),
                      "max_simultaneous_acquisition_shunts":max(len(row['loads']) for row in plan),
                      "max_simultaneous_validation_shunts":max(len(row['loads']) for row in heldout),
                      "observations":records,"max_voltage_error_V":max(r['max_voltage_error_V'] for r in records),
                      "recovered_K_float":recovered.tolist(),"max_K_coefficient_error":error,
                      "heldout_max_normalized_response_error":heldout_error,
                      "passed":bool(error<1e-6 and max(r['max_voltage_error_V'] for r in records)<1e-10)})
    return {"fixture_count":len(cases),"observation_count":total_observations,"cases":cases,
        "max_voltage_error_V":max(c['max_voltage_error_V'] for c in cases),
        "max_K_coefficient_error":max(c['max_K_coefficient_error'] for c in cases),
        "all_passed":all(c['passed'] for c in cases)}


def transient_campaign(binary, env, output):
    chosen=[fixtures()[i] for i in (0,1,5,7)]
    results=[]
    for name,K in chosen:
        n=len(K); A=np.array(K,dtype=float)*float(SCALE)/CAPACITANCE
        lamb,V=np.linalg.eigh(A)
        steady=np.linalg.solve(np.array(K,dtype=float)*float(SCALE),np.eye(n)[0]*float(CURRENT))
        tau=1/lamb[0]; maxstep=min(tau/200,1/lamb[-1]/25);stop=20*tau
        folder=output/'netlists'/name;folder.mkdir(parents=True,exist_ok=True)
        circuit,log,series=folder/'transient.cir',folder/'transient.log',folder/'transient.dat'
        lines=[f'Nonagon independent RC transient {name}',
            '.options reltol=1e-10 abstol=1e-16 vntol=1e-12 gmin=1e-18 method=gear maxord=2',
            f'Iinput 0 n1 DC {fmt(CURRENT)}']+resistor_lines(K,{})
        lines += [f'Cground_{i+1} n{i+1} 0 {fmt(CAPACITANCE)} IC=0' for i in range(n)]
        lines += ['.control','set numdgt=17','set wr_singlescale','set wr_vecnames',
                  f'tran {fmt(maxstep)} {fmt(stop)} 0 {fmt(maxstep)} uic',
                  'wrdata '+str(series)+' '+' '.join(f'v(n{i+1})' for i in range(n)),
                  'quit','.endc','.end']
        circuit.write_text('\n'.join(lines)+'\n')
        command,_,_,_=run_spice(binary,circuit,log,env)
        values=np.loadtxt(series,skiprows=1)
        times,observed=values[:,0],values[:,1:]
        ref=steady[None,:]-(np.exp(-times[:,None]*lamb[None,:])*(V.T@steady)[None,:])@V.T
        maxerror=float(np.max(abs(ref-observed)))
        checktimes=np.linspace(0,stop,17)
        eig_expm=max(float(np.max(abs(steady-expm(-A*t)@steady-(steady-V@(np.exp(-lamb*t)*(V.T@steady)))))) for t in checktimes)
        voltage_tolerance=1e-6
        fail=np.where(abs(observed[:,0]-steady[0])>voltage_tolerance)[0]
        settle_idx=(fail[-1]+1) if len(fail) else 0
        settled=bool(settle_idx<len(times))
        port_settle=float(times[settle_idx]) if settled else None
        bound_t=float(math.log(np.linalg.norm(steady)/voltage_tolerance)/lamb[0])
        results.append({'name':name,'n':n,'capacitance_each_F':CAPACITANCE,'current_A':str(CURRENT),
            'gscale_S':str(SCALE),'exact_K_normalized':strings(K),'command':command,
            'netlist':str(circuit.relative_to(output)),'raw_log':str(log.relative_to(output)),
            'raw_trace':str(series.relative_to(output)),'trace_rows':len(times),
            'maxstep_s':maxstep,'stop_s':stop,'slowest_time_constant_s':float(tau),
            'steady_voltage_V':steady.tolist(),'port_settling_tolerance_V':voltage_tolerance,
            'simulated_port_settling_time_s':port_settle,'all_node_eigenvalue_envelope_time_s':bound_t,
            'max_trace_error_against_matrix_exponential_V':maxerror,
            'eigen_exp_vs_scipy_expm_max_V':eig_expm,'settled':settled,
            'passed':bool(maxerror<1e-7 and eig_expm<1e-12 and settled)})
    return {'fixture_count':len(results),'cases':results,
        'max_trace_error_V':max(r['max_trace_error_against_matrix_exponential_V'] for r in results),
        'all_passed':all(r['passed'] for r in results)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ngspice',type=Path)
    ap.add_argument('--output',type=Path,default=HERE)
    ap.add_argument('--library-path',type=Path)
    args=ap.parse_args()
    binary=args.ngspice or shutil.which('ngspice')
    if not binary:raise SystemExit('ngspice is required; no circuit-solver substitute is used')
    env=os.environ.copy()
    if args.library_path:env['LD_LIBRARY_PATH']=str(args.library_path)
    scripts = args.output/'spice_config'
    scripts.mkdir(parents=True,exist_ok=True)
    (scripts/'spinit').write_text('* Reproducible passive R/C tests: no optional XSPICE/OSDI models are loaded.\n')
    env['SPICE_SCRIPTS']=str(scripts)
    version=subprocess.run([str(binary),'--version'],env=env,capture_output=True,text=True,check=True).stdout
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    dc=dc_campaign(binary,env,args.output)
    tr=transient_campaign(binary,env,args.output)
    result={'schema':'nonagon.ngspice.validation.v1','status':'PASSED' if dc['all_passed'] and tr['all_passed'] else 'FAILED',
        'scope':'Independent numerical circuit simulation; no physical apparatus connected',
        'ngspice_version':version,'binary':str(binary),'binary_sha256':hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        'conductance_scale_S':str(SCALE),'source_current_A':str(CURRENT),
        'solver_control':'DC reltol1e-12 abstol1e-15 vntol1e-12 gmin1e-18; transient gear2 reltol1e-10',
        'dc':dc,'transient':tr,'elapsed_seconds':time.monotonic()-started}
    write(args.output/'spice_results.json',result)
    print(json.dumps({'status':result['status'],'dc_fixtures':dc['fixture_count'],'dc_observations':dc['observation_count'],
        'max_dc_voltage_error_V':dc['max_voltage_error_V'],'max_reconstruction_coefficient_error':dc['max_K_coefficient_error'],
        'RC_transients':tr['fixture_count'],'max_transient_error_V':tr['max_trace_error_V'],'elapsed_seconds':result['elapsed_seconds']},indent=2))
    if result['status']!='PASSED':raise SystemExit(1)


if __name__=='__main__':main()
