#!/usr/bin/env python3
"""Find the sourcing adjustment cost phi_Z that delivers a target short-run trade elasticity.

The statistic is the short-run Armington (trade) elasticity: the US impact (t = 0) elasticity of
imported relative to domestic quantities with respect to their relative price, averaged across goods
sectors and sources with base-period trade weights (column `ae`; its long-run value is the
trade-weighted sigma, about 3.25). The experiment is an unanticipated, permanent, symmetric 25%
tariff on all goods: every country taxes imports from every other country
(-t 25 -s 6 -c 2 -r 2 -d 0 -a 4). Baseline phi_Z = 2.9 gives about 0.48.

(The import-value elasticity `te` in a unilateral US tariff is not usable: it is dominated by the
dollar's appreciation and barely moves with phi_Z.)

usage: python3 calibrate_phiz.py --target 0.5 [--target 0.2] [--bin ./bin/model.bin] [--rundir .]
  --rundir must contain output/ and ../python/output/iomat.txt (as programs/c does).
Each evaluation is one model solve, run with -z <phi>; its output files carry the suffix _pz<phi>.
Evaluations are cached in <rundir>/output/phiz_calibration_ae.csv, so reruns and later targets reuse
them. The values found for targets 1.0 and 0.2 are hard-coded as PHIZ_LO / PHIZ_HI (src/globals.h)
and PHIZ (../python/model_io.py); the paper's sensitivity runs use -z lo / -z hi (run_all.sh). The search is a secant method on (log phi, log elasticity), stopping at |te - target| < tol.
"""
import argparse
import csv
import math
import os
import subprocess

BASE_ARGS = ['-t', '25', '-s', '6', '-c', '2', '-r', '2', '-d', '0', '-a', '4']
STAT = 'ae'


def fmt(phi):
    return '%g' % phi            # same as the C code's "_pz%g" suffix


def run_model(phi, args):
    out = os.path.join(args.rundir, 'output')
    log = os.path.join(out, 'log_phiz_%s.txt' % fmt(phi))
    with open(log, 'w') as f:
        rc = subprocess.run([args.bin] + BASE_ARGS + ['-z', fmt(phi)], cwd=args.rundir,
                            stdout=f, stderr=subprocess.STDOUT).returncode
    path = os.path.join(out, 'vars1_usa_t25_s6_c2_r2_d0_a4_pz%s.csv' % fmt(phi))
    if rc != 0 or not os.path.exists(path):
        raise SystemExit('model failed at phi_Z = %s (see %s)' % (fmt(phi), log))
    with open(path) as f:
        row = next(csv.DictReader(f))           # t = 0: impact period
    return float(row[STAT])


def load_cache(path):
    cache = {}
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                cache[float(r['phiz'])] = float(r[STAT])
    return cache


def evaluate(phi, cache, args, cache_path):
    phi = float(fmt(phi))
    if phi not in cache:
        cache[phi] = run_model(phi, args)
        with open(cache_path, 'w') as f:
            f.write('phiz,%s\n' % STAT)
            for k in sorted(cache):
                f.write('%s,%.10f\n' % (fmt(k), cache[k]))
        print('  phi_Z = %-10s %s = %.4f' % (fmt(phi), STAT, cache[phi]), flush=True)
    return cache[phi]


def solve_target(target, cache, args, cache_path):
    print('target %s = %g' % (STAT, target), flush=True)
    if not cache:
        evaluate(2.9, cache, args, cache_path)
    for _ in range(args.maxit):
        pts = sorted(cache.items(), key=lambda kv: abs(math.log(kv[1] / target)))
        phi0, e0 = pts[0]
        if abs(e0 - target) < args.tol:
            return phi0, e0
        if len(pts) == 1:
            phi1 = phi0 * e0 / target           # first step: assume elasticity roughly proportional to 1/phi
        else:
            phi_b, e_b = pts[1]
            slope = (math.log(e0) - math.log(e_b)) / (math.log(phi0) - math.log(phi_b))
            if not slope < 0:                   # elasticity should fall with phi_Z
                slope = -1.0
            phi1 = math.exp(math.log(phi0) + (math.log(target) - math.log(e0)) / slope)
        phi1 = float('%.4g' % phi1)
        if phi1 in cache:
            phi1 = float('%.5g' % (phi1 * (1 + 1e-3)))
        evaluate(phi1, cache, args, cache_path)
    phi0, e0 = min(cache.items(), key=lambda kv: abs(kv[1] - target))
    if abs(e0 - target) >= args.tol:
        raise SystemExit('no convergence for target %g; closest phi_Z = %s (%s = %.4f)' % (target, fmt(phi0), STAT, e0))
    return phi0, e0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', type=float, action='append', required=True)
    ap.add_argument('--bin', default='./bin/model.bin')
    ap.add_argument('--rundir', default='.')
    ap.add_argument('--tol', type=float, default=0.002)
    ap.add_argument('--maxit', type=int, default=10)
    args = ap.parse_args()
    args.bin = os.path.abspath(args.bin)
    cache_path = os.path.join(args.rundir, 'output', 'phiz_calibration_%s.csv' % STAT)
    cache = load_cache(cache_path)
    results = [(t,) + solve_target(t, cache, args, cache_path) for t in args.target]
    print('\nresult:')
    for t, phi, e in results:
        print('  target %-5g phi_Z = %-10s achieved %s = %.4f' % (t, fmt(phi), STAT, e))


if __name__ == '__main__':
    main()
