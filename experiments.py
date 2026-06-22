# experiments.py
# Statistical experiment runner for "Attacks on Public-Key Encryption Schemes".
# Velibor Smilevski, University of Primorska, FAMNIT, 2026.
#
# Re-runs all seven attacks at higher trial counts, records the per-trial
# timing of every run, and computes per parameter point the mean, median,
# sample standard deviation, and inter-quartile range (IQR).
#
# Rationale for the trial counts (see thesis Section 8.1):
#   * 50 trials per point for six of the seven attacks.  These attacks are
#     deterministic (or, for BSGS, deterministic given the instance) and run
#     in micro- to milliseconds, so 50 seeded trials give a stable mean and
#     a meaningful median/IQR.
#   * 20 trials for Coppersmith's attack: each trial pays for a symbolic
#     resultant computation, which is far more expensive than the other six.
#   * Pollard's rho is a Las Vegas algorithm whose run-to-run variance comes
#     from the random walk and the random target, NOT from the curve.  The
#     curve (the costly part to construct) is therefore generated once per
#     parameter point and reused across the 50 trials with a fresh secret and
#     a fresh walk seed each time, which isolates and measures exactly that
#     intrinsic variance.  Smart's attack is handled the same way.
#
# Outputs (overwrites previous files):
#   results.csv        -- raw per-trial timings, one row per trial
#   results_stats.csv  -- per-point summary statistics
#   results.txt        -- human-readable summary table

import csv
import io
import random
import statistics
import time
from contextlib import redirect_stdout

from wiener import wiener_attack, generate_weak_rsa_keypair
from hastad import hastad_attack, generate_broadcast_setup
from coppersmith import (coppersmith_short_pad_sympy,
                         generate_rsa_keypair_fixed_e)
from dh import generate_dh_params, dh_keygen, generate_smooth_dh_params
from dh_attacks import bsgs, pohlig_hellman
from ecc import find_curve_with_prime_order, find_anomalous_curve
from ecc_attacks import pollard_rho_ecdlp, smart_attack

# ── Configuration ─────────────────────────────────────────────────────────────

SEED               = 0
TRIALS             = 50    # default trials per parameter point
COPPERSMITH_TRIALS = 20    # symbolic-resultant cost: fewer trials (justified)

RAW_CSV   = 'results.csv'
STATS_CSV = 'results_stats.csv'
TXT       = 'results.txt'

# perf_counter gives a high-resolution monotonic clock; this is essential for
# the sub-millisecond attacks, where the legacy time.time() resolution would
# dominate the measurement.
_timer = time.perf_counter


class _DevNull(io.IOBase):
    """A write sink that discards everything (suppresses noisy attack prints)."""
    def writable(self):
        return True

    def write(self, s):
        return len(s)


# ── Statistics helper ─────────────────────────────────────────────────────────

def _stats(times):
    """Return (mean, median, std, iqr) in seconds for a list of timings."""
    mean   = statistics.mean(times)
    median = statistics.median(times)
    std    = statistics.stdev(times) if len(times) > 1 else 0.0
    if len(times) >= 2:
        q1, _q2, q3 = statistics.quantiles(times, n=4, method='inclusive')
        iqr = q3 - q1
    else:
        iqr = 0.0
    return mean, median, std, iqr


# ── Per-point trial runners ───────────────────────────────────────────────────
# Each returns a list of (success: bool, elapsed_seconds: float), one per trial.

def run_wiener_point(prime_bits, trials):
    out = []
    for trial in range(trials):
        rng = random.Random(SEED + trial)
        N, e, d_real, p, q = generate_weak_rsa_keypair(prime_bits, rng)
        t0 = _timer()
        d_rec = wiener_attack(e, N)
        dt = _timer() - t0
        out.append((d_rec == d_real, dt))
    return out


def run_hastad_point(prime_bits, trials, e=3):
    out = []
    for trial in range(trials):
        rng = random.Random(SEED + trial)
        moduli, _, _ = generate_broadcast_setup(e, prime_bits, rng)
        M = rng.randint(2, min(moduli) - 1)
        cts = [pow(M, e, N) for N in moduli]
        t0 = _timer()
        M_rec = hastad_attack(cts, moduli, e)
        dt = _timer() - t0
        out.append((M_rec == M, dt))
    return out


def run_coppersmith_point(prime_bits, pad_bits, trials, e=3):
    out = []
    sink = _DevNull()
    for trial in range(trials):
        rng = random.Random(SEED + trial)
        N, e_key, d, p, q = generate_rsa_keypair_fixed_e(e, prime_bits, rng)
        max_M = N >> (pad_bits + 4)
        M  = rng.randint(2, max(2, max_M))
        r1 = rng.randint(0, 2 ** pad_bits - 1)
        r2 = rng.randint(0, 2 ** pad_bits - 1)
        while r2 == r1:
            r2 = rng.randint(0, 2 ** pad_bits - 1)
        M1, M2 = M + r1, M + r2
        C1 = pow(M1, e_key, N)
        C2 = pow(M2, e_key, N)
        with redirect_stdout(sink):
            t0 = _timer()
            M1_rec = coppersmith_short_pad_sympy(N, e_key, C1, C2, pad_bits)
            dt = _timer() - t0
        out.append((M1_rec == M1, dt))
    return out


def run_bsgs_point(q_bits, trials):
    out = []
    for trial in range(trials):
        rng = random.Random(SEED + trial)
        p, q, g = generate_dh_params(q_bits, rng)
        x_real, h = dh_keygen(p, q, g, rng)
        t0 = _timer()
        x_rec = bsgs(g, h, p, q)
        dt = _timer() - t0
        out.append((x_rec is not None and pow(g, x_rec, p) == h, dt))
    return out


def run_pohlig_point(q_bits, trials):
    out = []
    for trial in range(trials):
        rng = random.Random(SEED + trial)
        p, q, g, fact = generate_smooth_dh_params(q_bits, rng)
        x_real = rng.randint(2, q - 1)
        h = pow(g, x_real, p)
        t0 = _timer()
        x_rec = pohlig_hellman(g, h, p, q, fact)
        dt = _timer() - t0
        out.append((x_rec is not None and pow(g, x_rec, p) == h, dt))
    return out


def run_pollard_point(p_bits, trials):
    # One prime-order curve per point (not timed); fresh secret + walk per trial.
    rng_curve = random.Random(SEED + p_bits * 1000)
    curve, G, q = find_curve_with_prime_order(p_bits, rng_curve)
    out = []
    for trial in range(trials):
        rng_t = random.Random(SEED + trial)
        k_real = rng_t.randint(2, q - 1)
        Q = curve.scalar_mult(k_real, G)
        t0 = _timer()
        k_rec = pollard_rho_ecdlp(curve, G, Q, q, seed=SEED + trial)
        dt = _timer() - t0
        out.append((k_rec is not None and curve.scalar_mult(k_rec, G) == Q, dt))
    return out


def run_smart_point(p_bits, trials):
    # One anomalous curve per point (not timed); fresh secret per trial.
    rng_curve = random.Random(SEED + p_bits * 1000)
    acurve, AG, ap = find_anomalous_curve(p_bits, rng_curve)
    out = []
    for trial in range(trials):
        rng_t = random.Random(SEED + trial)
        k_real = rng_t.randint(2, ap - 2)
        Q = acurve.scalar_mult(k_real, AG)
        t0 = _timer()
        k_rec = smart_attack(acurve, AG, Q)
        dt = _timer() - t0
        out.append((k_rec is not None and acurve.scalar_mult(k_rec, AG) == Q, dt))
    return out


# ── Experiment specification ──────────────────────────────────────────────────
# (attack, parameter_name, [(parameter_value, runner, trials), ...])

SPECS = [
    ('wiener', 'key_bits', [
        (256,  lambda: run_wiener_point(128, TRIALS), TRIALS),
        (512,  lambda: run_wiener_point(256, TRIALS), TRIALS),
        (1024, lambda: run_wiener_point(512, TRIALS), TRIALS),
    ]),
    ('hastad', 'key_bits', [
        (256,  lambda: run_hastad_point(128, TRIALS), TRIALS),
        (512,  lambda: run_hastad_point(256, TRIALS), TRIALS),
        (1024, lambda: run_hastad_point(512, TRIALS), TRIALS),
    ]),
    ('coppersmith', 'key_bits', [
        (256,  lambda: run_coppersmith_point(128, 8, COPPERSMITH_TRIALS),
         COPPERSMITH_TRIALS),
    ]),
    ('bsgs', 'q_bits', [
        (16, lambda: run_bsgs_point(16, TRIALS), TRIALS),
        (20, lambda: run_bsgs_point(20, TRIALS), TRIALS),
        (24, lambda: run_bsgs_point(24, TRIALS), TRIALS),
    ]),
    ('pohlig_hellman', 'q_bits', [
        (20, lambda: run_pohlig_point(20, TRIALS), TRIALS),
        (25, lambda: run_pohlig_point(25, TRIALS), TRIALS),
        (30, lambda: run_pohlig_point(30, TRIALS), TRIALS),
    ]),
    ('pollard_rho', 'p_bits', [
        (16, lambda: run_pollard_point(16, TRIALS), TRIALS),
        (18, lambda: run_pollard_point(18, TRIALS), TRIALS),
        (20, lambda: run_pollard_point(20, TRIALS), TRIALS),
    ]),
    ('smart', 'p_bits', [
        (12, lambda: run_smart_point(12, TRIALS), TRIALS),
        (14, lambda: run_smart_point(14, TRIALS), TRIALS),
        (16, lambda: run_smart_point(16, TRIALS), TRIALS),
    ]),
]

ATTACK_LABELS = {
    'wiener':          "Wiener's Attack",
    'hastad':          "Håstad's Broadcast Attack",
    'coppersmith':     "Coppersmith's Short Pad Attack",
    'bsgs':            "Baby-step Giant-step (BSGS)",
    'pohlig_hellman':  "Pohlig-Hellman",
    'pollard_rho':     "Pollard's Rho (ECDLP)",
    'smart':           "Smart's Attack",
}
PARAM_LABELS = {
    'key_bits': 'Key size (bits)',
    'q_bits':   'Group order (bits)',
    'p_bits':   'Curve/prime size (bits)',
}


# ── Output writers ────────────────────────────────────────────────────────────

def write_results_txt(stats_rows, txt_path=TXT):
    """Write a human-readable summary table of the per-point statistics."""
    lines = [
        "EXPERIMENTAL RESULTS",
        "Attacks on Public-Key Encryption Schemes",
        "Velibor Smilevski, FAMNIT, 2026",
        "=" * 78,
    ]
    by_attack = {}
    for row in stats_rows:
        by_attack.setdefault(row[0], []).append(row)

    for attack, _pname, _points in SPECS:
        if attack not in by_attack:
            continue
        label = ATTACK_LABELS.get(attack, attack)
        plabel = PARAM_LABELS.get(by_attack[attack][0][1], 'parameter')
        lines.append("")
        lines.append(label)
        lines.append("-" * len(label))
        lines.append(f"  {plabel:<22} {'N':>4} {'Succ':>5} "
                     f"{'Mean (s)':>12} {'Median (s)':>12} "
                     f"{'Std (s)':>12} {'IQR (s)':>12}")
        for (_a, _pn, pval, n, succ, mean, median, std, iqr) in by_attack[attack]:
            lines.append(f"  {pval:<22} {n:>4} {succ:>5} "
                         f"{mean:>12.6f} {median:>12.6f} "
                         f"{std:>12.6f} {iqr:>12.6f}")
    lines.append("")
    lines.append("=" * 78)
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  STATISTICAL EXPERIMENT SUITE")
    print("  Attacks on Public-Key Encryption Schemes")
    print(f"  {TRIALS} trials per point ({COPPERSMITH_TRIALS} for Coppersmith)")
    print("=" * 60)

    raw_rows   = []
    stats_rows = []

    for attack, pname, points in SPECS:
        for pval, runner, _ntr in points:
            print(f"\n>>> {ATTACK_LABELS[attack]}  ({pname}={pval}) ...", flush=True)
            trials_data = runner()
            times = [t for (_ok, t) in trials_data]
            succ  = sum(1 for (ok, _t) in trials_data if ok)
            for i, (ok, t) in enumerate(trials_data):
                raw_rows.append([attack, pname, pval, i, int(ok), f"{t:.9f}"])
            mean, median, std, iqr = _stats(times)
            stats_rows.append([attack, pname, pval, len(times), succ,
                               mean, median, std, iqr])
            print(f"    n={len(times)} success={succ}/{len(times)}  "
                  f"mean={mean:.6f}s  median={median:.6f}s  iqr={iqr:.6f}s")

    with open(RAW_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['attack', 'parameter_name', 'parameter_value',
                    'trial', 'success', 'time_seconds'])
        w.writerows(raw_rows)

    with open(STATS_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['attack', 'parameter_name', 'parameter_value',
                    'n_trials', 'successes',
                    'mean_s', 'median_s', 'std_s', 'iqr_s'])
        for r in stats_rows:
            w.writerow(r[:5] + [f"{x:.9f}" for x in r[5:]])

    write_results_txt(stats_rows)

    print(f"\n[+] Wrote {RAW_CSV} ({len(raw_rows)} raw trial rows), "
          f"{STATS_CSV}, and {TXT}.")
    return stats_rows


if __name__ == "__main__":
    main()
