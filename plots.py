# plots.py
# Generate thesis figures from experimental results.
#
# Reads the raw per-trial timings in results.csv (written by experiments.py),
# and for each parameter point plots the median runtime with error bars showing
# the inter-quartile range (IQR).  The IQR makes Pollard's rho's run-to-run
# variance visible while remaining negligible for the deterministic attacks.

import csv
import sys
import statistics

import matplotlib
matplotlib.use('Agg')           # headless backend (no display needed)
import matplotlib.pyplot as plt

# ── Styling ──────────────────────────────────────────────────────────────────

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'axes.grid': True,
    'grid.alpha': 0.3,
})


# ── Raw-data loader ───────────────────────────────────────────────────────────

def load_results(csv_path='results.csv'):
    """Read the raw per-trial results.csv and return a dict keyed by attack.

    Each value is a list of per-point dicts with keys:
        param, n, median, q1, q3, mean
    sorted by parameter value ascending.
    """
    try:
        with open(csv_path, newline='') as f:
            rows = list(csv.DictReader(f))
    except FileNotFoundError:
        print(f"[!] '{csv_path}' not found.")
        print("    Run  python experiments.py  first to generate the data file.")
        sys.exit(1)

    if not rows or 'time_seconds' not in rows[0]:
        print(f"[!] '{csv_path}' is empty or not in the expected per-trial format.")
        print("    Run  python experiments.py  first to (re)generate it.")
        sys.exit(1)

    grouped = {}   # (attack, param) -> list of times
    for row in rows:
        key = (row['attack'], int(row['parameter_value']))
        grouped.setdefault(key, []).append(float(row['time_seconds']))

    results = {}
    for (attack, param), times in grouped.items():
        median = statistics.median(times)
        if len(times) >= 2:
            q1, _q2, q3 = statistics.quantiles(times, n=4, method='inclusive')
        else:
            q1 = q3 = median
        results.setdefault(attack, []).append({
            'param':  param,
            'n':      len(times),
            'median': median,
            'q1':     q1,
            'q3':     q3,
            'mean':   statistics.mean(times),
        })

    for attack in results:
        results[attack].sort(key=lambda r: r['param'])
    return results


def _xyerr(rows):
    """Return (xs, medians, asymmetric_yerr) for an errorbar plot."""
    xs      = [r['param']  for r in rows]
    medians = [r['median'] for r in rows]
    lower   = [r['median'] - r['q1'] for r in rows]
    upper   = [r['q3'] - r['median'] for r in rows]
    return xs, medians, [lower, upper]


def _ntrials(rows):
    return rows[0]['n'] if rows else 0


# ── Figure generation ─────────────────────────────────────────────────────────

def generate_figures(csv_path='results.csv'):
    """Load results.csv and write all figure PDFs/PNGs."""
    results = load_results(csv_path)

    wiener      = results['wiener']
    hastad      = results['hastad']
    coppersmith = results['coppersmith']
    bsgs        = results['bsgs']
    ph          = results['pohlig_hellman']
    pollard     = results['pollard_rho']
    smart       = results['smart']

    def save(name):
        plt.tight_layout()
        plt.savefig(f"{name}.pdf", bbox_inches='tight')
        plt.savefig(f"{name}.png", bbox_inches='tight')
        print(f"[+] Saved {name}.pdf")
        plt.close()

    # Figure 1: Wiener
    xs, ys, err = _xyerr(wiener)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='o', color='steelblue', linewidth=2,
                markersize=7, capsize=3, label="Wiener's attack")
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"Wiener's Attack: Runtime vs.\\ Key Size "
                 f"(median of {_ntrials(wiener)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_wiener_runtime")

    # Figure 2: Hastad
    xs, ys, err = _xyerr(hastad)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='s', color='darkorange', linewidth=2,
                markersize=7, capsize=3, label="Håstad's broadcast attack")
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"Håstad's Broadcast Attack: Runtime vs. Key Size "
                 f"(median of {_ntrials(hastad)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_hastad_runtime")

    # Figure 3: RSA comparison (log scale)
    wx, wy, werr = _xyerr(wiener)
    hx, hy, herr = _xyerr(hastad)
    cx, cy, cerr = _xyerr(coppersmith)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(wx, wy, yerr=werr, marker='o', color='steelblue', linewidth=2,
                markersize=7, capsize=3, label="Wiener's attack")
    ax.errorbar(hx, hy, yerr=herr, marker='s', color='darkorange', linewidth=2,
                markersize=7, capsize=3, label="Håstad's broadcast attack")
    ax.errorbar(cx, cy, yerr=cerr, marker='^', color='green',
                markersize=10, capsize=3, linestyle='none',
                label="Coppersmith's short pad attack")
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title("RSA Attack Runtime Comparison")
    ax.set_xticks([256, 512, 1024])
    ax.set_yscale('log')
    ax.legend()
    save("fig_comparison")

    # Figure 4: BSGS
    xs, ys, err = _xyerr(bsgs)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='o', color='mediumseagreen', linewidth=2,
                markersize=7, capsize=3, label="Baby-step Giant-step")
    ax.set_xlabel("Subgroup order size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"BSGS: Runtime vs. Subgroup Order Size "
                 f"(median of {_ntrials(bsgs)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_bsgs_runtime")

    # Figure 5: Pohlig-Hellman
    xs, ys, err = _xyerr(ph)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='s', color='mediumpurple', linewidth=2,
                markersize=7, capsize=3, label="Pohlig-Hellman")
    ax.set_xlabel("Smooth group order size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"Pohlig-Hellman: Runtime vs. Group Order Size "
                 f"(median of {_ntrials(ph)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_pohlig_hellman_runtime")

    # Figure 6: Pollard's rho (variance is the point here)
    xs, ys, err = _xyerr(pollard)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='^', color='tomato', linewidth=2,
                markersize=7, capsize=4, elinewidth=1.5,
                label="Pollard's rho (ECDLP)")
    ax.set_xlabel("Curve order size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"Pollard's Rho (ECDLP): Runtime vs. Curve Order Size "
                 f"(median of {_ntrials(pollard)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_pollard_runtime")

    # Figure 7: DH comparison (log scale)
    bx, by, berr = _xyerr(bsgs)
    px, py, perr = _xyerr(ph)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(bx, by, yerr=berr, marker='o', color='mediumseagreen', linewidth=2,
                markersize=7, capsize=3, label="Baby-step Giant-step")
    ax.errorbar(px, py, yerr=perr, marker='s', color='mediumpurple', linewidth=2,
                markersize=7, capsize=3, label="Pohlig-Hellman")
    ax.set_xlabel("Group order size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title("DH/ElGamal Attack Runtime Comparison")
    ax.set_yscale('log')
    ax.legend()
    save("fig_dh_comparison")

    # Figure 8: Smart's attack
    xs, ys, err = _xyerr(smart)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(xs, ys, yerr=err, marker='D', color='slategray', linewidth=2,
                markersize=7, capsize=3, label="Smart's attack")
    ax.set_xlabel("Prime field size (bits)")
    ax.set_ylabel("Median runtime (seconds)")
    ax.set_title(f"Smart's Attack: Runtime vs. Prime Field Size "
                 f"(median of {_ntrials(smart)} trials, error bars: IQR)")
    ax.set_xticks(xs)
    ax.legend()
    save("fig_smart_runtime")

    print("\n[+] All figures generated successfully.")
    print("    Copy the .pdf files to your Overleaf project.")


if __name__ == "__main__":
    generate_figures('results.csv')
