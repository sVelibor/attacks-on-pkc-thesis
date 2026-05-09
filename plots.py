# plots.py
# Generate thesis figures from experimental results (reads results.csv)

import csv
import sys
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


# ── CSV loader ────────────────────────────────────────────────────────────────

def load_results(csv_path='results.csv'):
    """Read results.csv and return a dict keyed by attack name.

    Each value is a list of row dicts (parameter_value, trials, successes,
    avg_time_seconds), sorted by parameter_value ascending.

    Exits with an error message if the file is missing or contains no data rows.
    """
    try:
        with open(csv_path, newline='') as f:
            rows = list(csv.DictReader(f))
    except FileNotFoundError:
        print(f"[!] '{csv_path}' not found.")
        print("    Run  python main.py  first to generate the data file.")
        sys.exit(1)

    if not rows:
        print(f"[!] '{csv_path}' exists but contains no data rows.")
        print("    Run  python main.py  first to populate it.")
        sys.exit(1)

    results = {}
    for row in rows:
        attack = row['attack']
        entry = {
            'parameter_value':   int(row['parameter_value']),
            'trials':            int(row['trials']),
            'successes':         int(row['successes']),
            'avg_time_seconds':  float(row['avg_time_seconds']),
        }
        results.setdefault(attack, []).append(entry)

    for attack in results:
        results[attack].sort(key=lambda r: r['parameter_value'])

    return results


def _extract(rows):
    """Return (parameter_values, avg_times) lists from a list of row dicts."""
    xs = [r['parameter_value']  for r in rows]
    ys = [r['avg_time_seconds'] for r in rows]
    return xs, ys


# ── write_results_txt ─────────────────────────────────────────────────────────

def write_results_txt(csv_path, txt_path):
    """Write a human-readable summary of csv_path to txt_path."""
    results = load_results(csv_path)

    attack_labels = {
        'wiener':          "Wiener's Attack",
        'hastad':          "Håstad's Broadcast Attack",
        'coppersmith':     "Coppersmith's Short Pad Attack",
        'bsgs':            "Baby-step Giant-step (BSGS)",
        'pohlig_hellman':  "Pohlig-Hellman",
        'pollard_rho':     "Pollard's Rho (ECDLP)",
        'smart':           "Smart's Attack",
    }
    param_labels = {
        'key_bits': 'Key size (bits)',
        'q_bits':   'Subgroup order (bits)',
        'p_bits':   'Curve/prime size (bits)',
    }

    param_name_map = {}
    try:
        with open(csv_path, newline='') as f:
            for row in csv.DictReader(f):
                param_name_map.setdefault(row['attack'], row['parameter_name'])
    except FileNotFoundError:
        pass

    order = ['wiener', 'hastad', 'coppersmith',
             'bsgs', 'pohlig_hellman',
             'pollard_rho', 'smart']

    lines = []
    lines.append("EXPERIMENTAL RESULTS")
    lines.append("Attacks on Public-Key Encryption Schemes")
    lines.append("Velibor Smilevski, FAMNIT, 2026")
    lines.append("=" * 60)

    for attack in order:
        if attack not in results:
            continue
        label  = attack_labels.get(attack, attack)
        pname  = param_name_map.get(attack, 'parameter')
        plabel = param_labels.get(pname, pname)

        lines.append("")
        lines.append(label)
        lines.append("-" * len(label))
        lines.append(
            f"  {plabel:<24}  {'Trials':>6}  {'Success':>7}  {'Avg time (s)':>14}"
        )

        for r in results[attack]:
            lines.append(
                f"  {r['parameter_value']:<24}  "
                f"{r['trials']:>6}  {r['successes']:>7}  "
                f"{r['avg_time_seconds']:>14.6f}"
            )

    lines.append("")
    lines.append("=" * 60)

    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


# ── Figure generation (only when run directly) ────────────────────────────────

def generate_figures(csv_path='results.csv'):
    """Load results.csv and write all figure PDFs/PNGs."""
    results = load_results(csv_path)

    wiener_xs,      wiener_ys      = _extract(results['wiener'])
    hastad_xs,      hastad_ys      = _extract(results['hastad'])
    coppersmith_xs, coppersmith_ys = _extract(results['coppersmith'])
    bsgs_xs,        bsgs_ys        = _extract(results['bsgs'])
    ph_xs,          ph_ys          = _extract(results['pohlig_hellman'])
    pollard_xs,     pollard_ys     = _extract(results['pollard_rho'])
    smart_xs,       smart_ys       = _extract(results['smart'])

    # Figure 1: Wiener
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(wiener_xs, wiener_ys,
            marker='o', color='steelblue', linewidth=2,
            markersize=7, label="Wiener's attack")
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("Wiener's Attack: Runtime vs. Key Size (5 trials each)")
    ax.set_xticks(wiener_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_wiener_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_wiener_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_wiener_runtime.pdf")
    plt.close()

    # Figure 2: Hastad
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(hastad_xs, hastad_ys,
            marker='s', color='darkorange', linewidth=2,
            markersize=7, label="Håstad's broadcast attack")
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("Håstad's Broadcast Attack: Runtime vs. Key Size (5 trials each)")
    ax.set_xticks(hastad_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_hastad_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_hastad_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_hastad_runtime.pdf")
    plt.close()

    # Figure 3: RSA comparison
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(wiener_xs, wiener_ys,
            marker='o', color='steelblue', linewidth=2,
            markersize=7, label="Wiener's attack")
    ax.plot(hastad_xs, hastad_ys,
            marker='s', color='darkorange', linewidth=2,
            markersize=7, label="Håstad's broadcast attack")
    ax.plot(coppersmith_xs, coppersmith_ys,
            marker='^', color='green', linewidth=2,
            markersize=10, label="Coppersmith's short pad attack",
            linestyle='none')
    ax.set_xlabel("RSA key size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("RSA Attack Runtime Comparison")
    ax.set_xticks([256, 512, 1024])
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_comparison.pdf", bbox_inches='tight')
    plt.savefig("fig_comparison.png", bbox_inches='tight')
    print("[+] Saved fig_comparison.pdf")
    plt.close()

    # Figure 4: BSGS
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(bsgs_xs, bsgs_ys,
            marker='o', color='mediumseagreen', linewidth=2,
            markersize=7, label="Baby-step Giant-step")
    ax.set_xlabel("Subgroup order size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("BSGS: Runtime vs. Subgroup Order Size (5 trials each)")
    ax.set_xticks(bsgs_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_bsgs_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_bsgs_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_bsgs_runtime.pdf")
    plt.close()

    # Figure 5: Pohlig-Hellman
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ph_xs, ph_ys,
            marker='s', color='mediumpurple', linewidth=2,
            markersize=7, label="Pohlig-Hellman")
    ax.set_xlabel("Smooth group order size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("Pohlig-Hellman: Runtime vs. Group Order Size (5 trials each)")
    ax.set_xticks(ph_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_pohlig_hellman_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_pohlig_hellman_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_pohlig_hellman_runtime.pdf")
    plt.close()

    # Figure 6: Pollard's rho
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(pollard_xs, pollard_ys,
            marker='^', color='tomato', linewidth=2,
            markersize=7, label="Pollard's rho (ECDLP)")
    ax.set_xlabel("Curve order size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("Pollard's Rho (ECDLP): Runtime vs. Curve Order Size (5 trials each)")
    ax.set_xticks(pollard_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_pollard_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_pollard_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_pollard_runtime.pdf")
    plt.close()

    # Figure 7: DH comparison
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bsgs_xs, bsgs_ys,
            marker='o', color='mediumseagreen', linewidth=2,
            markersize=7, label="Baby-step Giant-step")
    ax.plot(ph_xs, ph_ys,
            marker='s', color='mediumpurple', linewidth=2,
            markersize=7, label="Pohlig-Hellman")
    ax.set_xlabel("Group order size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("DH/ElGamal Attack Runtime Comparison")
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_dh_comparison.pdf", bbox_inches='tight')
    plt.savefig("fig_dh_comparison.png", bbox_inches='tight')
    print("[+] Saved fig_dh_comparison.pdf")
    plt.close()

    # Figure 8: Smart's attack
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(smart_xs, smart_ys,
            marker='D', color='slategray', linewidth=2,
            markersize=7, label="Smart's attack")
    ax.set_xlabel("Prime field size (bits)")
    ax.set_ylabel("Average runtime (seconds)")
    ax.set_title("Smart's Attack: Runtime vs. Prime Field Size (5 trials each)")
    ax.set_xticks(smart_xs)
    ax.legend()
    plt.tight_layout()
    plt.savefig("fig_smart_runtime.pdf", bbox_inches='tight')
    plt.savefig("fig_smart_runtime.png", bbox_inches='tight')
    print("[+] Saved fig_smart_runtime.pdf")
    plt.close()

    print("\n[+] All figures generated successfully.")
    print("    Copy the .pdf files to your Overleaf project.")


if __name__ == "__main__":
    generate_figures('results.csv')
