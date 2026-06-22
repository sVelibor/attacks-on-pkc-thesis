# main.py
# Main entry point - reproduces all experiments and regenerates the figures.
# Attacks on Public-Key Encryption Schemes
# Velibor Smilevski, University of Primorska FAMNIT, 2026
#
# This runs the full statistical experiment suite (experiments.py) and then
# regenerates the eight thesis figures (plots.py).  The individual per-module
# demonstrations remain available by running each module directly, e.g.
#     python wiener.py
#     python ecc_attacks.py

import experiments
import plots


def print_banner():
    print("\n" + "=" * 60)
    print("  ATTACKS ON PUBLIC-KEY ENCRYPTION SCHEMES")
    print("  Velibor Smilevski")
    print("  University of Primorska, FAMNIT, 2026")
    print("=" * 60)


if __name__ == "__main__":
    print_banner()

    # 1. Run all seven attacks, collect per-trial timings and statistics.
    #    Writes results.csv (raw per-trial), results_stats.csv, results.txt.
    experiments.main()

    # 2. Regenerate the eight figures (median runtime with IQR error bars).
    print("\n" + "=" * 60)
    print("  REGENERATING FIGURES")
    print("=" * 60)
    plots.generate_figures(experiments.RAW_CSV)

    print("\n[+] Done. Data in results.csv / results_stats.csv / results.txt,")
    print("    figures in fig_*.pdf (copy the PDFs to the Overleaf project).")
