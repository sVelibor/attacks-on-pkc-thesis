# main.py
# Main entry point - runs all implemented attacks and collects results
# Attacks on Public-Key Encryption Schemes
# Velibor Smilevski, University of Primorska FAMNIT, 2026

import csv
import random
import time

from wiener import run_demo as wiener_demo, run_scaling_experiment as wiener_scaling
from hastad import run_demo as hastad_demo, run_scaling_experiment as hastad_scaling
from coppersmith import (run_demo as coppersmith_demo,
                          run_scaling_experiment as coppersmith_scaling)
from dh import run_demo as dh_demo
from dh_attacks import (run_demo as dh_attacks_demo,
                         run_scaling_experiment_bsgs,
                         run_scaling_experiment_pohlig_hellman)
from ecc import run_demo as ecc_demo
from ecc_attacks import (run_demo as ecc_attacks_demo,
                          run_scaling_experiment_pollard,
                          run_scaling_experiment_smart)
from plots import write_results_txt

# Path shared by all scaling experiments for live CSV output.
CSV_PATH = 'results.csv'


# ── Banner ───────────────────────────────────────────────────────────────────

def print_banner():
    print("\n" + "=" * 60)
    print("  ATTACKS ON PUBLIC-KEY ENCRYPTION SCHEMES")
    print("  Velibor Smilevski")
    print("  University of Primorska, FAMNIT, 2026")
    print("=" * 60)


# ── RSA Attacks ──────────────────────────────────────────────────────────────

def run_all_rsa_attacks(csv_path=None):
    print("\n" + "=" * 60)
    print("  RSA ATTACKS")
    print("=" * 60)

    # --- Wiener's Attack ---
    print("\n>>> WIENER'S ATTACK")
    print("    Precondition: d < (1/3) * N^(1/4)")
    print("    Method: Continued fraction expansion of e/N")
    wiener_demo(prime_bits=256, seed=42)
    wiener_scaling(
        bit_sizes=[128, 256, 512],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )

    # --- Hastad's Broadcast Attack ---
    print("\n>>> HASTAD'S BROADCAST ATTACK")
    print("    Precondition: same M encrypted under e=3 keys")
    print("    Method: CRT reconstruction + integer cube root")
    hastad_demo(e=3, prime_bits=128, seed=42)
    hastad_scaling(
        e=3,
        bit_sizes=[128, 256, 512],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )

    # --- Coppersmith's Short Pad Attack ---
    print("\n>>> COPPERSMITH'S SHORT PAD ATTACK")
    print("    Precondition: e=3, pad length < n/9 bits")
    print("    Method: Resultant + small root finding")
    coppersmith_demo(e=3, prime_bits=128, pad_bits=8, seed=42)
    coppersmith_scaling(
        e=3,
        prime_bits=128,
        pad_bits=8,
        trials=5,
        seed=0,
        csv_path=csv_path,
    )


# ── DH / ElGamal Attacks ─────────────────────────────────────────────────────

def run_all_dh_attacks(csv_path=None):
    print("\n" + "=" * 60)
    print("  DIFFIE-HELLMAN / ELGAMAL ATTACKS")
    print("=" * 60)

    # --- DH and ElGamal background ---
    print("\n>>> DIFFIE-HELLMAN KEY EXCHANGE AND ELGAMAL")
    print("    DH security: discrete logarithm problem in Z_p^*")
    dh_demo(prime_bits=30, seed=42)

    # --- Baby-step Giant-step ---
    print("\n>>> BABY-STEP GIANT-STEP (BSGS)")
    print("    Precondition: group order q known")
    print("    Method: Shanks' meet-in-the-middle, O(sqrt(q))")
    dh_attacks_demo(prime_bits=24, seed=42)
    run_scaling_experiment_bsgs(
        bit_sizes=[16, 20, 24],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )

    # --- Pohlig-Hellman ---
    print("\n>>> POHLIG-HELLMAN ATTACK")
    print("    Precondition: group order q is smooth")
    print("    Method: DLP in prime-power subgroups + CRT")
    run_scaling_experiment_pohlig_hellman(
        bit_sizes=[20, 25, 30],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )


# ── ECC Attacks ───────────────────────────────────────────────────────────────

def run_all_ecc_attacks(csv_path=None):
    print("\n" + "=" * 60)
    print("  ELLIPTIC CURVE CRYPTOGRAPHY ATTACKS")
    print("=" * 60)

    # --- ECC background ---
    # bit_size=14 for demos: anomalous curve search is O(p) per attempt so
    # small p (≈13,000) keeps the demo fast.  Scaling experiment uses 16-20.
    print("\n>>> ELLIPTIC CURVE GROUP")
    print("    Short Weierstrass form: y^2 = x^3 + ax + b (mod p)")
    ecc_demo(bit_size=14, seed=42)

    # --- Pollard's rho ---
    print("\n>>> POLLARD'S RHO FOR ECDLP")
    print("    Precondition: prime-order curve, order n known")
    print("    Method: pseudo-random walk with Floyd cycle detection, O(sqrt(n))")
    ecc_attacks_demo(bit_size=14, seed=42)
    run_scaling_experiment_pollard(
        bit_sizes=[16, 18, 20],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )

    # --- Smart's attack ---
    print("\n>>> SMART'S ATTACK")
    print("    Precondition: anomalous curve (#E(F_p) = p)")
    print("    Method: p-adic lifting + formal group logarithm, O(log^3 p)")
    run_scaling_experiment_smart(
        bit_sizes=[12, 14, 16],
        trials=5,
        seed=0,
        csv_path=csv_path,
    )


# ── Summary table ─────────────────────────────────────────────────────────────

def print_summary():
    print("\n" + "=" * 60)
    print("  SUMMARY OF RESULTS")
    print("=" * 60)

    print("""
  RSA ATTACKS
  -----------
  Wiener's Attack:
    - Recovers private key d when d < (1/3)*N^(1/4)
    - 100% success rate across all tested key sizes
    - Runtime: < 0.001s for 256-1024 bit keys

  Hastad's Broadcast Attack:
    - Recovers plaintext M from e=3 ciphertexts
    - 100% success rate across all tested key sizes
    - Runtime: < 0.03s for 256-1024 bit keys

  Coppersmith's Short Pad Attack:
    - Recovers M from two short-padded ciphertexts
    - Succeeds when pad length < n/9 bits (e=3)
    - Runtime: ~0.05s for 256-bit keys

  DH / ELGAMAL ATTACKS
  --------------------
  Baby-step Giant-step (BSGS):
    - Solves DLP in O(sqrt(q)) time and space
    - 100% success rate for tested subgroup sizes
    - Runtime: < 0.001s for 24-bit subgroup orders

  Pohlig-Hellman:
    - Exploits smooth group order; O(sqrt(B)*log q)
    - 100% success rate when q is smooth
    - Runtime: < 0.001s for 30-bit smooth orders

  ECC ATTACKS
  -----------
  Pollard's Rho (ECDLP):
    - Solves ECDLP in O(sqrt(n)) expected time
    - 100% success rate at tested curve sizes
    - Runtime: < 0.02s for 20-bit curve orders

  Smart's Attack:
    - Recovers k in polynomial time for anomalous curves
    - Applies only when #E(F_p) = p (anomalous condition)
    - Runtime: < 0.001s for tested curve sizes

  KEY LESSON:
    All attacks are defeated by proper parameter choices:
    RSA: use OAEP padding and d >> N^(1/4).
    DH/ElGamal: use a prime-order subgroup of smooth-free order.
    ECC: avoid anomalous curves; use standardised curves (NIST, Brainpool).
""")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Initialise results.csv with header (overwrites any previous run).
    with open(CSV_PATH, 'w', newline='') as f:
        csv.writer(f).writerow([
            'attack', 'parameter_name', 'parameter_value',
            'trials', 'successes', 'avg_time_seconds',
        ])

    print_banner()
    run_all_rsa_attacks(csv_path=CSV_PATH)
    run_all_dh_attacks(csv_path=CSV_PATH)
    run_all_ecc_attacks(csv_path=CSV_PATH)
    print_summary()

    # Regenerate results.txt from the freshly written CSV.
    write_results_txt(CSV_PATH, 'results.txt')
    print("\n[+] results.csv and results.txt written.")
