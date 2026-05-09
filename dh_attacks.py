# dh_attacks.py
# Discrete logarithm attacks on Diffie-Hellman / ElGamal
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999
#           and Menezes, van Oorschot, Vanstone, "Handbook of Applied Cryptography"
#
# Attacks implemented:
#   1. Baby-step Giant-step (Shanks, 1971) -- O(sqrt(q)) time and space
#   2. Pohlig-Hellman (1978) -- O(sqrt(B) * log q) for smooth-order groups

from math import isqrt, gcd
from utils import modinv, crt, generate_prime


# ── Baby-step Giant-step ──────────────────────────────────────────────────────

def bsgs(g, h, p, q):
    """Baby-step Giant-step algorithm for the discrete logarithm problem.

    Given g, h in Z_p^* and subgroup order q, find x such that
    g^x ≡ h (mod p) for 0 <= x < q.

    Algorithm (Shanks):
      Let m = ceil(sqrt(q)).
      Baby steps: store {g^j mod p : j = 0, ..., m-1} in a hash table.
      Giant steps: compute h * (g^(-m))^i mod p for i = 0, 1, ..., m
        and look up in the table. A collision gives x = j + i*m.

    Parameters:
        g : generator (integer)
        h : target   (integer, h = g^x mod p)
        p : prime modulus
        q : order of g in Z_p^*

    Returns:
        x such that g^x ≡ h (mod p), or None if not found.
    """
    m = isqrt(q) + 1

    # Baby steps: table[g^j mod p] = j
    baby = {}
    gj = 1
    for j in range(m):
        baby[gj] = j
        gj = gj * g % p

    # Giant-step factor: g^(-m) mod p
    g_inv_m = modinv(pow(g, m, p), p)

    # Giant steps
    gamma = h % p
    for i in range(m + 1):
        if gamma in baby:
            x = (baby[gamma] + i * m) % q
            if pow(g, x, p) == h % p:
                return x
        gamma = gamma * g_inv_m % p

    return None


# ── Pohlig-Hellman ─────────────────────────────────────────────────────────────

def pohlig_hellman(g, h, p, q, factorization):
    """Pohlig-Hellman algorithm for the discrete logarithm problem.

    Exploits a smooth group order q = prod(l_i^e_i) to reduce the DLP to
    smaller subproblems of order l_i^e_i, then combines via CRT.

    For each prime power l^e in the factorization:
      * Work in the subgroup of order l by raising elements to q/l^e.
      * Recover the discrete log digit-by-digit in base l using BSGS
        on a subgroup of order l.
      * Combine the residues x_i ≡ x (mod l_i^e_i) via CRT.

    Parameters:
        g             : generator
        h             : target (h = g^x mod p)
        p             : prime modulus
        q             : order of g (integer)
        factorization : list of (prime, exponent) pairs for q

    Returns:
        x such that g^x ≡ h (mod p), or None on failure.
    """
    residues = []
    moduli   = []

    for l, e in factorization:
        l_e = l ** e

        # Work in the subgroup of order l: precompute g^(q/l) which has order l
        q_over_l  = q // l
        g_l       = pow(g, q_over_l, p)   # order l

        x_k = 0  # accumulate digits of x mod l^e
        for k in range(e):
            # Current "base" element contribution
            # h_k = (h * g^(-x_k))^(q / l^(k+1)) has order dividing l
            l_k1      = l ** (k + 1)
            exponent  = q // l_k1
            h_k_base  = h * modinv(pow(g, x_k, p), p) % p
            h_k       = pow(h_k_base, exponent, p)

            # DLP in subgroup of order l: find d_k such that g_l^d_k = h_k
            d_k = bsgs(g_l, h_k, p, l)
            if d_k is None:
                raise RuntimeError(
                    f"BSGS failed in subgroup of order {l}; "
                    f"this indicates a bug in subgroup setup."
                )

            x_k = (x_k + d_k * (l ** k)) % l_e

        residues.append(x_k)
        moduli.append(l_e)

    # CRT to combine residues
    try:
        x = crt(residues, moduli)
    except Exception:
        return None

    # Verify
    if pow(g, x, p) == h % p:
        return x

    return None


# ── Demo ─────────────────────────────────────────────────────────────────────

def run_demo(prime_bits=30, seed=42):
    """Demonstrate BSGS and Pohlig-Hellman discrete logarithm attacks."""
    import random
    import time
    from dh import (generate_dh_params, generate_smooth_dh_params,
                    dh_keygen, elgamal_keygen, elgamal_encrypt)

    rng = random.Random(seed)

    print("=" * 60)
    print("DISCRETE LOGARITHM ATTACKS DEMONSTRATION")
    print("=" * 60)

    # --- Baby-step Giant-step on standard group ---
    print(f"\n[*] Generating DH parameters ({prime_bits}-bit subgroup order)...")
    p, q, g = generate_dh_params(prime_bits, rng)
    print(f"    p = {p}, q = {q}, g = {g}")

    x_real, h = dh_keygen(p, q, g, rng)
    print(f"\n[*] Target: h = g^x mod p,  x = {x_real}  (secret)")

    print(f"\n[*] Running Baby-step Giant-step...")
    start = time.time()
    x_bsgs = bsgs(g, h, p, q)
    elapsed = time.time() - start

    if x_bsgs is not None:
        print(f"    [+] BSGS recovered x = {x_bsgs} in {elapsed:.6f} s")
        print(f"    Correct: {x_bsgs == x_real or pow(g, x_bsgs, p) == h}")
    else:
        print(f"    [-] BSGS failed")

    # --- Pohlig-Hellman on smooth-order group ---
    print(f"\n[*] Generating smooth-order DH parameters...")
    p_s, q_s, g_s, fact_s = generate_smooth_dh_params(prime_bits, rng)
    print(f"    p = {p_s}")
    print(f"    q = {q_s}  (smooth order)")
    print(f"    q = {' * '.join(f'{pr}^{ex}' for pr, ex in fact_s)}")

    x_s_real = rng.randint(2, q_s - 1)
    h_s = pow(g_s, x_s_real, p_s)
    print(f"\n[*] Target: h = g^x mod p,  x = {x_s_real}  (secret)")

    print(f"\n[*] Running Pohlig-Hellman...")
    start = time.time()
    x_ph = pohlig_hellman(g_s, h_s, p_s, q_s, fact_s)
    elapsed = time.time() - start

    if x_ph is not None:
        print(f"    [+] Pohlig-Hellman recovered x = {x_ph} in {elapsed:.6f} s")
        print(f"    Correct: {pow(g_s, x_ph, p_s) == h_s}")
    else:
        print(f"    [-] Pohlig-Hellman failed")

    # --- ElGamal key recovery via BSGS ---
    print(f"\n[*] ElGamal private-key recovery via BSGS...")
    from dh import elgamal_decrypt
    x_eg, h_eg = elgamal_keygen(p, q, g, rng)
    M = rng.randint(2, p - 2)
    c1, c2 = elgamal_encrypt(p, g, h_eg, M, rng, q)

    start = time.time()
    x_eg_rec = bsgs(g, h_eg, p, q)
    elapsed = time.time() - start

    if x_eg_rec is not None:
        M_dec = elgamal_decrypt(p, x_eg_rec, c1, c2)
        print(f"    [+] Recovered private key in {elapsed:.6f} s")
        print(f"    Decrypted M = {M_dec},  Correct: {M_dec == M}")
    else:
        print(f"    [-] Key recovery failed")


# ── Scaling experiments ────────────────────────────────────────────────────────

def run_scaling_experiment_bsgs(bit_sizes=None, trials=5, seed=0, csv_path=None):
    """Run BSGS across different subgroup orders.

    Parameters:
        bit_sizes : list of subgroup order sizes in bits
        trials    : independent trials per size
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time
    from dh import generate_dh_params, dh_keygen

    if bit_sizes is None:
        bit_sizes = [16, 20, 24]

    print("\n" + "=" * 60)
    print("BABY-STEP GIANT-STEP - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'q bits':>8} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 44)

    results = []

    for bits in bit_sizes:
        successes   = 0
        total_time  = 0.0

        for trial in range(trials):
            rng = random.Random(seed + trial)
            try:
                p, q, g = generate_dh_params(bits, rng)
                x_real, h = dh_keygen(p, q, g, rng)

                start      = time.time()
                x_rec      = bsgs(g, h, p, q)
                elapsed    = time.time() - start
                total_time += elapsed

                if x_rec is not None and pow(g, x_rec, p) == h:
                    successes += 1

            except Exception as ex:
                print(f"    [!] Trial {trial} failed: {ex}")

        avg_time = total_time / trials
        print(f"{bits:>8} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['bsgs', 'q_bits', bits, trials, successes, avg_time])

    return results


def run_scaling_experiment_pohlig_hellman(bit_sizes=None, trials=5, seed=0,
                                          csv_path=None):
    """Run Pohlig-Hellman across different (smooth) group orders.

    Parameters:
        bit_sizes : list of smooth group order sizes in bits
        trials    : independent trials per size
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time
    from dh import generate_smooth_dh_params

    if bit_sizes is None:
        bit_sizes = [20, 25, 30]

    print("\n" + "=" * 60)
    print("POHLIG-HELLMAN - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'q bits':>8} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 44)

    results = []

    for bits in bit_sizes:
        successes   = 0
        total_time  = 0.0

        for trial in range(trials):
            rng = random.Random(seed + trial)
            try:
                p, q, g, fact = generate_smooth_dh_params(bits, rng)
                x_real = rng.randint(2, q - 1)
                h      = pow(g, x_real, p)

                start      = time.time()
                x_rec      = pohlig_hellman(g, h, p, q, fact)
                elapsed    = time.time() - start
                total_time += elapsed

                if x_rec is not None and pow(g, x_rec, p) == h:
                    successes += 1

            except Exception as ex:
                print(f"    [!] Trial {trial} failed: {ex}")

        avg_time = total_time / trials
        print(f"{bits:>8} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['pohlig_hellman', 'q_bits', bits, trials, successes, avg_time])

    return results


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(prime_bits=24, seed=42)

    run_scaling_experiment_bsgs(
        bit_sizes=[16, 20, 24],
        trials=5,
        seed=0,
    )

    run_scaling_experiment_pohlig_hellman(
        bit_sizes=[20, 25, 30],
        trials=5,
        seed=0,
    )
