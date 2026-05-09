# hastad.py
# Hastad's Broadcast Attack on RSA with small public exponent e=3
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999

from math import gcd
from utils import modinv, crt, integer_nth_root, generate_prime


# ── Core attack ─────────────────────────────────────────────────────────────

def hastad_attack(ciphertexts, moduli, e):
    """Hastad's broadcast attack.

    Given e ciphertexts C_i = M^e mod N_i, all encrypting the same
    plaintext M with the same public exponent e, recover M.

    Parameters:
        ciphertexts : list of e ciphertexts [C_1, ..., C_e]
        moduli      : list of e moduli [N_1, ..., N_e]
        e           : public exponent (typically 3)

    Returns:
        M if the attack succeeds, None otherwise.
    """
    assert len(ciphertexts) >= e, "Need at least e ciphertexts"
    assert len(moduli) >= e, "Need at least e moduli"

    # Check moduli are pairwise coprime
    for i in range(e):
        for j in range(i + 1, e):
            if gcd(moduli[i], moduli[j]) != 1:
                raise ValueError(
                    f"Moduli N_{i} and N_{j} are not coprime. "
                    f"Factorization exposed: gcd = {gcd(moduli[i], moduli[j])}"
                )

    # Step 1: CRT to reconstruct M^e mod (N_1 * N_2 * ... * N_e)
    C_prime = crt(ciphertexts[:e], moduli[:e])

    # Step 2: Integer e-th root extraction
    # Since M < N_i for all i, we have M^e < N_1*N_2*...*N_e
    # so C_prime = M^e exactly as an integer
    M = integer_nth_root(C_prime, e)

    # Step 3: Verify
    if M ** e == C_prime:
        return M

    return None


# ── Key generation for broadcast scenario ───────────────────────────────────

def generate_broadcast_setup(e, prime_bits, rng):
    """Generate e independent RSA key pairs all with public exponent e.

    Returns:
        moduli     : list of e moduli N_i
        public_keys: list of (N_i, e) pairs
        private_keys: list of d_i values
    """
    from sympy import isprime

    moduli = []
    private_keys = []

    while len(moduli) < e:
        # Generate two distinct primes
        p = generate_prime(prime_bits, rng)
        q = generate_prime(prime_bits, rng)

        if p == q:
            continue

        N = p * q
        phi = (p - 1) * (q - 1)

        # e must be coprime with phi
        if gcd(e, phi) != 1:
            continue

        # Check this modulus is coprime with all previous moduli
        if any(gcd(N, Nj) != 1 for Nj in moduli):
            continue

        d = modinv(e, phi)
        moduli.append(N)
        private_keys.append(d)

    public_keys = [(N, e) for N in moduli]
    return moduli, public_keys, private_keys


# ── Demo and verification ────────────────────────────────────────────────────

def run_demo(e=3, prime_bits=128, seed=42):
    """Run a full demonstration of Hastad's broadcast attack."""
    import random
    import time

    rng = random.Random(seed)

    print("=" * 60)
    print("HASTAD'S BROADCAST ATTACK DEMONSTRATION")
    print("=" * 60)

    # Generate e independent RSA key pairs
    print(f"\n[*] Generating {e} independent RSA key pairs "
          f"(e={e}, {2*prime_bits}-bit moduli)...")

    moduli, public_keys, private_keys = generate_broadcast_setup(
        e, prime_bits, rng
    )

    for i, (N, _) in enumerate(public_keys):
        print(f"    N_{i+1} = {N}")

    # Choose a plaintext M smaller than all moduli
    M = rng.randint(2, min(moduli) - 1)
    print(f"\n[*] Plaintext M = {M}")
    print(f"    (M < all N_i: {all(M < N for N in moduli)})")

    # Encrypt M under each public key
    ciphertexts = [pow(M, e, N) for N in moduli]
    print(f"\n[*] Ciphertexts:")
    for i, C in enumerate(ciphertexts):
        print(f"    C_{i+1} = {C}")

    # Run attack
    print(f"\n[*] Running Hastad's broadcast attack...")
    start = time.time()
    M_recovered = hastad_attack(ciphertexts, moduli, e)
    elapsed = time.time() - start

    # Report result
    if M_recovered is not None:
        print(f"\n[+] Attack succeeded in {elapsed:.6f} seconds")
        print(f"    Recovered M  = {M_recovered}")
        print(f"    Correct      : {M_recovered == M}")
    else:
        print(f"\n[-] Attack failed")

    return M_recovered


# ── Scaling experiment ───────────────────────────────────────────────────────

def run_scaling_experiment(e=3, bit_sizes=None, trials=5, seed=0, csv_path=None):
    """Run Hastad's attack across different modulus sizes.

    Parameters:
        e         : public exponent (default 3)
        bit_sizes : list of prime half-sizes; modulus size = 2 * bits
        trials    : number of independent trials per size
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time

    if bit_sizes is None:
        bit_sizes = [128, 256, 512]

    print("\n" + "=" * 60)
    print("HASTAD'S BROADCAST ATTACK - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'Bits':>6} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 42)

    results = []

    for bits in bit_sizes:
        successes = 0
        total_time = 0.0

        for trial in range(trials):
            rng = random.Random(seed + trial)
            try:
                moduli, _, _ = generate_broadcast_setup(e, bits, rng)
                M = rng.randint(2, min(moduli) - 1)
                ciphertexts = [pow(M, e, N) for N in moduli]

                start = time.time()
                M_recovered = hastad_attack(ciphertexts, moduli, e)
                elapsed = time.time() - start
                total_time += elapsed

                if M_recovered == M:
                    successes += 1

            except Exception as ex:
                print(f"    [!] Trial {trial} failed: {ex}")

        avg_time = total_time / trials
        key_bits = 2 * bits
        print(f"{key_bits:>6} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((key_bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['hastad', 'key_bits', key_bits, trials, successes, avg_time])

    return results


# ── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Run demonstration with e=3
    run_demo(e=3, prime_bits=256, seed=42)

    # Run scaling experiment
    run_scaling_experiment(
        e=3,
        bit_sizes=[128, 256, 512],
        trials=5,
        seed=0
    )