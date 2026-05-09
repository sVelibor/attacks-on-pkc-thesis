# wiener.py
# Wiener's Attack on RSA with small private exponent
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999

from math import isqrt, gcd
from utils import modinv, generate_prime


# ── Continued fraction tools ────────────────────────────────────────────────

def convergents(n, d):
    """Generate convergents (numerator, denominator) of n/d.
    Uses the standard recurrence:
        h_i = a_i * h_{i-1} + h_{i-2}
        k_i = a_i * k_{i-1} + k_{i-2}
    with h_{-2}=0, h_{-1}=1, k_{-2}=1, k_{-1}=0.
    """
    h2, h1 = 0, 1
    k2, k1 = 1, 0

    while d:
        a = n // d
        h = a * h1 + h2
        k = a * k1 + k2
        yield h, k
        h2, h1 = h1, h
        k2, k1 = k1, k
        n, d = d, n % d


# ── Core attack ─────────────────────────────────────────────────────────────

def wiener_attack(e, N):
    """Wiener's attack: recover private exponent d from public key (N, e)
    when d < (1/3) * N^(1/4).

    Returns d if the attack succeeds, None otherwise.
    """
    for k, d in convergents(e, N):
        # Skip trivial cases
        if k == 0 or d == 0:
            continue

        # Check if (e*d - 1) is divisible by k
        if (e * d - 1) % k != 0:
            continue

        # Candidate phi(N) = (e*d - 1) / k
        phi = (e * d - 1) // k

        # phi must be positive and less than N
        if phi <= 0 or phi >= N:
            continue

        # Recover p and q by solving x^2 - (N - phi + 1)x + N = 0
        # From phi = (p-1)(q-1) = N - p - q + 1 => p + q = N - phi + 1
        b = N - phi + 1
        discriminant = b * b - 4 * N

        if discriminant < 0:
            continue

        sqrt_disc = isqrt(discriminant)

        # Check discriminant is a perfect square
        if sqrt_disc * sqrt_disc != discriminant:
            continue

        p = (b + sqrt_disc) // 2
        q = (b - sqrt_disc) // 2

        # Verify factorization
        if p > 1 and q > 1 and p * q == N:
            return d

    return None


# ── Key generation with weak d ──────────────────────────────────────────────

def generate_weak_rsa_keypair(prime_bits, rng):
    """Generate an RSA key pair with a deliberately small private exponent d,
    satisfying d < (1/3) * N^(1/4) to make Wiener's attack applicable.

    Returns (N, e, d, p, q).
    """
    from sympy import isprime

    while True:
        # Generate two distinct primes with q < p < 2q
        p = generate_prime(prime_bits, rng)
        q = generate_prime(prime_bits, rng)

        if p == q:
            continue

        if p < q:
            p, q = q, p

        if p >= 2 * q:
            continue

        N = p * q
        phi = (p - 1) * (q - 1)

        # Choose d well within Wiener's bound: use N^(1/4) / 12
        N_quarter = isqrt(isqrt(N))
        max_d = N_quarter // 12
        min_d = max(2, max_d // 2)

        if max_d < 3:
            continue

        # Pick a random odd d in valid range
        d = rng.randint(min_d, max_d)
        if d % 2 == 0:
            d -= 1
        if d < 2:
            continue

        if gcd(d, phi) != 1:
            continue

        e = modinv(d, phi)

        if e <= 1 or e >= phi:
            continue

        return N, e, d, p, q


# ── Demo and verification ───────────────────────────────────────────────────

def run_demo(prime_bits=256, seed=42):
    """Run a full demonstration of Wiener's attack."""
    import random
    import time

    rng = random.Random(seed)

    print("=" * 60)
    print("WIENER'S ATTACK DEMONSTRATION")
    print("=" * 60)

    print(f"\n[*] Generating weak RSA key pair ({2*prime_bits}-bit modulus)...")
    N, e, d_real, p, q = generate_weak_rsa_keypair(prime_bits, rng)

    print(f"    N  = {N}")
    print(f"    e  = {e}")
    print(f"    d  = {d_real}  (secret)")
    print(f"    p  = {p}")
    print(f"    q  = {q}")
    print(f"\n    Wiener bound : d < N^(1/4)/3  = {isqrt(isqrt(N))//3}")
    print(f"    Actual d     : {d_real}")
    print(f"    Condition met: {d_real < isqrt(isqrt(N))//3}")

    print(f"\n[*] Running Wiener's attack on public key (N, e)...")
    start = time.time()
    d_recovered = wiener_attack(e, N)
    elapsed = time.time() - start

    if d_recovered is not None:
        print(f"\n[+] Attack succeeded in {elapsed:.6f} seconds")
        print(f"    Recovered d  = {d_recovered}")
        print(f"    Correct      : {d_recovered == d_real}")

        # Verify decryption works
        M = rng.randint(2, N - 1)
        C = pow(M, e, N)
        M_dec = pow(C, d_recovered, N)
        print(f"\n[*] Encryption/decryption test:")
        print(f"    Original message   : {M}")
        print(f"    Decrypted message  : {M_dec}")
        print(f"    Decryption correct : {M == M_dec}")
    else:
        print(f"\n[-] Attack failed")

    return d_recovered


# ── Scaling experiment ──────────────────────────────────────────────────────

def run_scaling_experiment(bit_sizes, trials=5, seed=0, csv_path=None):
    """Run Wiener's attack across different key sizes.

    Parameters:
        bit_sizes : list of prime half-sizes; modulus size = 2 * bits
        trials    : number of independent trials per size
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time

    print("\n" + "=" * 60)
    print("WIENER'S ATTACK - SCALING EXPERIMENT")
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
                N, e, d_real, p, q = generate_weak_rsa_keypair(bits, rng)
                start = time.time()
                d_recovered = wiener_attack(e, N)
                elapsed = time.time() - start
                total_time += elapsed
                if d_recovered == d_real:
                    successes += 1
            except Exception:
                pass

        avg_time = total_time / trials
        key_bits = 2 * bits
        print(f"{key_bits:>6} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((key_bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['wiener', 'key_bits', key_bits, trials, successes, avg_time])

    return results


# ── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(prime_bits=256, seed=42)

    run_scaling_experiment(
        bit_sizes=[128, 256, 512],
        trials=5,
        seed=0
    )