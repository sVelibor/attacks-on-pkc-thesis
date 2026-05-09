# Coppersmith's Short Pad Attack on RSA
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999
#
# Implementation scope
# --------------------
# This implementation uses the resultant-based reduction of Coppersmith's
# short pad attack and recovers the small root delta = r2 - r1 by direct
# search within the bound [-2^pad_bits, 2^pad_bits]. This is feasible for
# the small pad sizes used in the thesis experiments. A full LLL-based
# small-root finder, which would extend the attack to larger pads, is
# left as future work (see thesis Section 9.4).

from math import gcd
from utils import modinv, generate_prime


# ── Polynomial GCD over Z_N ──────────────────────────────────────────────────

def poly_gcd_mod(f, g, N):
    """Compute GCD of two polynomials over Z_N using the Euclidean algorithm.

    Polynomials are represented as coefficient lists [a_0, a_1, ..., a_d]
    where the index equals the degree of the corresponding term.

    Returns the GCD polynomial (monic, normalized) or None if it fails.
    """
    def poly_mod(p, N):
        return [c % N for c in p]

    def poly_degree(p):
        p = list(p)
        while len(p) > 1 and p[-1] % N == 0:
            p.pop()
        return len(p) - 1

    def poly_div_mod(a, b, N):
        """Polynomial division a / b over Z_N."""
        a = list(a)
        b = list(b)
        deg_a = poly_degree(a)
        deg_b = poly_degree(b)

        if deg_a < deg_b:
            return [], a

        result = [0] * (deg_a - deg_b + 1)
        for i in range(deg_a - deg_b, -1, -1):
            if deg_a < deg_b:
                break
            lead_b = b[deg_b] % N
            if lead_b == 0:
                return None, None
            try:
                inv_lead = modinv(lead_b, N)
            except ValueError:
                return None, None
            coeff = (a[deg_a] * inv_lead) % N
            result[i] = coeff
            for j in range(deg_b + 1):
                a[i + j] = (a[i + j] - coeff * b[j]) % N
            deg_a = poly_degree(a)

        remainder = poly_mod(a[:deg_b], N)
        return result, remainder

    f = poly_mod(f, N)
    g = poly_mod(g, N)

    for _ in range(100):
        if poly_degree(g) == 0:
            break
        _, r = poly_div_mod(f, g, N)
        # Division fails when the leading coefficient is not invertible mod N,
        # which happens when it shares a factor with the composite modulus N.
        # In that case the GCD cannot be computed and None is returned.
        if r is None:
            return None
        f, g = g, poly_mod(r, N)

    # Normalize to monic
    deg = poly_degree(f)
    lead = f[deg] % N
    if lead == 0:
        return None
    try:
        inv_lead = modinv(lead, N)
    except ValueError:
        return None

    f = [(c * inv_lead) % N for c in f]
    return f


# ── Short pad attack ─────────────────────────────────────────────────────────

def coppersmith_short_pad_sympy(N, e, C1, C2, pad_bits):
    """Coppersmith's short pad attack via polynomial resultant + direct search.

    Given two ciphertexts C1 = (M + r1)^e mod N and C2 = (M + r2)^e mod N
    with short random pads r1, r2, recover M1 = M + r1.

    The attack proceeds in two stages:
      (1) Form g1(x) = x^e - C1 and g2(x, y) = (x + y)^e - C2 and compute
          h(y) = res_x(g1, g2). The pad difference delta = r2 - r1 is a
          small root of h(y) mod N.
      (2) Find delta by direct search over [-2^pad_bits, 2^pad_bits], then
          recover M1 = M + r1 via the Franklin-Reiter step: a polynomial
          GCD of x^e - C1 and (x + delta)^e - C2 modulo N yields a linear
          factor (x - M1).

    Note: stage (2) uses direct search rather than LLL-based small-root
    finding. This is feasible for the small pad sizes used in the thesis
    experiments but does not extend to larger parameters; a full LLL
    implementation is future work.

    Parameters:
        N        : RSA modulus
        e        : public exponent
        C1, C2   : two ciphertexts of the same message with different pads
        pad_bits : number of bits in the pad (search bound)

    Returns:
        M1 = M + r1 if attack succeeds, None otherwise.
    """
    from sympy import symbols, resultant, Poly, expand

    x, y = symbols('x y')

    # Stage 1: form polynomials with M1 = M + r1 as a common root mod N
    g1 = x**e - C1
    g2 = (x + y)**e - C2

    print("    [*] Computing resultant res_x(g1, g2)...")
    h = resultant(g1, g2, x)
    print(f"    [*] Resultant computed, degree in y = {e**2}")

    h_poly = Poly(h, y)
    coeffs_raw = h_poly.all_coeffs()

    # Stage 2: direct search for delta = r2 - r1 in [-2^pad_bits, 2^pad_bits]
    print(f"    [*] Searching for small root delta (pad_bits={pad_bits})...")
    bound = 2 ** pad_bits
    direct_roots = []
    for delta in range(-bound, bound + 1):
        val = sum(int(coeffs_raw[-(i + 1)]) * pow(delta, i, N)
                  for i in range(len(coeffs_raw))) % N
        if val == 0:
            direct_roots.append(delta)

    if not direct_roots:
        print("    [-] No small root found for delta")
        return None

    print(f"    [+] Found candidate delta values: {direct_roots}")

    # Franklin-Reiter recovery: for each candidate delta, recover M1 via
    # polynomial GCD of (x^e - C1) and ((x + delta)^e - C2) modulo N.
    for delta in direct_roots:
        if delta == 0:
            continue

        f1 = [(-C1) % N] + [0] * (e - 1) + [1]  # x^e - C1

        sp = Poly(expand((x + delta)**e - C2), x)
        f2_coeffs_high = [int(c) % N for c in sp.all_coeffs()]
        f2 = list(reversed(f2_coeffs_high))  # [a_0, a_1, ..., a_e]

        result = poly_gcd_mod(f1, f2, N)

        if result is not None and len(result) >= 2:
            deg = len(result) - 1
            if deg == 1:
                # Linear factor (x - root): root = -result[0] / result[1]
                try:
                    M1 = (-result[0] * modinv(result[1], N)) % N
                    if pow(M1, e, N) == C1:
                        print(f"    [+] Recovered M1 = M + r1 = {M1}")
                        return M1
                except ValueError:
                    continue

    return None


# ── RSA key generation ───────────────────────────────────────────────────────

def generate_rsa_keypair_fixed_e(e, prime_bits, rng):
    """Generate an RSA key pair with a fixed public exponent e."""
    while True:
        p = generate_prime(prime_bits, rng)
        q = generate_prime(prime_bits, rng)
        if p == q:
            continue
        N = p * q
        phi = (p - 1) * (q - 1)
        if gcd(e, phi) != 1:
            continue
        d = modinv(e, phi)
        return N, e, d, p, q


# ── Demo ─────────────────────────────────────────────────────────────────────

def run_demo(e=3, prime_bits=128, pad_bits=8, seed=42):
    """Run Coppersmith's short pad attack demonstration.

    Uses small parameters for a clear, fast demonstration.
    """
    import random
    import time

    rng = random.Random(seed)

    print("=" * 60)
    print("COPPERSMITH'S SHORT PAD ATTACK DEMONSTRATION")
    print("=" * 60)

    print(f"\n[*] Parameters: e={e}, modulus={2*prime_bits} bits, "
          f"pad={pad_bits} bits")
    print(f"    Condition: pad_bits < modulus_bits / {e**2} = "
          f"{2*prime_bits // e**2} bits")
    print(f"    Condition met: {pad_bits < 2*prime_bits // e**2}")

    # Generate RSA key pair
    print(f"\n[*] Generating RSA key pair...")
    N, e_key, d, p, q = generate_rsa_keypair_fixed_e(e, prime_bits, rng)
    print(f"    N = {N}")
    print(f"    e = {e_key}")

    # Choose message and pads
    max_M = N >> (pad_bits + 4)
    M = rng.randint(2, max(2, max_M))
    r1 = rng.randint(0, 2**pad_bits - 1)
    r2 = rng.randint(0, 2**pad_bits - 1)

    while r2 == r1:
        r2 = rng.randint(0, 2**pad_bits - 1)

    M1 = M + r1
    M2 = M + r2

    print(f"\n[*] Message M  = {M}")
    print(f"    Pad r1     = {r1}")
    print(f"    Pad r2     = {r2}")
    print(f"    delta      = r2 - r1 = {r2 - r1}")
    print(f"    M1 = M+r1  = {M1}")
    print(f"    M2 = M+r2  = {M2}")

    # Encrypt
    C1 = pow(M1, e, N)
    C2 = pow(M2, e, N)
    print(f"\n[*] C1 = M1^e mod N = {C1}")
    print(f"    C2 = M2^e mod N = {C2}")

    # Run attack
    print(f"\n[*] Running Coppersmith's short pad attack...")
    start = time.time()
    M1_recovered = coppersmith_short_pad_sympy(N, e, C1, C2, pad_bits)
    elapsed = time.time() - start

    if M1_recovered is not None:
        print(f"\n[+] Attack succeeded in {elapsed:.4f} seconds")
        print(f"    Recovered M1 = {M1_recovered}")
        print(f"    Correct M1   : {M1_recovered == M1}")
    else:
        print(f"\n[-] Attack failed in {elapsed:.4f} seconds")


# ── Scaling experiment ────────────────────────────────────────────────────────

def run_scaling_experiment(e=3, prime_bits=128, pad_bits=8, trials=5,
                           seed=0, csv_path=None):
    """Run Coppersmith's attack for a fixed parameter set across multiple trials.

    Coppersmith's attack is currently practical only for small pad sizes
    (direct root search), so a single operating point is measured rather
    than a sweep.  The modulus size is 2 * prime_bits bits.

    Parameters:
        e          : public exponent
        prime_bits : half-size of the RSA modulus in bits
        pad_bits   : pad length in bits
        trials     : number of independent trials
        seed       : base random seed
        csv_path   : if given, append one CSV row to this file
    """
    import csv
    import random
    import time

    key_bits = 2 * prime_bits

    print("\n" + "=" * 60)
    print("COPPERSMITH'S SHORT PAD ATTACK - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'Key bits':>10} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 46)

    successes  = 0
    total_time = 0.0

    for trial in range(trials):
        rng = random.Random(seed + trial)
        N, e_key, d, p, q = generate_rsa_keypair_fixed_e(e, prime_bits, rng)
        max_M = N >> (pad_bits + 4)
        M  = rng.randint(2, max(2, max_M))
        r1 = rng.randint(0, 2 ** pad_bits - 1)
        r2 = rng.randint(0, 2 ** pad_bits - 1)
        while r2 == r1:
            r2 = rng.randint(0, 2 ** pad_bits - 1)
        M1 = M + r1
        M2 = M + r2
        C1 = pow(M1, e_key, N)
        C2 = pow(M2, e_key, N)

        start   = time.time()
        M1_rec  = coppersmith_short_pad_sympy(N, e_key, C1, C2, pad_bits)
        elapsed = time.time() - start
        total_time += elapsed

        if M1_rec == M1:
            successes += 1

    avg_time = total_time / trials
    print(f"{key_bits:>10} {trials:>7} {successes:>9} {avg_time:>14.6f}")

    if csv_path is not None:
        with open(csv_path, 'a', newline='') as f:
            csv.writer(f).writerow(
                ['coppersmith', 'key_bits', key_bits, trials, successes, avg_time])

    return [(key_bits, trials, successes, avg_time)]


# ── Entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(
        e=3,
        prime_bits=128,
        pad_bits=8,
        seed=42,
    )