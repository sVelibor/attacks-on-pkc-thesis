# ecc_attacks.py
# Elliptic Curve Discrete Logarithm attacks
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999
#           Smart, "The discrete logarithm problem on elliptic curves of trace one", 1999
#
# Attacks implemented:
#   1. Pollard's rho for ECDLP -- O(sqrt(N)) expected, Floyd cycle detection
#   2. Smart's attack           -- polynomial time for anomalous curves (#E=p)

from utils import modinv


# ── Pollard's rho for ECDLP ───────────────────────────────────────────────────

def pollard_rho_ecdlp(curve, G, Q, n, max_iter=500_000, seed=0):
    """Pollard's rho algorithm for the ECDLP.

    Given a prime-order curve, generator G of order n, and point Q = kG,
    find k using a pseudo-random walk with Floyd's cycle detection.

    Partition: the walk is split into three sets S1, S2, S3 based on the
    x-coordinate of the current point mod 3.  Each step:
      S1: X -> X + G,     a -> a+1,   b -> b
      S2: X -> 2X,        a -> 2a,    b -> 2b
      S3: X -> X + Q,     a -> a,     b -> b+1

    Collision (X_i = X_{2i}) gives  a_i + b_i*k ≡ a_{2i} + b_{2i}*k (mod n)
    => k = (a_i - a_{2i}) * modinv(b_{2i} - b_i, n) (mod n).

    Parameters:
        curve   : EllipticCurve instance
        G       : generator point, order n
        Q       : target point (Q = kG)
        n       : prime order of G
        max_iter: maximum walk length before restarting
        seed    : starting seed for random restarts

    Returns:
        k such that kG = Q, or None on failure.
    """
    import random

    rng = random.Random(seed)

    def step(X, a, b):
        if X is None:
            s = 0
        else:
            s = X[0] % 3
        if s == 0:
            return curve.add(X, G), (a + 1) % n, b
        elif s == 1:
            return curve.double(X), (2 * a) % n, (2 * b) % n
        else:
            return curve.add(X, Q), a, (b + 1) % n

    # Multiple restarts with varied starting points
    for restart in range(20):
        a0 = rng.randint(1, n - 1)
        b0 = rng.randint(1, n - 1)
        X0 = curve.add(curve.scalar_mult(a0, G), curve.scalar_mult(b0, Q))

        # Floyd's tortoise and hare
        tort_X, tort_a, tort_b = X0, a0, b0
        hare_X, hare_a, hare_b = X0, a0, b0

        for _ in range(max_iter):
            tort_X, tort_a, tort_b = step(tort_X, tort_a, tort_b)
            hare_X, hare_a, hare_b = step(hare_X, hare_a, hare_b)
            hare_X, hare_a, hare_b = step(hare_X, hare_a, hare_b)

            if tort_X == hare_X:
                # Collision found: a_t + b_t*k = a_h + b_h*k (mod n)
                db = (hare_b - tort_b) % n
                da = (tort_a - hare_a) % n
                if db == 0:
                    break  # degenerate collision; restart
                try:
                    k = da * modinv(db, n) % n
                except ValueError:
                    break
                if curve.scalar_mult(k, G) == Q:
                    return k
                break  # wrong k; restart

    return None


# ── Smart's attack ────────────────────────────────────────────────────────────

def smart_attack(curve, G, Q):
    """Smart's attack on the ECDLP for anomalous curves (#E(F_p) = p).

    Lifts G and Q from E(F_p) to E(Z/p^2Z) via Hensel's lemma, then
    multiplies by p using projective coordinates (no modular inversion
    required), obtaining points in the formal-group kernel K(p).
    The formal logarithm t = -X/Y then reduces the ECDLP to division in Z/pZ.

    Algorithm:
      1. Hensel-lift y-coordinates of G, Q to Z/p^2Z.
      2. Embed in homogeneous projective: G' = (x_G : y_G_lift : 1).
      3. Compute p*G' and p*Q' via projective double-and-add mod p^2.
         Projective formulas (no inversion) avoid the singular denominator
         that arises in affine arithmetic when x_1 ≡ x_2 (mod p) during
         the scalar multiplication toward the formal-group kernel.
      4. For result (X:Y:Z), the formal-group parameter is t = -X/Y mod p^2
         (invariant under projective scaling).  Since p*G' ∈ K(p), t is
         divisible by p; the reduced log is t/p mod p.
      5. k = (t(p*Q')/p) * modinv(t(p*G')/p, p)  mod p.

    Parameters:
        curve : EllipticCurve instance (must be anomalous: #E(F_p) = p)
        G     : generator point of order p
        Q     : target point (Q = kG)

    Returns:
        k such that kG = Q (mod p), or None on failure.
    """
    p  = curve.p
    a  = curve.a
    b  = curve.b
    pp = p * p

    # ── Hensel lift ───────────────────────────────────────────────────────────

    def lift_y(x_int, y_int):
        """Lift y from F_p to Z/p^2Z for fixed x_int.

        Newton step: c = (y^2 - rhs)/p mod p, then y_lift = y - c*inv(2y)*p.
        """
        rhs   = (pow(x_int, 3) + a * x_int + b) % pp
        f_val = (y_int * y_int - rhs) % pp
        c     = f_val // p          # exact since f_val ≡ 0 (mod p)
        f_der = (2 * y_int) % p
        if f_der == 0:
            return None
        try:
            inv_der = modinv(f_der, p)
        except ValueError:
            return None
        return (y_int - c * inv_der * p) % pp

    # ── Homogeneous projective group law over Z/p^4Z ─────────────────────────
    # (X:Y:Z) represents affine (X/Z, Y/Z); O = (0:1:0).
    # All formulas are polynomial — no inversion needed.
    #
    # We work mod p^4 (not p^2) to avoid a precision issue: for the kernel
    # point p*G', the denominator V = X2*Z1 - X1*Z2 satisfies v_p(V) = 1
    # (the two inputs are "near inverses" mod p), so V^3 has v_p = 3.
    # Working mod p^2 would collapse Z3 = V^3*Z1*Z2 to 0 and falsely
    # signal the identity; mod p^4 it is correctly nonzero (v_p = 3 < 4).
    # The formal log t = -X/Y is computed mod p^2 at the end.

    pp4    = p * p * p * p
    O_PROJ = (0, 1, 0)

    def proj_double(P):
        X1, Y1, Z1 = P
        W  = (a * Z1 * Z1 + 3 * X1 * X1) % pp4
        S  = Y1 * Z1 % pp4
        B  = X1 * Y1 * S % pp4
        H  = (W * W - 8 * B) % pp4
        X3 = 2 * H * S % pp4
        Y3 = (W * (4 * B - H) - 8 * Y1 * Y1 * S * S) % pp4
        Z3 = 8 * S * S * S % pp4
        return (X3, Y3, Z3) if Z3 % pp4 != 0 else O_PROJ

    def proj_add(P, R):
        if P == O_PROJ:
            return R
        if R == O_PROJ:
            return P
        X1, Y1, Z1 = P
        X2, Y2, Z2 = R
        U  = (Y2 * Z1 - Y1 * Z2) % pp4
        V  = (X2 * Z1 - X1 * Z2) % pp4
        if V % pp4 == 0:
            return proj_double(P) if U % pp4 == 0 else O_PROJ
        V2 = V * V % pp4
        V3 = V2 * V % pp4
        A  = (U * U * Z1 * Z2 - V3 - 2 * V2 * X1 * Z2) % pp4
        X3 = V * A % pp4
        Y3 = (U * (V2 * X1 * Z2 - A) - V3 * Y1 * Z2) % pp4
        Z3 = V3 * Z1 * Z2 % pp4
        return (X3, Y3, Z3) if Z3 % pp4 != 0 else O_PROJ

    def proj_scalar(k, P):
        R = O_PROJ
        T = P
        while k:
            if k & 1:
                R = proj_add(R, T)
            T = proj_double(T)
            k >>= 1
        return R

    # ── Main computation ──────────────────────────────────────────────────────

    try:
        x_G, y_G = int(G[0]), int(G[1])
        x_Q, y_Q = int(Q[0]), int(Q[1])

        y_G_lift = lift_y(x_G, y_G)
        y_Q_lift = lift_y(x_Q, y_Q)
        if y_G_lift is None or y_Q_lift is None:
            return None

        pG = proj_scalar(p, (x_G, y_G_lift, 1))
        pQ = proj_scalar(p, (x_Q, y_Q_lift, 1))

        if pG == O_PROJ or pQ == O_PROJ:
            return None

        XG, YG, _ = pG
        XQ, YQ, _ = pQ

        # t = -X/Y (formal group parameter); divisible by p for kernel points
        try:
            t_G = (-XG * modinv(YG, pp)) % pp
            t_Q = (-XQ * modinv(YQ, pp)) % pp
        except ValueError:
            return None

        log_G = (t_G // p) % p
        log_Q = (t_Q // p) % p

        if log_G == 0:
            return None

        k = log_Q * modinv(log_G, p) % p

        if curve.scalar_mult(k, G) == Q:
            return k

    except (ValueError, ZeroDivisionError):
        pass

    return None


# ── Demo ─────────────────────────────────────────────────────────────────────

def run_demo(bit_size=18, seed=42):
    """Demonstrate Pollard's rho and Smart's attack on the ECDLP."""
    import random
    import time
    from ecc import find_curve_with_prime_order, find_anomalous_curve

    rng = random.Random(seed)

    print("=" * 60)
    print("ELLIPTIC CURVE DISCRETE LOGARITHM ATTACKS DEMONSTRATION")
    print("=" * 60)

    # --- Pollard's rho ---
    print(f"\n[*] Finding prime-order curve ({bit_size}-bit p)...")
    curve, G, q = find_curve_with_prime_order(bit_size, rng)
    p = curve.p
    print(f"    Curve : y^2 = x^3 + {curve.a}*x + {curve.b}  (mod {p})")
    print(f"    #E(F_p) = {q}  (prime, {q.bit_length()} bits)")
    print(f"    G = {G}")

    k_real = rng.randint(2, q - 1)
    Q = curve.scalar_mult(k_real, G)
    print(f"\n[*] Target: Q = k*G,  k = {k_real}  (secret)")
    print(f"    Q = {Q}")

    print(f"\n[*] Running Pollard's rho...")
    start    = time.time()
    k_rho    = pollard_rho_ecdlp(curve, G, Q, q, seed=seed)
    elapsed  = time.time() - start

    if k_rho is not None:
        print(f"    [+] Pollard's rho succeeded in {elapsed:.6f} s")
        print(f"    Recovered k = {k_rho}")
        print(f"    Correct     : {curve.scalar_mult(k_rho, G) == Q}")
    else:
        print(f"    [-] Pollard's rho failed in {elapsed:.6f} s")

    # --- Smart's attack ---
    print(f"\n[*] Finding anomalous curve ({bit_size}-bit p)...")
    try:
        acurve, AG, ap = find_anomalous_curve(bit_size, rng)
        print(f"    Curve : y^2 = x^3 + {acurve.a}*x + {acurve.b}  (mod {ap})")
        print(f"    #E(F_p) = {ap} = p  (anomalous)")
        print(f"    G = {AG}")

        k_smart_real = rng.randint(2, ap - 2)
        AQ = acurve.scalar_mult(k_smart_real, AG)
        print(f"\n[*] Target: Q = k*G,  k = {k_smart_real}  (secret)")
        print(f"    Q = {AQ}")

        print(f"\n[*] Running Smart's attack...")
        start       = time.time()
        k_smart_rec = smart_attack(acurve, AG, AQ)
        elapsed     = time.time() - start

        if k_smart_rec is not None:
            print(f"    [+] Smart's attack succeeded in {elapsed:.6f} s")
            print(f"    Recovered k = {k_smart_rec}")
            print(f"    Correct     : {acurve.scalar_mult(k_smart_rec, AG) == AQ}")
        else:
            print(f"    [-] Smart's attack failed in {elapsed:.6f} s")

    except RuntimeError as exc:
        print(f"    [!] {exc}")


# ── Scaling experiment ────────────────────────────────────────────────────────

def run_scaling_experiment_pollard(bit_sizes=None, trials=5, seed=0,
                                    csv_path=None):
    """Run Pollard's rho across different curve sizes.

    Parameters:
        bit_sizes : list of prime field sizes in bits
        trials    : independent trials per size
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time
    from ecc import find_curve_with_prime_order

    if bit_sizes is None:
        bit_sizes = [16, 18, 20]

    print("\n" + "=" * 60)
    print("POLLARD'S RHO (ECDLP) - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'p bits':>8} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 44)

    results = []

    for bits in bit_sizes:
        successes   = 0
        total_time  = 0.0

        for trial in range(trials):
            rng = random.Random(seed + trial)
            try:
                curve, G, q = find_curve_with_prime_order(bits, rng)
                k_real      = rng.randint(2, q - 1)
                Q           = curve.scalar_mult(k_real, G)

                start      = time.time()
                k_rec      = pollard_rho_ecdlp(curve, G, Q, q, seed=seed + trial)
                elapsed    = time.time() - start
                total_time += elapsed

                if k_rec is not None and curve.scalar_mult(k_rec, G) == Q:
                    successes += 1

            except Exception as ex:
                print(f"    [!] Trial {trial} failed: {ex}")

        avg_time = total_time / trials
        print(f"{bits:>8} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['pollard_rho', 'p_bits', bits, trials, successes, avg_time])

    return results


# ── Smart's attack scaling experiment ────────────────────────────────────────

def run_scaling_experiment_smart(bit_sizes=None, trials=5, seed=0,
                                  csv_path=None):
    """Run Smart's attack across different anomalous curve sizes.

    Anomalous curve search (O(p) per attempt) is the bottleneck, so p is
    kept small.  To separate attack timing from curve-finding, one anomalous
    curve is located per bit size (amortised pre-computation) and then all
    trials reuse that curve with fresh random secrets.  Only the attack
    itself is timed.

    Parameters:
        bit_sizes : list of prime field sizes in bits (default [12, 14, 16])
        trials    : independent trials per size (each with a fresh secret k)
        seed      : base random seed
        csv_path  : if given, append one CSV row per size to this file
    """
    import csv
    import random
    import time
    from ecc import find_anomalous_curve

    if bit_sizes is None:
        bit_sizes = [12, 14, 16]

    print("\n" + "=" * 60)
    print("SMART'S ATTACK - SCALING EXPERIMENT")
    print("=" * 60)
    print(f"\n{'p bits':>8} {'Trials':>7} {'Success':>9} {'Avg Time (s)':>14}")
    print("-" * 44)

    results = []

    for bits in bit_sizes:
        # Find one anomalous curve per bit size (pre-computation, not timed).
        rng_search = random.Random(seed + bits * 1000)
        try:
            acurve, AG, ap = find_anomalous_curve(bits, rng_search)
        except RuntimeError as exc:
            print(f"    [!] {bits}-bit: {exc}; skipping")
            continue

        successes  = 0
        total_time = 0.0

        for trial in range(trials):
            rng_t  = random.Random(seed + trial)
            k_real = rng_t.randint(2, ap - 2)
            Q      = acurve.scalar_mult(k_real, AG)

            start   = time.time()
            k_rec   = smart_attack(acurve, AG, Q)
            elapsed = time.time() - start
            total_time += elapsed

            if k_rec is not None and acurve.scalar_mult(k_rec, AG) == Q:
                successes += 1

        avg_time = total_time / trials
        print(f"{bits:>8} {trials:>7} {successes:>9} {avg_time:>14.6f}")
        results.append((bits, trials, successes, avg_time))

        if csv_path is not None:
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(
                    ['smart', 'p_bits', bits, trials, successes, avg_time])

    return results


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(bit_size=14, seed=42)

    run_scaling_experiment_pollard(
        bit_sizes=[16, 18, 20],
        trials=5,
        seed=0,
    )
