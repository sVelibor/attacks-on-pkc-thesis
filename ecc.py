# ecc.py
# Elliptic Curve Cryptography over prime fields
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999
#           and Hankerson, Menezes, Vanstone, "Guide to Elliptic Curve Cryptography"

from math import gcd
from utils import modinv, generate_prime


# ── Elliptic curve group over F_p ────────────────────────────────────────────

class EllipticCurve:
    """Short Weierstrass elliptic curve: y^2 ≡ x^3 + ax + b (mod p).

    The point at infinity is represented as None (the group identity).
    All coordinates are integers modulo p.
    """

    def __init__(self, p, a, b):
        """Initialise curve E: y^2 = x^3 + ax + b over F_p.

        Parameters:
            p : prime field characteristic
            a : curve coefficient a (mod p)
            b : curve coefficient b (mod p)
        """
        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
            raise ValueError("Discriminant is zero: curve is singular")
        self.p = p
        self.a = a % p
        self.b = b % p

    def is_on_curve(self, P):
        """Return True if P lies on the curve (or is the point at infinity)."""
        if P is None:
            return True
        x, y = P
        return (pow(y, 2, self.p) - pow(x, 3, self.p) - self.a * x - self.b) % self.p == 0

    def add(self, P, Q):
        """Add two points P and Q on the curve."""
        if P is None:
            return Q
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q
        p = self.p

        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None  # P + (-P) = point at infinity
            return self.double(P)

        lam = (y2 - y1) * modinv(x2 - x1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def double(self, P):
        """Double a point P on the curve."""
        if P is None:
            return None

        x1, y1 = P
        p = self.p

        if y1 == 0:
            return None  # tangent is vertical

        lam = (3 * x1 * x1 + self.a) * modinv(2 * y1, p) % p
        x3 = (lam * lam - 2 * x1) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def scalar_mult(self, k, P):
        """Compute kP using the double-and-add algorithm.

        Parameters:
            k : scalar (non-negative integer)
            P : curve point or None (point at infinity)

        Returns:
            kP as a curve point or None.
        """
        if k < 0:
            # k*P = (-k)*(-P); negate by flipping y
            k = -k
            if P is not None:
                P = (P[0], (-P[1]) % self.p)

        result = None  # additive identity
        addend = P

        while k:
            if k & 1:
                result = self.add(result, addend)
            addend = self.double(addend)
            k >>= 1

        return result

    def neg(self, P):
        """Return the additive inverse of P."""
        if P is None:
            return None
        return (P[0], (-P[1]) % self.p)

    def point_count_naive(self, limit=None):
        """Count curve points by exhaustive enumeration over F_p.

        Complexity O(p) — only suitable for very small p.
        Returns #E(F_p) including the point at infinity.

        Parameters:
            limit : if given, stop early once the count exceeds limit
        """
        count = 1  # point at infinity
        for x in range(self.p):
            rhs = (pow(x, 3, self.p) + self.a * x + self.b) % self.p
            # Count y such that y^2 ≡ rhs (mod p)
            if rhs == 0:
                count += 1
            else:
                # Euler criterion: rhs^((p-1)/2) ≡ 1 => two solutions
                if pow(rhs, (self.p - 1) // 2, self.p) == 1:
                    count += 2
            if limit is not None and count > limit:
                return count
        return count

    def __repr__(self):
        return f"EllipticCurve(p={self.p}, a={self.a}, b={self.b})"


# ── Curve constructors ────────────────────────────────────────────────────────

def find_curve_with_prime_order(bit_size, rng, max_tries=2000):
    """Find an elliptic curve E over a random prime field F_p whose group
    order #E(F_p) is prime (a prime-order curve).

    Uses small primes p (up to ~30 bits) and exhaustive point counting so
    that this remains fast for demonstration purposes.

    Parameters:
        bit_size  : bit length of the prime p (recommended 14–24)
        rng       : seeded random.Random instance
        max_tries : abort after this many candidate curves

    Returns:
        (curve, G, q) where G is a generator of order q = #E(F_p).
    """
    from sympy import isprime

    for _ in range(max_tries):
        p = generate_prime(bit_size, rng)
        a = rng.randint(0, p - 1)
        b = rng.randint(0, p - 1)

        # Reject singular curves
        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
            continue

        try:
            curve = EllipticCurve(p, a, b)
        except ValueError:
            continue

        q = curve.point_count_naive()

        if not isprime(q):
            continue

        # Find a generator: any affine point has order q (prime-order group)
        G = _find_affine_point(curve, rng)
        if G is None:
            continue

        return curve, G, q

    raise RuntimeError("Could not find a prime-order curve in the given number of tries")


def find_anomalous_curve(bit_size, rng, max_tries=5000):
    """Find an anomalous elliptic curve: #E(F_p) = p.

    Smart's attack applies to exactly this class of curve.

    Parameters:
        bit_size  : bit length of the prime p (recommended 14–22)
        rng       : seeded random.Random instance
        max_tries : abort after this many candidate primes/curves

    Returns:
        (curve, G, p) where G is an affine point of order p = #E(F_p).
    """
    from sympy import isprime

    for _ in range(max_tries):
        p = generate_prime(bit_size, rng)

        # Try several (a, b) pairs for this p
        for __ in range(50):
            a = rng.randint(1, p - 1)
            b = rng.randint(1, p - 1)

            if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
                continue

            try:
                curve = EllipticCurve(p, a, b)
            except ValueError:
                continue

            cnt = curve.point_count_naive(limit=p + int(2 * p**0.5) + 2)

            if cnt == p:
                G = _find_affine_point(curve, rng)
                if G is not None:
                    return curve, G, p

    raise RuntimeError("Could not find an anomalous curve in the given number of tries")


def _find_affine_point(curve, rng):
    """Return a random affine point on curve by trial over random x values."""
    p = curve.p
    for _ in range(200):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + curve.a * x + curve.b) % p
        if rhs == 0:
            return (x, 0)
        if pow(rhs, (p - 1) // 2, p) == 1:
            # Compute square root via Tonelli-Shanks (simplified: p | 4k+3 case)
            y = _sqrt_mod(rhs, p)
            if y is not None:
                return (x, y)
    return None


def _sqrt_mod(n, p):
    """Compute a square root of n modulo prime p using Tonelli-Shanks."""
    if n % p == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None  # not a QR

    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # General Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1

    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    M = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    R = pow(n, (q + 1) // 2, p)

    while True:
        if t == 1:
            return R
        i, tmp = 1, (t * t) % p
        while tmp != 1:
            tmp = (tmp * tmp) % p
            i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p


# ── Demo ─────────────────────────────────────────────────────────────────────

def run_demo(bit_size=18, seed=42):
    """Demonstrate elliptic curve group operations over a small prime field."""
    import random

    rng = random.Random(seed)

    print("=" * 60)
    print("ELLIPTIC CURVE GROUP DEMONSTRATION")
    print("=" * 60)

    print(f"\n[*] Searching for prime-order curve ({bit_size}-bit p)...")
    curve, G, q = find_curve_with_prime_order(bit_size, rng)
    p = curve.p

    print(f"    Curve : y^2 = x^3 + {curve.a}*x + {curve.b}  (mod {p})")
    print(f"    #E(F_p) = {q}  (prime order)")
    print(f"    Generator G = {G}")

    k = rng.randint(2, q - 1)
    kG = curve.scalar_mult(k, G)
    print(f"\n[*] Scalar multiplication: k = {k}")
    print(f"    kG = {kG}")
    print(f"    G on curve: {curve.is_on_curve(G)}")
    print(f"    kG on curve: {curve.is_on_curve(kG)}")

    # Verify group order: qG = O
    qG = curve.scalar_mult(q, G)
    print(f"\n[*] Order check: q*G = {qG}  (should be None/infinity)")

    # Demonstrate anomalous curve search
    print(f"\n[*] Searching for anomalous curve (anomalous => #E(F_p) = p)...")
    try:
        acurve, AG, ap = find_anomalous_curve(bit_size, rng)
        print(f"    Anomalous curve: y^2 = x^3 + {acurve.a}*x + {acurve.b}  (mod {ap})")
        print(f"    #E(F_p) = {ap} = p  (anomalous confirmed)")
    except RuntimeError as exc:
        print(f"    [!] {exc}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(bit_size=14, seed=42)
