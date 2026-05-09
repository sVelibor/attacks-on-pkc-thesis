# Shared utility functions used by all attack implementations

from math import gcd, isqrt


# ── Modular arithmetic ──────────────────────────────────────────────────────

def extended_gcd(a, b):
    """Extended Euclidean algorithm.
    Returns (g, x, y) such that a*x + b*y = g = gcd(a, b)."""
    if b == 0:
        return a, 1, 0
    g, x, y = extended_gcd(b, a % b)
    return g, y, x - (a // b) * y


def modinv(a, n):
    """Modular inverse of a modulo n.
    Raises ValueError if gcd(a, n) != 1."""
    g, x, _ = extended_gcd(a, n)
    if g != 1:
        raise ValueError(f"No modular inverse: gcd({a}, {n}) = {g}")
    return x % n


# ── Integer roots ───────────────────────────────────────────────────────────

def integer_nth_root(x, n):
    """Integer n-th root of x using binary search.
    Returns floor(x^(1/n)) exactly, without floating point errors."""
    if x < 0:
        raise ValueError("x must be non-negative")
    if x == 0:
        return 0
    if n == 1:
        return x

    # Initial estimate
    lo, hi = 0, min(x, 10 ** ((x.bit_length() // n) + 2))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** n <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


def is_perfect_nth_power(x, n):
    """Returns (True, root) if x is a perfect n-th power, else (False, None)."""
    root = integer_nth_root(x, n)
    if root ** n == x:
        return True, root
    return False, None


# ── Chinese Remainder Theorem ───────────────────────────────────────────────

def crt(residues, moduli):
    """Chinese Remainder Theorem.
    Given lists of residues [r1, r2, ...] and pairwise coprime
    moduli [m1, m2, ...], returns the unique x such that
    x ≡ ri (mod mi) for all i, with 0 <= x < prod(moduli)."""
    N = 1
    for m in moduli:
        N *= m

    result = 0
    for r, m in zip(residues, moduli):
        Ni = N // m
        result += r * Ni * modinv(Ni, m)

    return result % N


# ── RSA key generation ──────────────────────────────────────────────────────

def generate_prime(bits, rng):
    """Generate a random prime of the given bit length using rng."""
    from sympy import isprime, nextprime
    import random
    while True:
        candidate = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if isprime(candidate):
            return candidate


def generate_rsa_keypair(bits, rng):
    """Generate a standard RSA key pair (N, e, d, p, q).
    bits is the bit length of each prime factor."""
    from sympy import isprime
    e = 65537
    while True:
        p = generate_prime(bits, rng)
        q = generate_prime(bits, rng)
        if p == q:
            continue
        N = p * q
        phi = (p - 1) * (q - 1)
        if gcd(e, phi) == 1:
            d = modinv(e, phi)
            return N, e, d, p, q


# ── Primality check ─────────────────────────────────────────────────────────

def is_prime(n):
    """Simple primality check using sympy."""
    from sympy import isprime
    return isprime(n)