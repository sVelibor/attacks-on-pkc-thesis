# dh.py
# Diffie-Hellman key exchange and ElGamal encryption
# Based on: D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", 1999
#           and Menezes, van Oorschot, Vanstone, "Handbook of Applied Cryptography"

from math import gcd
from utils import modinv, generate_prime


# ── DH parameter generation ──────────────────────────────────────────────────

def generate_dh_params(prime_bits, rng):
    """Generate Diffie-Hellman parameters (p, q, g) using a safe prime approach.

    Generates a Sophie Germain prime pair: q prime, p = 2q + 1 prime.
    The generator g has order q in Z_p^*.

    Parameters:
        prime_bits : bit length of the subgroup order q
        rng        : seeded random.Random instance

    Returns:
        (p, q, g) where p is the modulus, q is the prime subgroup order,
        and g is a generator of the order-q subgroup.
    """
    from sympy import isprime

    for _ in range(10000):
        q = generate_prime(prime_bits, rng)
        p = 2 * q + 1
        if not isprime(p):
            continue

        # Find generator of the order-q subgroup
        # For safe prime p = 2q+1, elements of order q satisfy g^q ≡ 1 (mod p)
        # and g != 1. Any h with h^2 != 1 and h^q != 1 gives g = h^2 mod p
        # or g = h directly if h^q = 1.
        for _ in range(200):
            h = rng.randint(2, p - 2)
            g = pow(h, 2, p)
            if g != 1 and pow(g, q, p) == 1:
                return p, q, g

    raise RuntimeError("Could not generate DH parameters")


def generate_smooth_dh_params(prime_bits, rng):
    """Generate DH parameters where q - 1 is smooth (has only small prime factors).

    This makes the group order vulnerable to Pohlig-Hellman.
    Constructs q as a product of small primes bounded by 2^smooth_bits,
    then finds a prime p with q | p - 1 and a generator g of order q.

    Parameters:
        prime_bits : approximate bit length of the resulting modulus p
        rng        : seeded random.Random instance

    Returns:
        (p, q, g, factorization) where factorization is the prime factorisation
        of q as a list of (prime, exponent) pairs.
    """
    from sympy import isprime, factorint

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    for _ in range(5000):
        # Build q as a product of small primes so that q has ~prime_bits bits
        q = 1
        factors = {}
        candidates = list(small_primes)
        rng.shuffle(candidates)

        for sp in candidates:
            exp = rng.randint(1, 3)
            for _ in range(exp):
                if q.bit_length() >= prime_bits:
                    break
                if not isprime(sp):
                    continue
                q *= sp
                factors[sp] = factors.get(sp, 0) + 1
            if q.bit_length() >= prime_bits:
                break

        if q < 2:
            continue

        # Find prime p = k*q + 1 for some k
        for k in range(2, 3000, 2):
            p = k * q + 1
            if p.bit_length() > prime_bits + 10:
                break
            if not isprime(p):
                continue

            # Find a generator of order exactly q in Z_p^*
            # Check g^(q/l) != 1 for every prime l | q to rule out sub-orders.
            prime_factors = list(factors.keys())
            for _ in range(200):
                h = rng.randint(2, p - 2)
                g = pow(h, (p - 1) // q, p)
                if g == 1 or pow(g, q, p) != 1:
                    continue
                if all(pow(g, q // l, p) != 1 for l in prime_factors):
                    factorization = sorted(factors.items())
                    return p, q, g, factorization

    raise RuntimeError("Could not generate smooth-order DH parameters")


# ── DH key exchange ───────────────────────────────────────────────────────────

def dh_keygen(p, q, g, rng):
    """Generate a DH key pair.

    Parameters:
        p, q, g : DH parameters
        rng     : seeded random.Random instance

    Returns:
        (private_key, public_key) where public_key = g^private_key mod p.
    """
    private = rng.randint(2, q - 1)
    public = pow(g, private, p)
    return private, public


def dh_shared_secret(p, their_public, my_private):
    """Compute the Diffie-Hellman shared secret.

    Returns their_public^my_private mod p.
    """
    return pow(their_public, my_private, p)


# ── ElGamal encryption ────────────────────────────────────────────────────────

def elgamal_keygen(p, q, g, rng):
    """Generate an ElGamal key pair.

    Returns:
        (x, h) where x is the private key and h = g^x mod p is the public key.
    """
    x = rng.randint(2, q - 1)
    h = pow(g, x, p)
    return x, h


def elgamal_encrypt(p, g, h, M, rng, q=None):
    """Encrypt plaintext M under ElGamal public key (p, g, h).

    Ciphertext: (c1, c2) = (g^k mod p, M * h^k mod p).

    Parameters:
        p, g, h : public key parameters
        M       : plaintext integer, 1 <= M <= p-2
        rng     : seeded random.Random instance
        q       : subgroup order (used to bound k; defaults to p-2)

    Returns:
        (c1, c2) ciphertext pair.
    """
    bound = q if q is not None else p - 2
    k = rng.randint(2, bound - 1)
    c1 = pow(g, k, p)
    c2 = (M * pow(h, k, p)) % p
    return c1, c2


def elgamal_decrypt(p, x, c1, c2):
    """Decrypt ElGamal ciphertext (c1, c2) with private key x.

    M = c2 * (c1^x)^(-1) mod p.
    """
    s = pow(c1, x, p)
    return (c2 * modinv(s, p)) % p


# ── Demo ─────────────────────────────────────────────────────────────────────

def run_demo(prime_bits=30, seed=42):
    """Demonstrate DH key exchange and ElGamal encryption."""
    import random

    rng = random.Random(seed)

    print("=" * 60)
    print("DIFFIE-HELLMAN KEY EXCHANGE AND ELGAMAL DEMONSTRATION")
    print("=" * 60)

    # --- DH key exchange ---
    print(f"\n[*] Generating DH parameters ({prime_bits}-bit subgroup order)...")
    p, q, g = generate_dh_params(prime_bits, rng)
    print(f"    p = {p}  ({p.bit_length()} bits)")
    print(f"    q = {q}  (subgroup order, prime)")
    print(f"    g = {g}  (generator, order q)")

    print(f"\n[*] DH key exchange (Alice and Bob)...")
    a_priv, a_pub = dh_keygen(p, q, g, rng)
    b_priv, b_pub = dh_keygen(p, q, g, rng)

    alice_secret = dh_shared_secret(p, b_pub, a_priv)
    bob_secret   = dh_shared_secret(p, a_pub, b_priv)

    print(f"    Alice private key  : {a_priv}")
    print(f"    Alice public key   : {a_pub}")
    print(f"    Bob   private key  : {b_priv}")
    print(f"    Bob   public key   : {b_pub}")
    print(f"    Alice shared secret: {alice_secret}")
    print(f"    Bob   shared secret: {bob_secret}")
    print(f"    Secrets match      : {alice_secret == bob_secret}")

    # --- ElGamal ---
    print(f"\n[*] ElGamal encryption...")
    x, h = elgamal_keygen(p, q, g, rng)
    M = rng.randint(2, p - 2)
    c1, c2 = elgamal_encrypt(p, g, h, M, rng, q)
    M_dec = elgamal_decrypt(p, x, c1, c2)

    print(f"    Private key x   : {x}")
    print(f"    Public key h    : {h}")
    print(f"    Plaintext M     : {M}")
    print(f"    Ciphertext      : (c1={c1}, c2={c2})")
    print(f"    Decrypted M     : {M_dec}")
    print(f"    Correct         : {M_dec == M}")

    # --- Smooth-order group (vulnerable to Pohlig-Hellman) ---
    print(f"\n[*] Generating smooth-order DH parameters for Pohlig-Hellman demo...")
    p_s, q_s, g_s, fact_s = generate_smooth_dh_params(prime_bits, rng)
    print(f"    p = {p_s}")
    print(f"    q = {q_s}  (smooth order)")
    print(f"    Factorisation q = {' * '.join(f'{pr}^{ex}' for pr, ex in fact_s)}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_demo(prime_bits=30, seed=42)
