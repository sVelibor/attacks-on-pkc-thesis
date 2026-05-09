# Attacks on Public-Key Encryption Schemes

This repository contains the complete implementation accompanying the bachelor's thesis **"Attacks on Public-Key Encryption Schemes"** by **Velibor Smilevski**, submitted to the Faculty of Mathematics, Natural Sciences and Information Technologies (FAMNIT) at the University of Primorska, 2026. The codebase provides working demonstrations and scaling experiments for seven classical cryptanalytic attacks targeting RSA, Diffie–Hellman/ElGamal, and elliptic-curve discrete-logarithm problems. It is the authoritative source for the experimental results and figures presented in the thesis; the printed appendix reproduces only the core attack functions.

---

## Disclaimer

> **This code is for educational and research purposes only.**
> All attacks are demonstrated against deliberately weak or specially constructed parameters (small private exponents, smooth group orders, anomalous curves). None of the implementations constitute a threat to correctly parameterised cryptographic systems, and none should be used in or adapted for production software. The repository is published solely to make the thesis experiments reproducible and to illustrate the mathematical ideas described in the thesis.

---

## Attacks

| # | Attack | Scheme | Precondition | Complexity |
|---|--------|--------|--------------|------------|
| 1 | Wiener's attack | RSA | Private exponent d < (1/3) N^(1/4) | O(log N) — continued-fraction convergents of e/N |
| 2 | Håstad's broadcast attack | RSA | Same plaintext M encrypted under e = 3 independent public keys | O(e log N) — CRT reconstruction + integer e-th root |
| 3 | Coppersmith's short pad attack | RSA | e = 3, pad length < n/9 bits | O(2^pad_bits) in this implementation (direct search); O(poly log N) with LLL [future work] |
| 4 | Baby-step Giant-step (BSGS) | DH / ElGamal | Group order q known | O(√q) time and space — Shanks meet-in-the-middle |
| 5 | Pohlig–Hellman | DH / ElGamal | Group order q is B-smooth | O(√B · log q) — DLP in prime-power subgroups + CRT |
| 6 | Pollard's rho for ECDLP | ECC | Prime-order curve, generator order n known | O(√n) expected — pseudo-random walk with Floyd cycle detection |
| 7 | Smart's attack | ECC | Anomalous curve: #E(F_p) = p | O(log³ p) — p-adic Hensel lift + formal-group logarithm |

For the mathematical background and proofs see the thesis chapters; this README focuses on running the code.

---

## Requirements

- **Python 3.10 or later** (developed and tested on Python 3.14)
- [sympy](https://www.sympy.org/) >= 1.14 — symbolic polynomial arithmetic (Coppersmith resultant, primality)
- [matplotlib](https://matplotlib.org/) >= 3.7 — figure generation

Install dependencies:

```
pip install sympy matplotlib
```

or use the pinned versions from `requirements.txt`:

```
pip install -r requirements.txt
```

---

## Usage

### Reproduce all experiments

```
python main.py
```

This runs all seven attacks across their respective parameter ranges (5 independent trials per data point, fixed random seeds throughout). On completion it writes:

- **`results.csv`** — raw timing and success-rate data, one row per (attack, parameter value)
- **`results.txt`** — human-readable summary of the same data

All random seeds are set deterministically (`seed=42` for demos, `seed=0` for scaling experiments), so repeated runs produce identical output.

### Regenerate figures

```
python plots.py
```

Reads `results.csv` and writes eight figures as both `.pdf` (for inclusion in the LaTeX thesis) and `.png`:

```
fig_wiener_runtime.{pdf,png}
fig_hastad_runtime.{pdf,png}
fig_comparison.{pdf,png}          # RSA attack comparison
fig_bsgs_runtime.{pdf,png}
fig_pohlig_hellman_runtime.{pdf,png}
fig_dh_comparison.{pdf,png}       # DH attack comparison
fig_pollard_runtime.{pdf,png}
fig_smart_runtime.{pdf,png}
```

`python plots.py` requires a populated `results.csv`; run `python main.py` first if you do not have one.

---

## Repository structure

```
.
├── main.py            Entry point: runs all attacks, writes results.csv and results.txt
├── plots.py           Reads results.csv, generates the eight thesis figures
├── utils.py           Shared utilities: modular arithmetic, CRT, integer nth-root, RSA key generation
├── wiener.py          Wiener's attack on RSA (small private exponent)
├── hastad.py          Håstad's broadcast attack on RSA (small public exponent e = 3)
├── coppersmith.py     Coppersmith's short pad attack via polynomial resultant + direct root search
├── dh.py              Diffie–Hellman and ElGamal primitives: parameter generation, key exchange, encryption
├── dh_attacks.py      Baby-step Giant-step and Pohlig–Hellman attacks on the DLP
├── ecc.py             Elliptic-curve group arithmetic and prime-order / anomalous curve generation
├── ecc_attacks.py     Pollard's rho (ECDLP) and Smart's attack on anomalous curves
├── results.csv        Experimental data produced by main.py (committed as part of the published artifact)
├── fig_*.pdf          Eight runtime figures (PDF, for LaTeX)
└── fig_*.png          Eight runtime figures (PNG)
```

---

## Citation

If you use this code or the experimental results in your own work, please cite the thesis:

```bibtex
@mastersthesis{smilevski2026attacks,
  author  = {Velibor Smilevski},
  title   = {Attacks on Public-Key Encryption Schemes},
  school  = {University of Primorska, Faculty of Mathematics,
             Natural Sciences and Information Technologies (FAMNIT)},
  year    = {2026},
  type    = {Bachelor's thesis},
  note    = {Source code available at
             \url{https://github.com/sVelibor/attacks-on-pkc-thesis}},
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for the full text.
