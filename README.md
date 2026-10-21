# shamir-secret-share

Shamir's Secret Sharing (t-of-n threshold scheme) over GF(256) — in pure Python, zero dependencies.

Split any byte string — a wallet seed, a private key, a password — into `n` shares such that **any** `t` of them reconstruct the secret, while fewer than `t` reveal nothing about it at all. This is the classic 1979 Adi Shamir scheme, the same math behind seed-phrase backup splitting (e.g. SLIP-39 style social recovery).

## How it works

For every byte of the secret, a random polynomial `f(x)` of degree `t − 1` is built with the secret byte as `f(0)`; share `i` is the point `(i, f(i))`. Reconstruction is Lagrange interpolation evaluated at `x = 0`. All arithmetic happens in the finite field GF(256) (AES's field, irreducible polynomial `x⁸+x⁴+x³+x+1`), with O(1) multiply/divide via log/antilog tables.

## Usage

```python
from shamir import split, combine, encode_share, decode_share

# split a 32-byte seed into 5 shares; any 3 recover it
shares = split(seed_bytes, n=5, threshold=3)

# hand out portable hex-encoded shares
for s in shares:
    print(encode_share(s))   # e.g. "01-447b8841..."

# later: any 3 shares rebuild the secret (order doesn't matter)
secret = combine([shares[0], shares[2], shares[4]])
```

```bash
python3 shamir.py            # demo: split, reconstruct, show t-1 shares give garbage
python3 shamir.py --test     # run the 27 doctests
```

## Notes

- Educational implementation, not audited crypto. Randomness comes from `secrets` (CSPRNG).
- Share x-coordinates are `1..255`, so `n ≤ 255`.
- Corrupting even one byte of a share silently yields a wrong secret — in production, pair shares with a checksum/MAC if you need tamper *detection*.
