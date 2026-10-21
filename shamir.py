"""Shamir's Secret Sharing (t-of-n threshold scheme) over GF(256) — pure Python.

Split any byte string (a wallet seed, a private key, a password) into ``n``
shares such that ANY ``t`` of them reconstruct the secret, while fewer than
``t`` reveal nothing about it at all.

This is the classic 1979 Adi Shamir scheme: for every byte of the secret a
random polynomial ``f(x)`` of degree ``t - 1`` is built with the secret byte
as ``f(0)``, and share ``i`` is the point ``(i, f(i))``. Reconstruction is
Lagrange interpolation evaluated at ``x = 0``. All arithmetic happens in the
finite field GF(256) (AES's field, irreducible polynomial ``x^8+x^4+x^3+x+1``).

Example:
    >>> shares = split(b"correct horse battery staple", n=5, threshold=3)
    >>> len(shares)
    5
    >>> # any 3 shares reconstruct the secret ...
    >>> combine(shares[:3]) == b"correct horse battery staple"
    True
    >>> combine([shares[0], shares[2], shares[4]]) == b"correct horse battery staple"
    True
    >>> # ... but 2 shares give back nothing useful (a random-looking blob)
    >>> combine(shares[:2]) == b"correct horse battery staple"
    False
"""

from __future__ import annotations

import secrets

__all__ = ["split", "combine", "encode_share", "decode_share", "FIELD_SIZE", "MAX_SHARES"]

# Irreducible polynomial for GF(2^8): x^8 + x^4 + x^3 + x + 1 (0x11B), the AES field.
_IRREDUCIBLE = 0x11B
FIELD_SIZE = 256
MAX_SHARES = 255  # share x-coordinates are 1..255 (x = 0 is the secret itself)


def _xtime(a: int) -> int:
    """Multiply a GF(256) element by x (0x02), reducing mod the irreducible poly."""
    return ((a << 1) ^ _IRREDUCIBLE) & 0xFF if a & 0x80 else (a << 1) & 0xFF


# Log/antilog tables for the generator 3 — O(1) field multiply/divide.
_EXP = [0] * 512
_LOG = [0] * 256
_acc = 1
for _i in range(511):
    _EXP[_i] = _acc
    if _i < 255:
        _LOG[_acc] = _i
    _acc = _gf_mul_raw = _xtime(_acc) ^ _acc  # multiply by 3 = x+1


def _gf_add(a: int, b: int) -> int:
    """Field addition == subtraction == XOR."""
    return a ^ b


def _gf_mul(a: int, b: int) -> int:
    """Field multiplication via log tables.

    >>> _gf_mul(0x53, 0xCA)  # AES spec example
    1
    >>> _gf_mul(0, 123)
    0
    >>> _gf_mul(200, 0)
    0
    """
    if a == 0 or b == 0:
        return 0
    return _EXP[_LOG[a] + _LOG[b]]


def _gf_div(a: int, b: int) -> int:
    """Field division.

    >>> _gf_mul(_gf_div(1, 0x53), 0x53)
    1
    >>> _gf_mul(_gf_div(0x53, 0xCA), 0xCA) == 0x53
    True
    """
    if b == 0:
        raise ZeroDivisionError("division by zero in GF(256)")
    if a == 0:
        return 0
    return _EXP[(_LOG[a] - _LOG[b]) % 255]


def _eval_poly(coeffs: list[int], x: int) -> int:
    """Evaluate polynomial (coeffs[0] + coeffs[1]*x + ...) at x via Horner."""
    result = 0
    for c in reversed(coeffs):
        result = _gf_add(_gf_mul(result, x), c)
    return result


def split(secret: bytes, n: int, threshold: int) -> list[tuple[int, bytes]]:
    """Split ``secret`` into ``n`` shares; any ``threshold`` of them rebuild it.

    Returns a list of ``(x, y_bytes)`` tuples with ``x`` in 1..n.

    >>> shares = split(b"abc", n=3, threshold=2)
    >>> [x for x, _ in shares]
    [1, 2, 3]
    >>> all(len(y) == 3 for _, y in shares)
    True
    >>> split(b"", n=3, threshold=2)
    Traceback (most recent call last):
    ...
    ValueError: secret must be non-empty
    >>> split(b"abc", n=5, threshold=6)
    Traceback (most recent call last):
    ...
    ValueError: threshold must satisfy 1 <= threshold <= n <= 255
    >>> split(b"abc", n=300, threshold=2)
    Traceback (most recent call last):
    ...
    ValueError: threshold must satisfy 1 <= threshold <= n <= 255
    """
    if not secret:
        raise ValueError("secret must be non-empty")
    if not (1 <= threshold <= n <= MAX_SHARES):
        raise ValueError("threshold must satisfy 1 <= threshold <= n <= 255")

    # One random polynomial of degree threshold-1 per secret byte, fixed across
    # all shares: f_j(x) = secret_j + a_1j*x + ... (f_j(0) == secret byte j).
    polys = [
        [byte] + [secrets.randbelow(256) for _ in range(threshold - 1)]
        for byte in secret
    ]
    shares: list[tuple[int, bytes]] = []
    for x in range(1, n + 1):
        y = bytes(_eval_poly(coeffs, x) for coeffs in polys)
        shares.append((x, y))
    return shares


def combine(shares: list[tuple[int, bytes]]) -> bytes:
    """Reconstruct the secret from shares via Lagrange interpolation at x = 0.

    Needs at least ``threshold`` DISTINCT shares; with fewer shares the output
    is unrelated to the secret (this is the security guarantee, not a bug).

    >>> s = b"the eagle has landed"
    >>> shares = split(s, n=5, threshold=3)
    >>> combine(shares[1:4]) == s
    True
    >>> combine([shares[4], shares[0], shares[3]]) == s  # order doesn't matter
    True
    >>> combine([shares[0], shares[0], shares[1]])  # duplicate x rejected
    Traceback (most recent call last):
    ...
    ValueError: share x-coordinates must be distinct
    >>> combine([])
    Traceback (most recent call last):
    ...
    ValueError: need at least one share
    """
    if not shares:
        raise ValueError("need at least one share")
    xs = [x for x, _ in shares]
    if len(set(xs)) != len(xs):
        raise ValueError("share x-coordinates must be distinct")
    if any(len(y) != len(shares[0][1]) for _, y in shares):
        raise ValueError("all shares must have the same length")

    secret = bytearray()
    for byte_idx in range(len(shares[0][1])):
        acc = 0
        for j, (xj, yj) in enumerate(shares):
            # Lagrange basis l_j(0) = prod_{m != j} x_m / (x_m - x_j); note
            # x_m - x_j == x_m + x_j in characteristic 2.
            num, den = 1, 1
            for m, (xm, _) in enumerate(shares):
                if m == j:
                    continue
                num = _gf_mul(num, xm)
                den = _gf_mul(den, _gf_add(xm, xj))
            acc = _gf_add(acc, _gf_mul(yj[byte_idx], _gf_div(num, den)))
        secret.append(acc)
    return bytes(secret)


def encode_share(share: tuple[int, bytes]) -> str:
    """Encode a share as a compact portable hex string ``<t>-<xx>-<data>``.

    >>> shares = split(b"hi", n=3, threshold=2)
    >>> enc = encode_share(shares[0])
    >>> enc.split("-")[0]
    '01'
    >>> decode_share(enc) == shares[0]
    True
    """
    x, y = share
    return f"{x:02x}-{y.hex()}"


def decode_share(text: str) -> tuple[int, bytes]:
    """Inverse of :func:`encode_share`.

    >>> decode_share("zz-00")
    Traceback (most recent call last):
    ...
    ValueError: malformed share string
    """
    try:
        x_hex, y_hex = text.split("-")
        x = int(x_hex, 16)
        y = bytes.fromhex(y_hex)
    except ValueError:
        raise ValueError("malformed share string") from None
    if not (1 <= x <= MAX_SHARES) or not y:
        raise ValueError("malformed share string")
    return x, y


def _demo() -> None:
    seed = bytes.fromhex(
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    )
    print(f"secret seed (32 bytes): {seed.hex()}\n")

    n, t = 5, 3
    shares = split(seed, n=n, threshold=t)
    print(f"split into {n} shares, threshold {t}:")
    for s in shares:
        print("  ", encode_share(s))

    print(f"\nreconstruct from shares 1,3,5: {combine([shares[0], shares[2], shares[4]]).hex()}")
    wrong = combine(shares[:2])
    print(f"reconstruct from only 2 shares: {wrong.hex()}")
    print(f"  -> matches secret? {wrong == seed} (fewer than t shares reveal nothing)")


def _self_test() -> None:
    import doctest

    results = doctest.testmod(verbose=False)
    print(f"doctests: {results.attempted - results.failed}/{results.attempted} passed")
    if results.failed:
        raise SystemExit(1)


if __name__ == "__main__":
    import sys

    if "--test" in sys.argv:
        _self_test()
    else:
        _demo()
