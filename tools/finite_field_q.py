#!/usr/bin/env python3
"""Finite fields F_q = F_p[t]/(m), q = p^e, with the coefficient-list
convention of the certificates (elements encoded as sum c_i p^i, digits
low to high; moduli as documented in FORMAT.md).  Used by the
rank-generic verifier."""
import itertools

MODULI = {  # low-to-high coefficient lists of monic irreducible polynomials
    (3, 1): [0, 1], (3, 2): [2, 2, 1], (3, 3): [1, 2, 0, 1], (3, 4): [2, 1, 0, 0, 1],
    (5, 1): [0, 1], (5, 2): [3, 0, 1],            # t^2 + 3 = t^2 - 2, 2 a non-square mod 5
    (5, 3): [2, 1, 0, 1],                          # t^3 + t + 2 (checked irreducible at load)
    (7, 1): [0, 1], (7, 2): [1, 0, 1],             # t^2 + 1, -1 a non-square mod 7
}


class FiniteField:
    def __init__(self, p, e, modulus=None):
        self.p, self.e = p, e
        self.q = p ** e
        self.modulus = modulus or MODULI[(p, e)]
        assert len(self.modulus) == e + 1 and self.modulus[-1] == 1
        self.zero, self.one = 0, 1
        # multiplication via polynomial arithmetic, cached in a table for small q
        self._mul = {}
        self.digits = [self._dig(x) for x in range(self.q)]
        self.enc = {tuple(d): x for x, d in enumerate(self.digits)}
        for x in range(self.q):
            for y in range(x, self.q):
                z = self._polymul(x, y)
                self._mul[(x, y)] = z; self._mul[(y, x)] = z
        self.neg = [self.enc[tuple((-a) % p for a in self.digits[x])] for x in range(self.q)]
        self.inv = [0] * self.q
        for x in range(1, self.q):
            for y in range(1, self.q):
                if self._mul[(x, y)] == 1:
                    self.inv[x] = y; break
            assert self.inv[x], f"modulus {self.modulus} not irreducible (no inverse of {x})"
        self.irreducible_checked = True

    def _dig(self, x):
        return [(x // self.p ** i) % self.p for i in range(self.e)]

    def _polymul(self, x, y):
        p, e = self.p, self.e
        a, b = self.digits[x], self.digits[y]
        prod = [0] * (2 * e - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    prod[i + j] = (prod[i + j] + ai * bj) % p
        for k in range(2 * e - 2, e - 1, -1):  # reduce t^k using modulus
            c = prod[k]
            if c:
                prod[k] = 0
                for i in range(e):
                    prod[k - e + i] = (prod[k - e + i] - c * self.modulus[i]) % p
        return self.enc[tuple(prod[:e])]

    def add(self, x, y):
        return self.enc[tuple((a + b) % self.p for a, b in zip(self.digits[x], self.digits[y]))]

    def sub(self, x, y):
        return self.add(x, self.neg[y])

    def mul(self, x, y):
        return self._mul[(x, y)]

    def frobenius(self, x):
        r = 1
        for _ in range(self.p):
            r = self.mul(r, x)
        return r

    def show(self, x):
        return self.digits[x]


def normalised(v, F):
    lead = next(x for x in v if x)
    inv = F.inv[lead]
    return tuple(F.mul(x, inv) for x in v)


