# p = 3 certificates and standalone verification

This directory contains the arithmetic certificates for the eleven
fields of the eleven-fields theorem of the p = 3 companion paper, plus
the pilot field `D = -3640387`, together with a standalone C verifier
and the exact strong-freeness checks.  The certificate format is
documented in `FORMAT.md`.

Build the verifier and verify all twelve certificates:

```sh
make -C verifier PARI=/usr/local
python3 verifier/verify_all.py
```

Each certificate is accepted only after exact field, integral-basis,
Artin, ideal, norm, class-coordinate, reconstruction, and shuffle
checks; the reconstructed tensor is also compared against the stored
source tensors under `source-tensors/`.  Results are written to
`results/`.

Re-check the strong-freeness certificates on the verifier output:

```sh
python3 strong-freeness/verify_strong_freeness.py
```

This confirms the exact Anick witnesses (rational and over F9 with
descent) and the terminating Groebner/automaton certificate, closing
the chain from certificate to STRONGLY_FREE for the eleven theorem
fields; the pilot field remains UNDECIDED.  The canonical summaries
are `results/verification.json` and `results/final-results.json`.

The negative tests in `verifier/test_rejections.py` check that
malformed certificates, missing entries, and wrong expected tensors
are rejected.
