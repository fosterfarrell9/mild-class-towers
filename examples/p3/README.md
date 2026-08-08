# p = 3 certificates and standalone verification

This directory contains the arithmetic certificates for all 2497
imaginary quadratic fields of 3-class rank three with |D_K| < 2^28,
together with a standalone C verifier and the exact strong-freeness
checks.  The certificate format is documented in `FORMAT.md`.

Build the verifier and verify the certificates:

```sh
make -C verifier PARI=/usr/local
python3 verifier/verify_all.py
```

Each certificate is accepted only after exact field, integral-basis,
Artin, ideal, norm, class-coordinate, reconstruction, and shuffle
checks; the reconstructed tensor is also compared against the stored
source tensors under `source-tensors/`.  Results are written to
`results/`; the complete block verification is
`results/verification-004.json`.

## Strong freeness

Three drivers read the verified tensors and decide strong freeness of
the cubic relation space.  None of them recomputes any arithmetic.

```sh
python3 strong-freeness/verify_strong_freeness.py     # twelve fields
python3 strong-freeness/block_witnesses.py            # five block witnesses
python3 strong-freeness/block_sweep.py --degree-bound 13 --only <D> ...
```

The first confirms the Anick witnesses (rational and over F9 with
descent) and the terminating Groebner/automaton certificate for the
twelve fields of `results/verification-002.json`; the pilot field
`D = -3640387` remains UNDECIDED.  The second recomputes the five
witnesses that lie beyond the cone criterion.

`block_sweep.py` completes the two-sided Groebner basis of the cubic
relation ideal to a degree bound and compares the normal word counts
with the strongly free series 1/(1-3z+3z^3).  A field is
STRONGLY_FREE when the completion terminates by the diamond-lemma
bound and the series agrees, NOT_STRONGLY_FREE when a coefficient
deviates, and UNDECIDED otherwise.  Its verdicts for the 69 decided
fields of the block are `results/strong-freeness-block-001.json`; to
re-derive them,

```sh
python3 strong-freeness/block_sweep.py --degree-bound 13 \
    --recheck ../results/strong-freeness-block-001.json
```

which takes seconds, since a decided field either terminates with a
small basis or deviates in low degree.  Searching the whole block is
what is expensive.  With `--engine singular` the same driver computes
the bases with Singular's Letterplace subsystem instead; every
recorded verdict was produced by both engines.

The negative tests in `verifier/test_rejections.py` check that
malformed certificates, missing entries, and wrong expected tensors
are rejected.

## Rebuilding certificates

The builder under `builder/` regenerates certificates from scratch:

```sh
python3 builder/build_certificates.py --field 3640387
```

The relative class-field data are found by PARI as a search oracle;
every entry passes the exact audits before it is written.  A rebuild
against an existing certificate must reproduce it byte for byte
(status `UNCHANGED`); a difference aborts with a determinism failure.
