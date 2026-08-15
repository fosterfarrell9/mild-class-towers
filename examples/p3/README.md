# p = 3 certificates and standalone verification

This directory contains the arithmetic certificates for all 12749
imaginary quadratic fields of 3-class rank three with |D_K| < 2^30,
together with the exact strong-freeness
checks.  The certificate format is documented in `FORMAT.md`.

The certificates live under `certificates/<bucket>/K-<|D|>-p3/` and
the reference tensors under `source-tensors/<bucket>/D-<|D|>/`,
where `<bucket>` is floor(|D|/10^7) as three digits (`000`..`107`);
the sharding keeps every directory listing below GitHub's rendering
limit.  A few certificate directories carry a `hints.gp` sidecar of
proven prime factors that makes their verification deterministic;
wrong or missing hints can only slow a run down or lead to
rejection, never to a wrong acceptance.

Build the verifier and verify the certificates:

```sh
make -C ../../verifier PARI=/usr/local
python3 verify_all.py
```

The verifier is shared between the collections of all three primes and
lives in `verifier/` at the repository root.

Each certificate is accepted only after exact field, integral-basis,
Artin, ideal, norm, class-coordinate, reconstruction, and shuffle
checks; the reconstructed tensor is also compared against the stored
source tensors under `source-tensors/`.  Results are written to
`results/`; the verification records of the block
2^29 <= |D_K| < 2^30 are `results/block23-verification-starthinker.json`
and `results/block23-verification-independent.json`, its
strong-freeness verdicts `results/block23-strong-freeness.json` with
the second-engine record `results/block23-crosscheck-python.json`.

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

`make -C verifier check`, run from the repository root, first verifies
the unmodified certificates of one field per prime and then runs
`verifier/test_rejections.py`, which
generates temporary copies of that certificate and requires the
verifier to reject every one of them with the expected message: five
copies alter the arithmetic content --- the character vector, the
normalized automorphism, the multiplicity of an entry, a norm-class
vector, and the absolute field model, each separately --- and three
alter the container (an unsupported format version, a missing entry,
a wrong expected tensor).  The outcomes are recorded in
`results/rejection-tests.json`.

The verifier accepts an optional second argument, a GP vector of
primes read as factorization hints: each hint must pass a proven
primality test and is then added to PARI's prime table, which makes
the discriminant and index factorizations inside `nfinit` and
`rnfinit` deterministic for fields on which they would otherwise
depend on fortunate ECM draws.  Wrong or missing hints can only slow
a run down or lead to rejection, never to a wrong acceptance.  The
driver `verifier/verify_all.py` passes a `hints.gp` lying next to a
`certificate.gp` automatically.

## Rebuilding certificates

The builder under `builder/` regenerates certificates from scratch:

```sh
python3 builder/build_certificates.py --field 3640387
```

The relative class-field data are found by PARI as a search oracle;
every entry passes the exact audits before it is written.  A rebuild
against an existing certificate must reproduce it byte for byte
(status `UNCHANGED`); a difference aborts with a determinism failure.
