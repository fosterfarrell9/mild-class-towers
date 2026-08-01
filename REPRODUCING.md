# Reproducing the computations

Build against a local PARI 2.17.4 installation:

```sh
make clean
make PARI="$HOME/.local"
```

The PARI build must expose three routines that are `static` in stock
PARI 2.17.4 (`rnfcycaut`, `allauts`, `cyclicrelfrob`); the required
one-file patch, the rebuild steps, and a regression test are documented
in `doc/pari-2.17.4-patch.md`.  The whole build is scripted:

```sh
tools/build-patched-pari.sh          # installs into $HOME/.local
```

Run the fixed rank-three arithmetic and tensor computation for
\(K=\mathbf Q(\sqrt{-2800905})\) and \(p=5\):

```sh
MASSEY_RANK3_TEST=1 ./build/massey 5 "s^2+2800905"
```

Run the finite cubic relation-space and strong-freeness computation:

```sh
MASSEY_CERTIFICATE_TEST=1 ./build/massey
```

Run the direct arithmetic verification of the secondary norm computations:

```sh
MASSEY_ARITHMETIC_AUDIT=1 ./build/massey 5 "s^2+2800905"
```

The ordinary rank-two example can be used as a quick regression check:

```sh
./build/massey 5 "s^2+22717"
```

The rank-three and arithmetic-verification modes perform number-field
arithmetic and may take substantially longer than the finite relation test.

## Arithmetic certificates

A certificate is a per-field text file that pins down the expensive
arithmetic behind the six secondary-norm matrices: for every pair
(character, torsion basis element) it stores the relative and absolute
polynomials of the class field, the normalized automorphism, the pair
(a', J) representing the torsion class, the auxiliary ideal I', the
compactly represented element t_AC, a residue prime for the sign check,
and the expected norm class.  The standalone verifier recomputes and
certifies the base field, identifies the class field and the normalized
automorphism from the stored data, checks AC1 and AC2 by exact ideal
arithmetic, and recomputes each norm class -- all without repeating the
search that found the data.

Certificates do not cover the strong-freeness searches: a failed search
has no succinct witness.  They do not need to, because everything after
the certified matrices is fast finite algebra: the cubic matrix T
follows by the reconstruction formula, the exhaustive GL_3(F_5) search
takes under a minute, and the Groebner tools above re-verify the
Hilbert evidence in minutes.

The verifier is shared between all certified fields and lives in
`certificate/`, next to the per-field data directories.  Build it once
and verify a certificate (the optional second argument cross-checks the
certified matrices against a committed result record, reporting
RESULT_RECORD_MATCH=PASS):

```sh
cd certificate
make PARI="$HOME/.local"
./verify_certificate K-2800905-p5/certificate.gp
./verify_certificate K-51213139-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-51213139/result.gp
```

Export a certificate for a further field by running the audited
pipeline with the export path set (expensive; it repeats the search):

```sh
MASSEY_CERTIFICATE_EXPORT="$PWD/certificate/K-<n>-p5/certificate.gp" \
./build/massey --example-result /tmp/result.gp \
  --strong-search-limit exhaustive 5 '<polynomial>'
```

Pipeline-exported certificates contain the 18 entries for the six
characters.  The principal example's certificate additionally contains
the nine genuinely computed doubled-character entries (27 in total),
produced by the fixed audit mode described in
`certificate/K-2800905-p5/README.md`; the verifier accepts both forms
and performs the doubled-character checks only when the entries are
present.

## Transverse rank-one certificates

The transverse rank-one mildness criterion is certified for all nine
computed fields by

```sh
python3 tools/transverse_rank_one.py
```

which rebuilds the quadratic secondary-norm family of each field by
polarization from the six verified matrices (committed `result.gp`
records; for the principal example
`certificate/K-2800905-p5/secondary-norms.gp`, exported from the
arithmetic certificate by `gp -qf tools/export_secondary_norms.gp`),
locates the closed points of the norm-degeneracy scheme exactly (Singular
radical per affine chart, eliminant roots over explicit fields
F_5[t]/(f)), and certifies `rank(D_x)=1` and `det(B_x)!=0` by direct
linear algebra over the residue field, with a Jacobian tangent-space
cross-check at every point.  Results are committed under
`examples/p5/transverse-rank-one/`.

## Pure relation-algebra searches for a stored matrix

The relation-algebra fixture (lexicographic and weighted-degree Anick
searches plus the Efrat section-8.1 ordered-monoid families) runs on the
audited matrix of the principal example by default:

```sh
MASSEY_RELATION_TEST=1 ./build/massey
```

To run the same searches on any other cubic relation matrix, pass the
3 x 27 matrix over F_5 in GP syntax:

```sh
MASSEY_RELATION_TEST=1 \
MASSEY_RELATION_MATRIX="[...;...;...]" ./build/massey
```

The matrix of a computed example is the `cubic_relation_matrix` entry of
its `result.gp` record.

## Groebner test for strong freeness

Independently of the order searches, a degree-truncated noncommutative
Groebner (diamond-lemma) completion can decide strong freeness of the
cubic relation space when it terminates, and rigorously verifies the
Hilbert series of the relation algebra degree by degree when it does not:

```sh
python3 tools/strong_freeness_gb.py \
  --result examples/p5/D-18397407/result.gp --order xyz --max-degree 12
python3 tests/test_strong_freeness_gb.py
```

A verdict STRONGLY_FREE proves mildness of the corresponding tower group;
NOT_STRONGLY_FREE proves that it is not mild with respect to the
Zassenhaus filtration; INCONCLUSIVE_SERIES_MATCHES records that the
Hilbert series agrees with the strongly free prediction through the
processed degree.

The same computation can be driven through Singular's Letterplace
subsystem as an independent second engine (requires the `Singular`
executable; higher degree bounds become practical):

```sh
python3 tools/strong_freeness_singular.py \
  --result examples/p5/D-18397407/result.gp --order xyz --degree-bound 13
python3 tests/test_strong_freeness_singular.py
```

Only the truncated Groebner leading words are taken from Singular; the
normal-word counting and the series comparison reuse the routines of
`tools/strong_freeness_gb.py`.  `--emit-script` prints the generated
Singular input for archival instead of running it.
