# Reproducing the computations

Build against a local PARI 2.17.4 installation:

```sh
make clean
make PARI="$HOME/.local"
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
  --result examples/p5/K-18397407/result.gp --order xyz --max-degree 12
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
  --result examples/p5/K-18397407/result.gp --order xyz --degree-bound 13
python3 tests/test_strong_freeness_singular.py
```

Only the truncated Groebner leading words are taken from Singular; the
normal-word counting and the series comparison reuse the routines of
`tools/strong_freeness_gb.py`.  `--emit-script` prints the generated
Singular input for archival instead of running it.
