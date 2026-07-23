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
