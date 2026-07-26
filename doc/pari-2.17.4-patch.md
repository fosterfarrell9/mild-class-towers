# PARI 2.17.4 patch required by Massey-pari

Every computation in this repository that touches Artin symbols relies
on this patch: the original Artin-symbol implementation, the audited
rank-three pipeline, and the standalone certificate verifier in
`certificate/`.

`Massey-pari` uses three internal PARI functions from
`src/modules/algebras.c`:

- `cyclicrelfrob`
- `allauts`
- `rnfcycaut`

In stock PARI 2.17.4 these functions are declared `static`, so they are
not exported by `libpari` and cannot be called from `Massey-pari`.

The original `Massey-pari` source assumes that these functions are made
available externally.

## PARI source patch

In

```text
src/modules/algebras.c
```

change the following three definitions.

### 1. `cyclicrelfrob`

Change

```c
static long
cyclicrelfrob(GEN rnf, GEN auts, GEN pr)
```

to

```c
long
cyclicrelfrob(GEN rnf, GEN auts, GEN pr)
```

### 2. `allauts`

Change

```c
static GEN
allauts(GEN rnf, GEN aut)
```

to

```c
GEN
allauts(GEN rnf, GEN aut)
```

### 3. `rnfcycaut`

Change

```c
static GEN
rnfcycaut(GEN rnf)
```

to

```c
GEN
rnfcycaut(GEN rnf)
```

No other changes to the PARI source are required.

## Rebuild and install PARI

PARI was configured with

```bash
./Configure --prefix=$HOME/.local
```

After applying the patch, rebuild and reinstall:

```bash
make -j$(nproc)
make install
```

Check that the three symbols are now exported:

```bash
nm -D --defined-only ~/.local/lib/libpari.so \
  | grep -E ' (rnfcycaut|allauts|cyclicrelfrob)$'
```

All three functions should appear in the output.

## Massey-pari declarations

`Massey-pari` contains the local header

```text
headers/pari_internal.h
```

with the corresponding declarations:

```c
#ifndef MASSEY_PARI_INTERNAL_H
#define MASSEY_PARI_INTERNAL_H

#include <pari/pari.h>

GEN rnfcycaut(GEN rnf);
GEN allauts(GEN rnf, GEN aut);
long cyclicrelfrob(GEN rnf, GEN auts, GEN pr);

#endif
```

`src/artin_symbol.c` includes this header.

## Local PARI installation

The patched PARI 2.17.4 installation used for development is:

```text
$HOME/.local/include/pari/pari.h
$HOME/.local/lib/libpari.so
```

`Massey-pari` is built with

```bash
make PARI=$HOME/.local
```

and the Makefile uses an rpath pointing to

```text
$HOME/.local/lib
```

so that the executable does not accidentally load the older
system-wide Ubuntu PARI library.

Verify with:

```bash
ldd build/massey | grep pari
```

The library path should point to `$HOME/.local/lib`.

## Regression test

The following test reproduces the expected rank-one triple Massey product
for the imaginary quadratic field of discriminant `-90868`:

```bash
./build/massey 5 "s^2+22717"
```

Expected relevant output:

```text
K cyc: [10, 5]

p-rank: 2
r-rank: 2

Matrix:

[0, 2, 1, 0]
[0, 4, 2, 0]

rank: 1
```

This test was run successfully with patched PARI 2.17.4.
