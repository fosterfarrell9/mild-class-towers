#!/usr/bin/env bash
# Build the patched PARI 2.17.4 that Massey-pari requires.
#
# Stock PARI declares rnfcycaut, allauts and cyclicrelfrob static, so they
# are not exported by libpari; this script applies the one-file patch of
# doc/pari-2.17.4-patch.md and installs the result into a private prefix.
# Nothing outside that prefix and ~/src is touched, and no PARI source is
# redistributed with this repository.
#
# Usage:
#   tools/build-patched-pari.sh [PREFIX]        # default: $HOME/.local
#
# Environment:
#   PARI_MT=pthread    build with the multithreading engine (rarely useful:
#                      measured 14% slower on relation-collection-bound
#                      fields, and each thread needs its own stack)
#   PARI_SRC=<dir>     where to unpack sources (default: $HOME/src)
#   PARI_GMP=<dir>     prefix of a GMP outside the system paths; a system
#                      GMP is found without this, multiarch layout included
#   CFLAGS=<flags>     default -march=native, which ties the build to the
#                      instruction set of the machine it is built on; drop
#                      it for a binary that runs on older hardware too
#
# Afterwards, from the repository root:
#   make PARI=<PREFIX>
#   make -C verifier PARI=<PREFIX>
set -euo pipefail

PREFIX="${1:-$HOME/.local}"
SRC="${PARI_SRC:-$HOME/src}"
VERSION=2.17.4
TARBALL="pari-${VERSION}.tar.gz"
# Pinned upstream digest.  The patch below rewrites three declarations by
# exact string match, so it is worth knowing that the source tree is the
# one those patterns were read from rather than a truncated download or a
# mirror serving a different file.
SHA256=02651d99c391007d384b3fadbc20abc6916b77036f9e496c99e9ce8688ca4b53

echo "=== prerequisites"
missing=""
for tool in gcc make wget python3; do
    command -v "$tool" >/dev/null || missing="$missing $tool"
done
# Ask the compiler rather than guessing a path: on Debian and Ubuntu the
# header is architecture-dependent and lives under /usr/include/<triplet>,
# so a test for /usr/include/gmp.h reports it missing when it is installed.
echo '#include <gmp.h>' | cc -E ${PARI_GMP:+-I"$PARI_GMP/include"} - \
    >/dev/null 2>&1 || missing="$missing libgmp-dev"
if [ -n "$missing" ]; then
    echo "missing:$missing" >&2
    echo "on Debian/Ubuntu: sudo apt-get install -y build-essential wget libgmp-dev" >&2
    exit 1
fi

echo "=== fetching PARI $VERSION"
mkdir -p "$SRC" && cd "$SRC"
if [ ! -f "$TARBALL" ]; then
    wget -q "https://pari.math.u-bordeaux.fr/pub/pari/unix/$TARBALL" \
      || wget -q "https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.17/$TARBALL"
fi
# Checked on every run, not only after a fetch: the copy in $SRC may be
# left over from an earlier run that was interrupted mid-download.
digest=$(sha256sum "$TARBALL" | cut -d' ' -f1)
if [ "$digest" != "$SHA256" ]; then
    echo "checksum mismatch for $SRC/$TARBALL" >&2
    echo "  expected $SHA256" >&2
    echo "  found    $digest" >&2
    echo "delete the file and rerun to fetch it again" >&2
    exit 1
fi
echo "checksum ok ($SHA256)"
rm -rf "pari-$VERSION"
tar xf "$TARBALL"
cd "pari-$VERSION"

echo "=== applying the three-function patch"
cp src/modules/algebras.c src/modules/algebras.c.orig
python3 - <<'PATCH'
from pathlib import Path
path = Path("src/modules/algebras.c")
text = path.read_text()
replacements = [
    ("static long\ncyclicrelfrob(GEN rnf, GEN auts, GEN pr)",
     "long\ncyclicrelfrob(GEN rnf, GEN auts, GEN pr)"),
    ("static GEN\nallauts(GEN rnf, GEN aut)",
     "GEN\nallauts(GEN rnf, GEN aut)"),
    ("static GEN\nrnfcycaut(GEN rnf)",
     "GEN\nrnfcycaut(GEN rnf)"),
]
for old, new in replacements:
    if text.count(old) != 1:
        raise SystemExit(f"pattern not found exactly once: {old!r}")
    text = text.replace(old, new)
path.write_text(text)
print("patched cyclicrelfrob, allauts, rnfcycaut")
PATCH

echo "=== configuring and building"
configure_args=(--prefix="$PREFIX")
[ -n "${PARI_MT:-}" ] && configure_args+=(--mt="$PARI_MT")
# Only needed for a GMP outside the system paths; PARI's Configure finds a
# system GMP by itself, including under the Debian multiarch layout.
gmp_dir="${PARI_GMP:-$PREFIX}"
[ -f "$gmp_dir/include/gmp.h" ] && configure_args+=(--with-gmp="$gmp_dir")
CFLAGS="${CFLAGS:--march=native}" ./Configure "${configure_args[@]}"
nice -n 19 make -j"$(nproc)" gp
make install

echo
echo "=== verification"
grep -E "^(CFLAGS|gmp|thread_engine|asmarch)" "$PREFIX/lib/pari/pari.cfg"
found=$(nm -D --defined-only "$PREFIX"/lib/libpari*.so \
        | grep -cE ' (rnfcycaut|allauts|cyclicrelfrob)$' || true)
echo "exported patch symbols: $found of 3"
[ "$found" -eq 3 ] || { echo "patch symbols missing" >&2; exit 1; }
echo
echo "PARI $VERSION installed in $PREFIX"
echo "next: make PARI=$PREFIX  &&  make -C verifier PARI=$PREFIX"
