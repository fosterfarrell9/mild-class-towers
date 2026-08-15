# Mosunov--Jacobson imaginary-quadratic class-group tables

## Purpose and provenance

A. S. Mosunov and M. J. Jacobson Jr., *Unconditional Class Group
Tabulation of Imaginary Quadratic Fields to \(|\Delta|<2^{40}\)*,
Math. Comp. 85 (2016), 1983--2009, computed the class groups of all
imaginary quadratic fields in that range unconditionally.  LMFDB distributes
the raw tables separately from its ordinary searchable number-field database:

<https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups>

LMFDB data is licensed CC-BY-SA.  Cite both the paper and LMFDB when using the
tables.

These tables are used here only to discover candidate fields.  Every mildness
claim still depends on this repository's own exact, audited secondary-norm and
Massey-product pipeline.

## Block-zero validation files

The four files used for the independent 25-million comparison were retrieved
on 2026-07-25 through the individual-download form at the LMFDB URL above.
The corresponding query has `Fetch=fetch`, the filename base shown below,
and `k=0`.  Raw files are stored outside this repository.

| file | compressed bytes | uncompressed bytes | SHA-256 |
|---|---:|---:|---|
| `cl3mod8.0.gz` | 113740185 | 355146018 | `3abc1a9044e7e6f9bf6c52c94ad16325a05435fa6b76fbaca96587d03fa3914b` |
| `cl7mod8.0.gz` | 135449175 | 377294247 | `1e37f94365e6999a12840f96e80cb7a6aebe2078cc06c4d08d7fb724a8c0c61b` |
| `cl4mod16.0.gz` | 56849892 | 199066873 | `09c7a0b2683321fa7946413cb313928925c3542d8a08155f82eb9c4d0576b4fe` |
| `cl8mod16.0.gz` | 57165418 | 198158390 | `d2ef3bab41f24baeb01454fda34279fdd73dc39f7665b67a0313d86d73a3ab2d` |

For example, the first file's exact source URL is:

```text
https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups?Fetch=fetch&filenamebase=cl3mod8&k=0
```

Replace `cl3mod8` by the other filename base for the remaining files.

## Raw organization and format

Raw files must remain outside the Git repository.  There are four families,
each with blocks `k=0,...,4095`:

| filename base | magnitude residue |
|---|---|
| `cl3mod8` | `3 mod 8` |
| `cl7mod8` | `7 mod 8` |
| `cl4mod16` | `4 mod 16` |
| `cl8mod16` | `8 mod 16` |

File `cl{r}mod{m}.{k}.gz` covers
`k*2^28 <= |D| < (k+1)*2^28`.  Initialize the signed accumulator

```text
d_0 = -k*2^28-r.
```

Each decompressed line is

```text
a h c_1 c_2 ... c_t
```

and decodes by

```text
d_i = d_(i-1)-m*a
h(d_i) = h
Cl(d_i) = [c_1,c_2,...,c_t]
h = product(c_j).
```

Thus `a` is a delta, including on the first line.  The lines are ordered by
increasing absolute discriminant.  The invariant vector is preserved exactly;
in particular, the raw files represent a trivial group as `[1]`.

The published first lines of `cl4mod16.1.gz`,

```text
0 12160 380 4 4 2
2 4392 2196 2
```

decode to discriminants `-268435460` and `-268435492`.

## Parser

`tools/mj_classgroup_tables.py` streams gzip input without permanently
decompressing it.  It validates filenames, block bounds, residue classes,
integer fields, positivity, and the class-number/invariant product.  Candidate
records alone are retained and sorted for deterministic TSV output.

The default filter is `p=5`, `p-rank >= 3`.  Here p-rank is the number of
cyclic invariant factors divisible by p.

```sh
python3 tools/mj_classgroup_tables.py \
  --prime 5 --min-p-rank 3 \
  --min-abs-disc 1 --max-abs-disc 25000000 \
  --output /tmp/mj-p5-r3.tsv \
  /path/to/cl3mod8.0.gz /path/to/cl7mod8.0.gz \
  /path/to/cl4mod16.0.gz /path/to/cl8mod16.0.gz
```

TSV columns are:

```text
D_K
class_group_invariants
class_number
p
p_rank
source_filename
source_block
source_record
```

Run decoder tests alone, or include the four block-zero files for the
25-million acceptance test:

```sh
sh tests/test_mj_classgroup_tables.sh
sh tests/test_mj_classgroup_tables.sh /path/to/raw/files
```

The real-data test compares the imported candidates programmatically with
`census/p5-r3/scan-1-25000000.gp` and independently verifies both known
class groups using PARI `quadclassunit`.

Run the deterministic 20-record PARI audit with:

```sh
python3 tools/validate_mj_tables.py /path/to/raw/files
```

It selects the first record at or above each of five fixed magnitude targets
in every one of the four block-zero files.  This gives reproducible coverage
across residue families and discriminant ranges.

The validated output for `p=5`, p-rank at least 3, and
`1 <= |D| <= 25000000` is committed as
`p5-rank-ge3-1-25000000.tsv`.  It agrees exactly with the independent PARI
scanner artifact and contains two records.

## Scaling

LMFDB documents 16,384 files of roughly 50--200 MB each; the collection was
reported as approximately 2.1 TB.  The partition contains no class-group or
p-rank index, so no block can be skipped for a p-rank query.

A full local copy should be avoided.  The practical design is a resumable
stream-and-discard driver: fetch one file, record its URL/size/SHA-256 or ETag,
parse and retain candidates, atomically mark the manifest entry complete, then
discard the raw file.  Such a network driver is intentionally not included
until stable direct URLs and authoritative checksums or ETags are available.
