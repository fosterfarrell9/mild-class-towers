#!/bin/sh
set -eu

python3 tests/test_mj_classgroup_tables.py

if [ "$#" -eq 0 ]; then
    exit 0
fi
if [ "$#" -ne 1 ]; then
    echo "usage: $0 [raw-table-directory]" >&2
    exit 2
fi

data_dir=$1
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

python3 tools/mj_classgroup_tables.py \
    --prime 5 --min-p-rank 3 \
    --min-abs-disc 1 --max-abs-disc 25000000 \
    --output "$test_dir/imported.tsv" \
    "$data_dir/cl3mod8.0.gz" \
    "$data_dir/cl7mod8.0.gz" \
    "$data_dir/cl4mod16.0.gz" \
    "$data_dir/cl8mod16.0.gz"

sed -n '2,$p' "$test_dir/imported.tsv" |
    cut -f1-2 >"$test_dir/imported-key.tsv"
printf '%s\n' \
    '-11203620	[10,10,10]' \
    '-18397407	[40,10,5]' >"$test_dir/expected-key.tsv"
cmp "$test_dir/expected-key.tsv" "$test_dir/imported-key.tsv"

gp -q <<EOF
r=read("candidates/p5-r3/scan-1-25000000.gp");
c=r[10][2];
if(#c!=2,error("production scan candidate count"));
if(c[1][1]!=-11203620 || c[1][3]!=[10,10,10],error("first production candidate"));
if(c[2][1]!=-18397407 || c[2][3]!=[40,10,5],error("second production candidate"));
q1=quadclassunit(-11203620);
q2=quadclassunit(-18397407);
if(q1.no!=1000 || q1.cyc!=[10,10,10],error("PARI first candidate"));
if(q2.no!=2000 || q2.cyc!=[40,10,5],error("PARI second candidate"));
print("MJ_25M_ACCEPTANCE_TEST PASS");
quit
EOF
