#!/bin/sh
set -eu

binary=${1:-./build/massey}
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT

"$binary" --candidate-self-test

"$binary" --scan-candidates \
  --prime 5 --rank 3 \
  --min-abs-disc 11203618 --max-abs-disc 11203622 \
  --output "$test_dir/known1.gp" --checkpoint-every 2 \
  --progress-seconds 1 >"$test_dir/progress.log"
grep '^SCAN_START ' "$test_dir/progress.log"
grep '^SCAN_PROGRESS ' "$test_dir/progress.log"
grep '^CANDIDATE ' "$test_dir/progress.log"
"$binary" --scan-candidates \
  --prime 5 --rank 3 \
  --min-abs-disc 18397405 --max-abs-disc 18397409 \
  --output "$test_dir/known2.gp" --checkpoint-every 2

"$binary" --scan-candidates \
  --prime 5 --rank 1 --min-abs-disc 1 --max-abs-disc 500 \
  --output "$test_dir/full.gp" --checkpoint-every 73
"$binary" --scan-candidates \
  --prime 5 --rank 1 --min-abs-disc 1 --max-abs-disc 500 \
  --output "$test_dir/resumed.gp" --checkpoint-every 73 --stop-after 137
"$binary" --scan-candidates \
  --prime 5 --rank 1 --min-abs-disc 1 --max-abs-disc 500 \
  --output "$test_dir/resumed.gp" --checkpoint-every 73 --resume

cmp "$test_dir/full.gp" "$test_dir/resumed.gp"

gp -q <<EOF
r1=read("$test_dir/known1.gp");
r2=read("$test_dir/known2.gp");
if(r1[2][2]!="COMPLETE" || r2[2][2]!="COMPLETE", error("scan incomplete"));
if(#r1[10][2]!=1 || r1[10][2][1][1]!=-11203620, error("first known field missing"));
if(r1[10][2][1][3]!=[10,10,10] || r1[10][2][1][5]!=3, error("first class data"));
if(#r2[10][2]!=1 || r2[10][2][1][1]!=-18397407, error("second known field missing"));
if(r2[10][2][1][5]!=3, error("second p-rank"));
rr=read("$test_dir/resumed.gp");
for(i=2,#rr[10][2],if(abs(rr[10][2][i-1][1])>=abs(rr[10][2][i][1]),error("candidate order")));
print("CANDIDATE_SCANNER_INTEGRATION_TEST PASS");
quit
EOF
