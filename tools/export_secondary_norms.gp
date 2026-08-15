\\ Export the six verified secondary-norm matrices of the principal
\\ example from its arithmetic certificate into a small result-style
\\ data file, so that finite-algebra tools can read all nine fields
\\ uniformly (the other eight fields store the matrices in their
\\ committed result.gp records).
\\
\\ Run from the repository root:
\\   gp -qf tools/export_secondary_norms.gp
\\ writes certificates/p5/K-2800905-p5/secondary-norms.gp deterministically.

certpath = "certificates/p5/K-2800905-p5/certificate.gp";
outpath  = "certificates/p5/K-2800905-p5/secondary-norms.gp";

readcert() =
{
  my(rl = readstr(certpath), t = "");
  for (i = 1, #rl, t = concat(t, rl[i]));
  eval(t);
}

{
my(cert = readcert(), entries = cert[7],
   labels = ["a","b","c","a+b","a+c","b+c"], mats = vector(6));
for (q = 1, 6,
  my(cols = vector(3));
  for (i = 1, #entries,
    if (entries[i][1] == labels[q],
        cols[entries[i][2]] = Col(entries[i][13])));
  mats[q] = matconcat(cols));
system(concat(["rm -f '", outpath, "'"]));
write(outpath, Str([["base_discriminant", cert[5]],
                    ["secondary_norm_samples", mats]]));
print("written: ", outpath);
}
quit
