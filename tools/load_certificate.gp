\\ Pure-GP loader for a single arithmetic certificate.
\\
\\ It reads a certificate and leaves its data in the session, under the same
\\ names as tools/verify_certificate.gp, but performs no checks: no base field
\\ is built, no class field is reconstructed, no identity is tested.  Use it
\\ to inspect or reuse the stored data of a certificate; use
\\ verify_certificate.gp when you need the data proved.
\\
\\ Run from the repository root:
\\   CERT_DIR=certificates/p5/K-2800905-p5 gp -q tools/load_certificate.gp
\\
\\ Afterwards the session holds:
\\   p, basepol, basedisc                     the field and the prime
\\   Clcyc, h, clgen, torsion, zk             the stored base data of K
\\   entries                                  the raw certificate entries
\\   D (the six matrices), Dlabels            the secondary norm operators
\\   B11 B22 B33 B12 B13 B23                  the polarization
\\   tensor                                   the 3 x 27 cubic relation tensor
\\ and the functions
\\   Dx(v)      the secondary norm operator at an arbitrary character x
\\   ent(n)     the n-th entry as a named vector [see labels printed below]
\\   Kbnf()     builds and caches bnfinit(basepol) on demand (the slow step)

default(parisize, 800*10^6);
default(parisizemax, 8*10^9);

certdir = getenv("CERT_DIR");
if (certdir == 0 || #certdir == 0, certdir = "certificates/p5/K-2800905-p5");
certpath = concat(certdir, "/certificate.gp");

\\ Position of the word X_i X_m X_k in the 27 coordinates, k fastest.
widx(i, j, k) = ((i*3 + j)*3 + k) + 1;

\\ Character label -> slot 1..6 and coordinate vector, as in the certificate.
character_slot(lab, pp) =
{
  my(three = (pp == 3), labs, co, idx = 0);
  labs = if (three, ["x1","x2","x3","x1+x2+x3","x1+x2","x1+x3"],
                    ["x1","x2","x3","x2+x3","x1+x2","x1+x3"]);
  co   = if (three, [[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,1,0],[1,0,1]],
                    [[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1]]);
  for (i = 1, 6,
    if (lab == labs[i] || (!three && lab == ["a","b","c","b+c","a+b","a+c"][i]),
        idx = i; break()));
  if (idx, [idx, co[idx]], [0, 0]);
}

\\ The secondary-norm operator at an arbitrary character x, from the
\\ polarization.
Dx(v) =
{
  lift(Mod(v[1]^2*B11 + v[2]^2*B22 + v[3]^2*B33
         + v[1]*v[2]*B12 + v[1]*v[3]*B13 + v[2]*v[3]*B23, p));
}

\\ The n-th entry, with its fields named.  Returns a two-column matrix of
\\ [name, value] pairs for readable inspection.
ent(n) =
{
  my(e = entries[n], names);
  names = ["label","column","character","rel. polynomial","abs. polynomial",
           "sigma","a'","J","I'","t_AC","ell","prime","norm class",
           "integral basis"];
  matconcat([Col(names), Col(vector(#e, i, e[i]))]);
}

\\ bnfinit is the one slow step; build it only when asked, and cache it.
K = 0;
Kbnf() = { if (K == 0, K = bnfinit(basepol, 1)); K; }

cert     = eval(concat(readstr(certpath)));
p        = cert[3];
basepol  = cert[4];
basedisc = cert[5];
basedata = cert[6];
entries  = cert[7];

Clcyc   = basedata[1];
h       = basedata[2];
clgen   = basedata[3];
torsion = [basedata[4], basedata[5]];
zk      = basedata[6];
expected_tensor = if (#basedata >= 7, basedata[7], 0);

print("certificate: ", certpath);
print("p = ", p, "   base = ", basepol, "   disc = ", basedisc);
print("Cl(K) = ", Clcyc, "   h = ", h);
rank = 0;
for (i = 1, #Clcyc, if (Clcyc[i] % p == 0, rank = rank + 1));
print("p-class rank = ", rank, "   entries = ", #entries);
print("");

\\ The six secondary-norm matrices, from the stored norm classes.  No
\\ recomputation: column `col` of slot `slot` is the stored norm-class vector.
cols = vector(18);
seen = vector(18);
{Dlabels = if (p == 3, ["D_x1","D_x2","D_x3","D_(x1+x2+x3)","D_(x1+x2)","D_(x1+x3)"],
                       ["D_x1","D_x2","D_x3","D_(x2+x3)","D_(x1+x2)","D_(x1+x3)"]);}
{
for (n = 1, #entries,
  my(e = entries[n], lab = e[1], col = e[2], slot);
  slot = character_slot(lab, p);
  if (slot[1] == 0, error(concat("unknown character label ", lab)));
  if (col < 1 || col > 3, error("column out of range"));
  cols[(slot[1]-1)*3 + col] = Col(e[13]);
  seen[(slot[1]-1)*3 + col] = 1);
}

full = 1;
for (i = 1, 18, if (!seen[i], full = 0));
if (!full, print("*** WARNING: not all eighteen entries present; D is incomplete ***"));

D = vector(6, i, matconcat(vector(3, j, cols[(i-1)*3 + j])));

\\ Polarization.  B[i][j] is the mixed term of the quadratic family D.
{
B11 = D[1]; B22 = D[2]; B33 = D[3];
B12 = lift(Mod(D[5] - D[1] - D[2], p));
B13 = lift(Mod(D[6] - D[1] - D[3], p));
B23 = if (p == 3,
          lift(Mod(D[4] - D[1] - D[2] - D[3] - B12 - B13, p)),
          lift(Mod(D[4] - D[2] - D[3], p)));
}

\\ The 3 x 27 tensor from the polarization, with the signs of the verifier.
tensor = vector(3);
{
for (rel = 1, 3,
  my(row = vector(27), contr, prs, i, k, v);
  for (i = 0, 2,
    for (m = 0, 2,
      row[widx(i, m, i)] = lift(Mod(-2 * D[i+1][m+1, rel], p))));
  contr = [lift(Mod(-B12, p)), lift(Mod(-B13, p)), lift(Mod(-B23, p))];
  prs = [[0,1], [0,2], [1,2]];
  for (t = 1, 3,
    for (m = 0, 2,
      i = prs[t][1]; k = prs[t][2];
      v = lift(Mod(contr[t][m+1, rel], p));
      row[widx(i, m, k)] = v;
      row[widx(k, m, i)] = v));
  tensor[rel] = row);
}

for (i = 1, 6, print(Dlabels[i], " = ", D[i]));
print("");
print("TENSOR_3_BY_27 = ", tensor);
{
if (expected_tensor != 0,
  print("embedded expected tensor ",
        if (tensor == vector(3, i, vector(27, j, lift(Mod(expected_tensor[i][j], p)))),
            "matches", "DIFFERS")));
}
print("");
print("loaded (no checks performed).  Session holds:");
print("  p, basepol, basedisc, Clcyc, h, clgen, torsion, zk, entries,");
print("  D, Dlabels, B11 B22 B33 B12 B13 B23, tensor.");
print("Functions:  Dx([1,2,3])   ent(n)   Kbnf()");
