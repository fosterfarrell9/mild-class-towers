\\ Pure-GP verifier for a single arithmetic certificate.
\\
\\ An independent second implementation of the checks carried out by
\\ verifier/verify_certificate.c.  It needs no patched PARI: the Artin
\\ normalization goes through galoisinit and the definition of the Frobenius
\\ instead of the library-internal rnfcycaut, allauts and cyclicrelfrob.
\\ It reads one certificate at a time; the sweep over all certificates stays
\\ with the C program, which is orders of magnitude faster.
\\
\\ Run from the repository root:
\\   CERT_DIR=certificates/p5/K-2800905-p5 gp -q tools/verify_certificate.gp
\\
\\ Everything stays in the session afterwards: K, and for the last entry
\\ Lrel, Labs, sig, tAC, Iprime, JJ.

default(parisize, 800*10^6);
default(parisizemax, 8*10^9);

certdir = getenv("CERT_DIR");
if (certdir == 0 || #certdir == 0, certdir = "certificates/p5/K-2800905-p5");
certpath = concat(certdir, "/certificate.gp");

failures  = 0;
checks    = 0;
completed = 0;
report(lab, col, what, ok) =
{
  checks = checks + 1;
  if (!ok, failures = failures + 1);
  printf("%-5s e%d  %-36s %s\n", lab, col, what, if(ok, "PASS", "*** FAIL ***"));
}

\\ The principal ideal of a factored element.  Exponents are merged prime by
\\ prime before the product is formed; forming it factor by factor would raise
\\ prime ideals to the stored exponents, which run into the tens of millions.
compact_ideal(nf, cpt) =
{
  my(P = List(), E = List(), fa, hit);
  for (i = 1, #cpt[,1],
    fa = idealfactor(nf, cpt[i,1]);
    for (j = 1, #fa[,1],
      hit = 0;
      for (k = 1, #P,
        if (P[k] == fa[j,1], E[k] = E[k] + fa[j,2]*cpt[i,2]; hit = 1; break()));
      if (!hit, listput(P, fa[j,1]); listput(E, fa[j,2]*cpt[i,2]))));
  idealfactorback(nf, Vec(P), Vec(E));
}

\\ The secondary-norm operator at an arbitrary character x, from the
\\ polarization.  Available once the aggregate stage has run.
Dx(v) =
{
  lift(Mod(v[1]^2*B11 + v[2]^2*B22 + v[3]^2*B33
         + v[1]*v[2]*B12 + v[1]*v[3]*B13 + v[2]*v[3]*B23, p));
}

\\ Position of the word X_i X_m X_k in the 27 coordinates, k fastest.
widx(i, j, k) = ((i*3 + j)*3 + k) + 1;

one_minus_sigma(nf, s, id) = idealdiv(nf, id, nfgaloisapply(nf, s, id));

relative_norm_compact(rnf, base, cpt) =
{
  my(r = cpt);
  for (i = 1, #cpt[,1],
    r[i,1] = nfalgtobasis(base, rnfeltnorm(rnf, rnfeltabstorel(rnf, cpt[i,1]))));
  r;
}

\\ The factored element num/den.  Built entry by entry: matconcat would
\\ flatten the columns, whose entries are themselves coordinate vectors.
compact_quotient(num, den) =
{
  my(n = #num[,1], r = matrix(n+1, 2));
  for (i = 1, n, r[i,1] = num[i,1]; r[i,2] = num[i,2]);
  r[n+1,1] = den; r[n+1,2] = -1;
  r;
}

class_coords(bnf, id, pp) =
{
  my(ex = bnfisprincipal(bnf, id, 0), r = []);
  for (i = 1, #bnf.cyc,
    if (bnf.cyc[i] % pp == 0, r = concat(r, [lift(Mod(ex[i], pp))])));
  r~;
}

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

\\ The automorphism attached to a permutation of galoisinit, written the way
\\ the certificate writes automorphisms: as the image of the field generator
\\ in integral coordinates.
autcol(nf, gal, perm) = nfalgtobasis(nf, galoispermtopol(gal, perm));

\\ The Frobenius of a prime PK of K in Gal(L|K), as an exponent relative to
\\ the fixed generator gen0.  This is the definition, tested on the integral
\\ basis; idealfrobenius cannot be used because PK may be ramified over Q.
frobexp(PK) =
{
  my(PL, mpr, q, tau, ok, res = -1);
  PL = idealfactor(Labs, rnfidealup(Lrel, PK, 1))[1,1];
  mpr = nfmodprinit(Labs, PL);
  q = idealnorm(K, PK);
  tau = ag;
  for (k = 0, p-1,
    ok = 1;
    for (j = 1, #zk,
      if (nfmodpr(Labs, nfgaloisapply(Labs, tau, zk[j]), mpr)
          != nfmodpr(Labs, zk[j], mpr)^q, ok = 0; break()));
    if (ok, res = k; break());
    tau = nfgaloisapply(Labs, gen0, tau));
  res;
}

cert     = eval(concat(readstr(certpath)));
p        = cert[3];
basepol  = cert[4];
basedisc = cert[5];
entries  = cert[7];

print("certificate: ", certpath);
print("p = ", p, "   base = ", basepol, "   disc = ", basedisc);
K = bnfinit(basepol, 1);
if (K.disc != basedisc, error("stored discriminant does not match bnfinit"));
print("Cl(K) = ", K.cyc, "   h = ", K.no);
print("");

\\ The certified base data: everything the certificate asserts about K
\\ before any class field is built.
basedata = cert[6];
{
report("base", 0, "class group invariants", basedata[1] == K.cyc);
report("base", 0, "class number", basedata[2] == K.no);
report("base", 0, "class group generators", basedata[3] == K.gen);
report("base", 0, "torsion unit order is two", basedata[4] == 2 && K.tu[1] == 2);
report("base", 0, "torsion unit generator is -1", basedata[5] == -1 && K.tu[2] == -1);
report("base", 0, "unit rank is zero", #K.fu == 0);
report("base", 0, "integral basis", basedata[6] == K.zk);
rank = 0;
for (i = 1, #K.cyc, if (K.cyc[i] % p == 0, rank = rank + 1));
report("base", 0, "p-class rank is three", rank == 3);
print("");
}
expected_tensor = if (#basedata >= 7, basedata[7], 0);

\\ The six secondary-norm matrices, filled column by column as the entries
\\ are verified.  Slot order follows the character labels.
cols = vector(18);
seen = vector(18);
{Dlabels = if (p == 3, ["D_x1","D_x2","D_x3","D_(x1+x2+x3)","D_(x1+x2)","D_(x1+x3)"],
                       ["D_x1","D_x2","D_x3","D_(x2+x3)","D_(x1+x2)","D_(x1+x3)"]);}

{
for (n = 1, #entries,
  entry     = entries[n];
  if (#entry != 14, error("invalid entry schema"));
  lab       = entry[1];  col    = entry[2];  charvec = entry[3];
  relpol    = entry[4];  abspol = entry[5];  sig     = entry[6];
  aprime    = entry[7];  JJ     = entry[8];  Iprime  = entry[9];
  tAC       = entry[10]; ell    = entry[11]; prid    = entry[12];
  normclass = entry[13]; intbasis = entry[14];

  \\ (1) character label and vector
  slot = character_slot(lab, p);
  if (slot[1] == 0, error(concat("unknown character label ", lab)));
  chi = apply(c -> lift(Mod(c, p)), Vec(charvec));
  report(lab, col, "character label and vector", chi == slot[2] && chi != [0,0,0]);

  \\ (2) field model
  Lrel = rnfinit(K.nf, relpol);
  Labs = nfinit(abspol);
  report(lab, col, "absolute model matches relative", Lrel.polabs == abspol);
  report(lab, col, "relative degree equals p", poldegree(abspol)/poldegree(basepol) == p);
  report(lab, col, "unramified: disc = disc_K^p", Labs.disc == basedisc^p);

  \\ (3) stored integral basis, and the change into this session's basis
  M = matconcat(vector(#intbasis, i, nfalgtobasis(Labs, intbasis[i])));
  report(lab, col, "stored basis integral and unimodular",
         M == round(M) && abs(matdet(M)) == 1);
  sig    = M * sig;
  \\ A fractional ideal is cleared of its denominator before the Hermite form
  \\ is taken, and scaled back afterwards.
  den    = denominator(Iprime);
  Iprime = mathnf(M * (den * Iprime)) / den;
  \\ Rational generators of the factored element carry no basis and are left
  \\ alone; only coordinate vectors are transformed.
  for (i = 1, #tAC[,1],
    if (type(tAC[i,1]) == "t_COL", tAC[i,1] = M * tAC[i,1]));

  \\ (4) sigma fixes K and has order exactly p
  bg = rnfeltup(Lrel, variable(basepol), 1);
  ag = nfalgtobasis(Labs, variable(abspol));
  report(lab, col, "sigma fixes the base field", nfgaloisapply(Labs, sig, bg) == bg);
  cur = ag; ord = 0;
  for (k = 1, p, cur = nfgaloisapply(Labs, sig, cur); if (cur == ag && !ord, ord = k));
  report(lab, col, "sigma has order exactly p", ord == p);

  \\ (5) Artin character and normalization -- no patched PARI here
  gal = galoisinit(Labs);
  fixK = [];
  for (i = 1, #gal.group,
    aut = autcol(Labs, gal, gal.group[i]);
    if (nfgaloisapply(Labs, aut, bg) == bg, fixK = concat(fixK, [aut])));
  report(lab, col, "Gal(L|K) has order p", #fixK == p);
  gen0 = 0;
  for (i = 1, #fixK,
    c = ag; o = 0;
    for (k = 1, p, c = nfgaloisapply(Labs, fixK[i], c); if (c == ag && !o, o = k));
    if (o == p && gen0 == 0, gen0 = fixK[i]));
  zk = vector(#Labs.zk, j, nfalgtobasis(Labs, Labs.zk[j]));
  sigpow = 0; cur = ag;
  for (k = 1, p-1, cur = nfgaloisapply(Labs, gen0, cur); if (cur == sig, sigpow = k));
  report(lab, col, "sigma is a power of the generator", sigpow != 0);

  inv = lift(Mod(sigpow, p)^-1);
  canon = []; norml = [];
  for (i = 1, #K.cyc,
    if (K.cyc[i] % p == 0,
      fa = idealfactor(K, K.gen[i]); ex = 0;
      for (j = 1, #fa[,1],
        d = frobexp(fa[j,1]);
        if (d < 0, error("Frobenius not found among the automorphisms"));
        ex = lift(Mod(ex + d * fa[j,2], p)));
      canon = concat(canon, [ex]);
      norml = concat(norml, [lift(Mod(ex * inv, p))])));
  report(lab, col, "base p-rank is three", #canon == 3);
  piv = 0;
  for (i = 1, 3, if (chi[i] != 0 && !piv, piv = i));
  scale = lift(canon[piv] * Mod(chi[piv], p)^-1);
  report(lab, col, "Artin character cuts out ker(x)",
         scale != 0 && vector(3, i, lift(Mod(scale*chi[i], p))) == canon);
  report(lab, col, "sigma normalization matches x", norml == chi);

  \\ (6) the input class and the relation (a') J^p = O_K
  report(lab, col, "(a') J^p = O_K",
         idealmul(K, idealhnf(K, aprime), idealpow(K, JJ, p)) == matid(2));
  ex = bnfisprincipal(K, JJ, 0); pd = 0; ok = 1;
  for (i = 1, #K.cyc,
    want = 0;
    if (K.cyc[i] % p == 0, pd = pd + 1; if (pd == col, want = K.cyc[i]/p));
    if (Mod(ex[i], K.cyc[i]) != Mod(want, K.cyc[i]), ok = 0));
  report(lab, col, "[J] is the labelled torsion class", ok);

  \\ (7) AC1
  op = one_minus_sigma(Labs, sig, one_minus_sigma(Labs, sig, Iprime));
  ac1 = idealmul(Labs, idealmul(Labs, op, compact_ideal(Labs, tAC)),
                 rnfidealup(Lrel, JJ, 1));
  report(lab, col, "AC1: (1-s)^2 I' (t) i(J) = O_L", ac1 == matid(poldegree(abspol)));

  \\ (8) AC2 and the modular sign
  quo = compact_quotient(relative_norm_compact(Lrel, K, tAC), aprime);
  report(lab, col, "AC2: N(t)/a' generates O_K", compact_ideal(K, quo) == matid(2));
  dec = idealprimedec(K, ell); found = 0;
  for (i = 1, #dec, if (dec[i] == prid, found = 1));
  report(lab, col, "stored prime lies above ell", found && isprime(ell) && ell % 2 == 1);
  mp = nfmodprinit(K, prid); res = nfmodpr(K, 1, mp); unit = 1;
  for (i = 1, #quo[,1],
    if (nfeltval(K, quo[i,1], prid) != 0, unit = 0);
    res = res * nfmodpr(K, quo[i,1], mp)^quo[i,2]);
  report(lab, col, "AC2 factors are units at q", unit);
  report(lab, col, "N(t)/a' is +1 modulo q", res == nfmodpr(K, 1, mp));

  \\ (9) the norm class
  NI = rnfidealnormrel(Lrel, rnfidealabstorel(Lrel, Iprime));
  corr = if (p == 3, idealmul(K, NI, JJ), NI);
  coords = class_coords(K, corr, p);
  report(lab, col, "norm-class coordinates", coords == Col(normclass));
  if (seen[(slot[1]-1)*3 + col], error("duplicate certificate entry"));
  seen[(slot[1]-1)*3 + col] = 1;
  cols[(slot[1]-1)*3 + col] = coords;
  print("");
);
}

\\ ---------------------------------------------------------- aggregate stage
{
full = 1;
for (i = 1, 18, if (!seen[i], full = 0));
report("all", 0, "all eighteen entries present", full);

D = vector(6, i, matconcat(vector(3, j, cols[(i-1)*3 + j])));

\\ Polarization.  B[i][j] is the mixed term of the quadratic form D.
B11 = D[1]; B22 = D[2]; B33 = D[3];
B12 = lift(Mod(D[5] - D[1] - D[2], p));
B13 = lift(Mod(D[6] - D[1] - D[3], p));
B23 = if (p == 3,
          lift(Mod(D[4] - D[1] - D[2] - D[3] - B12 - B13, p)),
          lift(Mod(D[4] - D[2] - D[3], p)));

\\ Consistency of the polarization against the stored matrices themselves.
report("all", 0, "polarization reproduces D_(x1+x2)",
       lift(Mod(B11 + B22 + B12, p)) == D[5]);
report("all", 0, "polarization reproduces D_(x1+x3)",
       lift(Mod(B11 + B33 + B13, p)) == D[6]);
report("all", 0, "homogeneity D_(2x1) = 4 D_x1",
       lift(Mod(4*D[1], p)) == lift(Mod(4*B11, p)));

\\ The 3 x 27 tensor.  The entry at the word X_i X_m X_k is MINUS the
\\ polarization term, and the sign is forced rather than conventional:
\\ D_x[m] = M(x,x,x_m) puts the two characters in the first two slots, the
\\ tensor entry has them on the outside, and the cyclic shuffle identity
\\ together with the outer-reversal symmetry gives
\\     M_imk = -(M_ikm + M_kim) = -DeltaD(x_i,x_k)[m].
\\ The diagonal rule M_imi = -2 D_(x_i)[m] is the same statement, since
\\ DeltaD(x,x) = D_2x - 2 D_x = 2 D_x.  With the natural sign the cyclic
\\ identity fails.  At p = 3 the scalar -2 = 1 is invisible.
tensor = vector(3);
for (rel = 1, 3,
  row = vector(27);
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

\\ Shuffle identities: outer reversal and the cyclic sum.
okrev = 1; okcyc = 1;
for (rel = 1, 3,
  row = tensor[rel];
  for (i = 0, 2, for (j = 0, 2, for (k = 0, 2,
    if (row[widx(i,j,k)] != row[widx(k,j,i)], okrev = 0);
    if (Mod(row[widx(i,j,k)] + row[widx(j,k,i)] + row[widx(k,i,j)], p) != 0,
        okcyc = 0)))));
report("all", 0, "outer-reversal shuffle identity", okrev);
report("all", 0, "cyclic shuffle identity", okcyc);

if (expected_tensor != 0,
  report("all", 0, "embedded expected tensor",
         tensor == vector(3, i, vector(27, j, lift(Mod(expected_tensor[i][j], p)))))
, print("     embedded expected tensor           ABSENT"));

print("");
for (i = 1, 6, print(Dlabels[i], " = ", D[i]));
print("");
print("TENSOR_3_BY_27 = ", tensor);
completed = 1;
}

print("=================================================");
{
if (!completed,
    print("*** INCOMPLETE: the run stopped after ", checks, " checks ***"),
  if (failures == 0,
      print("CERTIFICATE VERIFIED   (", checks, " checks, ", #entries, " entries)"),
      print("*** ", failures, " of ", checks, " CHECKS FAILED ***")));
}
print("");
print("session holds:  D  (the six matrices),  Dlabels,  tensor,");
print("                B11 B22 B33 B12 B13 B23  (the polarization),");
print("                K,  and for the last entry Lrel, Labs, sig, tAC, Iprime, JJ.");
print("Evaluate D_x for any x:   Dx([1,2,3])");
