\\ Zeugen-Export-Orakel: rechnet wie orakel-capped-red.gp, exportiert aber
\\ pro Standardcharakter und Torsionsspalte den vollen Normzeugen (t, I')
\\ und schreibt ein certificate.gp im Format der Anhang-E-Schicht
\\ (Standardfamilie e_i, e_i+e_j; bei d=4 also 10 x 4 = 40 Eintraege).
\\
\\   POL="y^2 - y + 3036244875" P=3 gp -q zeugen-orakel.gp
\\
\\ Die einzige unbedingt gemachte Klassengruppe ist die von K (bnfcertify);
\\ alle Teilkoerperdaten beweisen sich im Verifier durch Idealarithmetik.

default(parisize, 6*10^9);
default(parisizemax, 20*10^9);

basepol = subst(eval(getenv("POL")), variable(eval(getenv("POL"))), y);
p = eval(getenv("P"));
print("Zeugen-Orakel   p = ", p, "   Grundkoerper ", basepol);
gettime();
K = bnfinit(basepol, 1);
print("bnfinit(K)  ", gettime(), " ms   cyc = ", K.cyc, "   h = ", K.no);
{if (bnfcertify(K) != 1, error("bnfcertify(K) fehlgeschlagen"));}
print("bnfcertify(K)  ", gettime(), " ms   UNBEDINGT");

tor = [];
{for (j = 1, #K.cyc, if (K.cyc[j] % p == 0, tor = concat(tor, [j])));}
rk = #tor;
print("p-Rang = ", rk);

\\ --- J_n und a'_n:  (a'_n) J_n^p = O_K
JJ = vector(rk); AP = vector(rk);
{
for (n = 1, rk,
  e = vector(#K.cyc);
  e[tor[n]] = K.cyc[tor[n]]/p;
  JJ[n] = idealred(K, idealfactorback(K, K.gen, e));
  pr = bnfisprincipal(K, idealpow(K, JJ[n], p), 3);
  if (pr[1] != vector(#K.cyc)~, error("J^p nicht prinzipal?"));
  AP[n] = nfeltdiv(K, 1, pr[2]);
  if (idealmul(K, idealhnf(K, AP[n]), idealpow(K, JJ[n], p)) != matid(poldegree(basepol)),
      error("(a') J^p != O_K")));
}
print("J_n und a'_n stehen");

\\ --- die Standardfamilie
need = List(); needlab = List();
{
for (i = 1, rk, v = vector(rk); v[i] = 1;
  listput(need, v); listput(needlab, concat("x", Str(i))));
for (i = 1, rk-1, for (j = i+1, rk,
  v = vector(rk); v[i] = 1; v[j] = 1;
  listput(need, v); listput(needlab, concat(["x", Str(i), "+x", Str(j)]))));
}
done = vector(#need);
entries = List();
Dmat = vector(#need);

bnr = bnrinit(K, 1);
subs = subgrouplist(bnr, [p]);
print("Untergruppen vom Index ", p, ": ", #subs);

frobexp(PK) =
{
  my(PL, mpr, q, tau, ok, res = -1);
  PL = idealfactor(nfL, rnfidealup(Lrel, PK, 1))[1,1];
  mpr = nfmodprinit(nfL, PL);
  q = idealnorm(K, PK);
  tau = ag;
  for (k = 0, p-1,
    ok = 1;
    for (j = 1, #zk,
      if (nfmodpr(nfL, nfgaloisapply(nfL, tau, zk[j]), mpr)
          != nfmodpr(nfL, zk[j], mpr)^q, ok = 0; break()));
    if (ok, res = k; break());
    tau = nfgaloisapply(nfL, s0, tau));
  res;
}

ccoords(id) = {vector(rk, n, lift(Mod(bnfisprincipal(K, id, 0)[tor[n]], p)))~;}

one_minus_sigma(nf, s, id) = idealdiv(nf, id, nfgaloisapply(nf, s, id));

\\ Relative Norm eines kompakten Elements, Zeile fuer Zeile
norm_compact(rnf, base, cpt) =
{
  my(r = cpt);
  for (i = 1, #cpt[,1],
    r[i,1] = nfalgtobasis(base, rnfeltnorm(rnf, rnfeltabstorel(rnf, cpt[i,1]))));
  r;
}

\\ --- Hauptschleife
{
for (s = 1, #subs,
  if (vecsum(done) == #need, break());
  gettime();
  relpol = bnrclassfield(bnr, subs[s])[1];
  Lrel = rnfinit(K.nf, relpol);
  Labs = Lrel.polabs;
  L = bnfinit(Labs, 1);
  nfL = L.nf;
  tbnf = gettime();

  gal = galoisinit(nfL);
  ag = nfalgtobasis(nfL, variable(Labs));
  bg = rnfeltup(Lrel, y, 1);
  fixK = [];
  for (i = 1, #gal.group,
    a = nfalgtobasis(nfL, galoispermtopol(gal, gal.group[i]));
    if (nfgaloisapply(nfL, a, bg) == bg, fixK = concat(fixK, [a])));
  s0 = 0;
  for (i = 1, #fixK,
    c = ag; o = 0;
    for (k = 1, p, c = nfgaloisapply(nfL, fixK[i], c); if (c == ag && !o, o = k));
    if (o == p && s0 == 0, s0 = fixK[i]));

  zk = vector(#nfL.zk, j, nfalgtobasis(nfL, nfL.zk[j]));
  canon = vector(rk);
  for (n = 1, rk,
    fa = idealfactor(K, K.gen[tor[n]]); ex = 0;
    for (j = 1, #fa[,1], ex = ex + frobexp(fa[j,1]) * fa[j,2]);
    canon[n] = lift(Mod(ex, p)));

  cycL = L.cyc;
  S = matconcat(vector(#cycL, i,
        bnfisprincipal(L, nfgaloisapply(nfL, s0, L.gen[i]), 0)));

  for (k = 1, p-1,
    chi = vector(rk, i, lift(Mod(canon[i], p) * Mod(k, p)^-1));
    idx = 0;
    for (t = 1, #need, if (need[t] == chi && !idx, idx = t));
    if (!idx || done[idx], next());
    \\ sigma = s0^k als Erzeugerbild
    sig = ag;
    for (j = 1, k, sig = nfgaloisapply(nfL, s0, sig));
    Mk = matid(#cycL);
    for (j = 1, k, Mk = Mk * S);
    M2 = (matid(#cycL) - Mk)^2;

    cols = vector(rk);
    for (n = 1, rk,
      iJ = rnfidealup(Lrel, JJ[n], 1);
      target = -bnfisprincipal(L, iJ, 0);
      X = matsolvemod(M2, cycL~, target);
      if (type(X) == "t_INT", error("keine H90-Loesung  s=", s, " n=", n));
      Ip = idealhnf(L, 1);
      for (ii = 1, #cycL,
        if (X[ii] != 0,
          Ip = idealmul(L, Ip, idealpow(L, L.gen[ii], X[ii], 1), 1)));
      Ip = idealhnf(L, Ip);

      \\ Zeuge t:  (t0) = (1-sigma)^2 Ip * iJ  (prinzipal),  t = 1/t0
      op = one_minus_sigma(nfL, sig, one_minus_sigma(nfL, sig, Ip));
      pid = idealmul(L, op, iJ);
      pr = bnfisprincipal(L, pid, 7);
      if (pr[1] != vector(#cycL)~, error("AC1-Ideal nicht prinzipal  s=", s, " n=", n));
      tc = pr[2];
      if (type(tc) != "t_MAT", tc = Mat([tc, 1]));
      for (i = 1, #tc[,1], tc[i,2] = -tc[i,2]);

      \\ Vorzeichen: N(t)/a' als Produkt an einer Primstelle q auswerten
      quo = norm_compact(Lrel, K, tc);
      quo = matconcat([quo; Mat([nfalgtobasis(K, AP[n]), -1])]);
      ell = 0; prid = 0;
      forprime (q = 3, 10^6,
        if (q == p || (poldisc(basepol) % q) == 0, next());
        my(dec = idealprimedec(K, q), pi = dec[1], good = 1);
        for (i = 1, #quo[,1],
          if (nfeltval(K, quo[i,1], pi) != 0, good = 0; break()));
        if (good, ell = q; prid = pi; break()));
      if (!ell, error("keine Vorzeichen-Primstelle gefunden"));
      mp = nfmodprinit(K, prid);
      res = nfmodpr(K, 1, mp);
      for (i = 1, #quo[,1], res = res * nfmodpr(K, quo[i,1], mp)^quo[i,2]);
      if (res != nfmodpr(K, 1, mp),
        if (res == nfmodpr(K, -1, mp),
          tc = matconcat([tc; Mat([nfalgtobasis(nfL, -1), 1])]),
          error("N(t)/a' ist keine Einheit +-1  s=", s, " n=", n)));

      \\ Normklasse
      NI = rnfidealnormrel(Lrel, rnfidealabstorel(Lrel, Ip));
      corr = if (p == 3, idealmul(K, NI, JJ[n]), NI);
      cols[n] = ccoords(corr);

      listput(entries, [needlab[idx], n, chi~, relpol, Labs, sig,
                        AP[n], idealhnf(K, JJ[n]), Ip, tc, ell, prid,
                        cols[n], nfL.zk]));
    Dmat[idx] = matconcat(cols);
    done[idx] = 1;
    printf("  s=%2d k=%d  %-8s  bnfinit %5d ms  gesamt %6d ms\n",
           s, k, needlab[idx], tbnf, gettime() + tbnf));
);
}
{if (vecsum(done) != #need, error("nicht alle Charaktere erreicht"));}
print("alle ", #need, " Charaktere mit Zeugen exportiert");

\\ --- D-matrizen.txt (Diff-Anker gegen die Ankuendigung)
{
out = "D-matrizen.txt";
system(concat(["rm -f ", out]));
write(out, "\\\\ D_x fuer die Standardcharaktere, p = ", p, ", Rang ", rk);
for (t = 1, #need,
  write(out, "[", need[t], ", ", lift(Mod(Dmat[t], p)), "]"));
}

\\ --- erwarteter Tensor (rk x rk^3), Vorzeichen wie im d=3-Verifier
widx(i, j, k) = ((i*rk + j)*rk + k) + 1;
{
B = matrix(rk, rk);
for (i = 1, rk, B[i,i] = 0);
for (i = 1, rk-1, for (j = i+1, rk,
  pos = 0;
  for (t = 1, #need, if (need[t][i] == 1 && need[t][j] == 1 && vecsum(need[t]) == 2, pos = t));
  Bij = lift(Mod(Dmat[pos] - Dmat[i] - Dmat[j], p));
  B[i,j] = pos; B[j,i] = pos));
tensor = vector(rk);
for (rel = 1, rk,
  row = vector(rk^3);
  for (i = 0, rk-1,
    for (m = 0, rk-1,
      row[widx(i, m, i)] = lift(Mod(-2 * Dmat[i+1][m+1, rel], p))));
  for (i = 0, rk-2, for (kk = i+1, rk-1,
    Bij = lift(Mod(Dmat[B[i+1,kk+1]] - Dmat[i+1] - Dmat[kk+1], p));
    for (m = 0, rk-1,
      v = lift(Mod(-Bij[m+1, rel], p));
      row[widx(i, m, kk)] = v;
      row[widx(kk, m, i)] = v)));
  tensor[rel] = row);
}

\\ --- certificate.gp schreiben
{
basedata = [K.cyc, K.no, K.gen, K.tu[1], lift(K.tu[2]), K.zk, tensor];
cert = [2, 0, p, basepol, K.disc, basedata, Vec(entries)];
system("rm -f certificate.gp");
write("certificate.gp", cert);
}
print("geschrieben: certificate.gp (", #entries, " Eintraege)");
quit
