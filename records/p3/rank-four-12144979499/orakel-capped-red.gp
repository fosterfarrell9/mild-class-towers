\\ Ein Orakel in reinem gp: berechnet die D_x-Matrizen von Grund auf.
\\ Rangunabhaengig geschrieben; validiert wird gegen ein gespeichertes Zertifikat.
\\
\\   CERT_DIR=certificates/p3/000/K-3321607-p3 gp -q orakel.gp
\\
\\ Der Grundkoerper lebt in y, damit bnrclassfield x benutzen darf.

default(parisize, 6*10^9);
default(parisizemax, 20*10^9);

\\ Zwei Betriebsarten: gegen ein Zertifikat validieren (CERT_DIR), oder
\\ einen Koerper ohne Vergleichswert rechnen (POL und P).
certdir = getenv("CERT_DIR");
polenv  = getenv("POL");
haveCert = (certdir != 0 && #certdir > 0);
{
if (haveCert,
  cert = eval(concat(readstr(concat(certdir, "/certificate.gp"))));
  p = cert[3];
  basepol = subst(cert[4], variable(cert[4]), y);
  print("Zertifikat: ", certdir)
,
  if (polenv == 0 || #polenv == 0,
      error("weder CERT_DIR noch POL gesetzt"));
  basepol = subst(eval(polenv), variable(eval(polenv)), y);
  p = eval(getenv("P"));
  print("ohne Zertifikat, freier Lauf"));
}
print("p = ", p, "   Grundkoerper ", basepol);
gettime();
K = bnfinit(basepol, 1);
print("bnfinit(K)  ", gettime(), " ms   cyc = ", K.cyc, "   h = ", K.no);

\\ --- die p-relevanten Erzeuger und die Torsionsbasis
tor = [];
{for (j = 1, #K.cyc, if (K.cyc[j] % p == 0, tor = concat(tor, [j])));}
rk = #tor;
print("p-Rang = ", rk, "   relevante Erzeuger ", tor);

\\ J_j im Klassenkoerper e_j = (cyc_j/p) g_j
JJ = vector(rk);
{
for (n = 1, rk,
  e = vector(#K.cyc);
  e[tor[n]] = K.cyc[tor[n]]/p;
  JJ[n] = idealred(K, idealfactorback(K, K.gen, e)));
}
print("Normen der J_j: ", vector(rk, n, idealnorm(K, JJ[n])));
print("");

\\ --- Klassenkoerper vom Index p
bnr = bnrinit(K, 1);
subs = subgrouplist(bnr, [p]);
print("Untergruppen vom Index ", p, ": ", #subs);

\\ Frobenius-Exponent einer Primstelle von K in Gal(L|K), relativ zu s0;
\\ benutzt die im Schleifenrumpf gesetzten Globalen nfL, Lrel, ag, s0, zk.
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

\\ Koordinaten einer Idealklasse in Cl(K)/p auf den relevanten Erzeugern
ccoords(id) = {vector(rk, n, lift(Mod(bnfisprincipal(K, id, 0)[tor[n]], p)))~;}

\\ --- Hauptschleife ueber die Klassenkoerper
Dfound = List(); Xfound = List();
{
for (s = 1, #subs,
  gettime();
  relpol = bnrclassfield(bnr, subs[s])[1];
  Lrel = rnfinit(K.nf, relpol);
  Labs = Lrel.polabs;
  L = bnfinit(Labs, 1);
  nfL = L.nf;
  tbnf = gettime();

  \\ Gal(L|K): die Automorphismen, die den Grundkoerper festhalten
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

  \\ Artin-Charakter von s0: Frobenius der Klassengruppenerzeuger
  zk = vector(#nfL.zk, j, nfalgtobasis(nfL, nfL.zk[j]));
  canon = vector(rk);
  for (n = 1, rk,
    fa = idealfactor(K, K.gen[tor[n]]); ex = 0;
    for (j = 1, #fa[,1], ex = ex + frobexp(fa[j,1]) * fa[j,2]);
    canon[n] = lift(Mod(ex, p)));

  \\ Matrix von (1-sigma)^2 auf Cl(L), fuer sigma = s0
  cycL = L.cyc;
  S = matconcat(vector(#cycL, i,
        bnfisprincipal(L, nfgaloisapply(nfL, s0, L.gen[i]), 0)));
  M1 = matid(#cycL) - S;
  M2 = M1 * M1;

  \\ die Spalten der Matrix: fuer jedes J_n den Zeugen loesen
  col = vector(rk);
  ok = 1;
  for (n = 1, rk,
    iJ = rnfidealup(Lrel, JJ[n], 1);
    target = -bnfisprincipal(L, iJ, 0);
    X = matsolvemod(M2, cycL~, target);
    \\ matsolvemod meldet Unloesbarkeit mit der ZAHL 0; ist die Loesung
    \\ selbst der Nullvektor, kommt eine Spalte zurueck -- und in gp ist
    \\ [0,..,0]~ == 0 wahr.  Deshalb am Typ unterscheiden, nicht am Wert.
    if (type(X) == "t_INT", ok = 0; break());
    \\ Reduzierte Potenzprodukte statt idealfactorback: downstream wird
    \\ nur die KLASSE von Ip (via ccoords der Norm) benutzt, also darf
    \\ jeder Faktor LLL-reduziert bleiben; Exponenten bis h bleiben so
    \\ speicherflach (binaere Potenzierung mit Reduktion je Schritt).
    Ip = idealhnf(L, 1);
    for (ii = 1, #cycL,
      if (X[ii] != 0,
        Ip = idealmul(L, Ip, idealpow(L, L.gen[ii], X[ii], 1), 1)));
    Ip = idealhnf(L, Ip);
    NI = rnfidealnormrel(Lrel, rnfidealabstorel(Lrel, Ip));
    corr = if (p == 3, idealmul(K, NI, JJ[n]), NI);
    col[n] = ccoords(corr));

  if (!ok,
    print("  s=", s, "  canon=", canon, "  KEINE LOESUNG fuer (1-sigma)^2 X = -[iJ]"),
    \\ sigma = s0^k realisiert den Charakter canon * k^-1
    for (k = 1, p-1,
      chi = vector(rk, i, lift(Mod(canon[i], p) * Mod(k, p)^-1));
      listput(Xfound, chi);
      listput(Dfound, [s, k, matconcat(col)]));
    printf("  s=%2d  canon=%s  bnfinit(L) %5d ms  gesamt %5d ms\n",
           s, canon, tbnf, gettime() + tbnf));
);
}
print("");
print("berechnete Charaktere: ", #Xfound);
print("");

\\ --- Vergleich mit dem gespeicherten Zertifikat
{
if (!haveCert,
  \\ Die d(d+1)/2 Matrizen fuer die Standardcharaktere herausschreiben:
  \\ zuerst e_1..e_d, dann e_i+e_j fuer i<j.
  need = List();
  for (i = 1, rk, v = vector(rk); v[i] = 1; listput(need, Vec(v)));
  for (i = 1, rk-1, for (j = i+1, rk,
    v = vector(rk); v[i] = 1; v[j] = 1; listput(need, Vec(v))));
  out = "D-matrizen.txt";
  system(concat(["rm -f ", out]));
  write(out, "\\\\ D_x fuer die Standardcharaktere, p = ", p, ", Rang ", rk);
  for (t = 1, #need,
    pos = 0;
    for (n = 1, #Xfound, if (Vec(Xfound[n]) == need[t] && !pos, pos = n));
    if (!pos,
      print("Charakter ", need[t], " nicht erreicht"),
      M = lift(Mod(Dfound[pos][3], p));
      print(need[t], "  ->  ", M);
      write(out, "[", need[t], ", ", M, "]")));
  print("");
  print("geschrieben: ", out);
  quit);
}
labs3 = ["x1","x2","x3","x1+x2+x3","x1+x2","x1+x3"];
labs5 = ["x1","x2","x3","x2+x3","x1+x2","x1+x3"];
let5  = ["a","b","c","b+c","a+b","a+c"];
chis3 = [[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,1,0],[1,0,1]];
chis5 = [[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1]];
labs = if (p == 3, labs3, labs5);
chis = if (p == 3, chis3, chis5);

\\ die sechs Matrizen aus dem Zertifikat zusammensetzen
Dcert = vector(6, i, matrix(rk, rk));
{
for (n = 1, #cert[7],
  e = cert[7][n];
  lab = e[1]; col = e[2];
  idx = 0;
  for (i = 1, 6,
    if (lab == labs[i] || (p != 3 && lab == let5[i]), idx = i));
  if (idx, Dcert[idx][,col] = Col(e[13])));
}

print("=== Vergleich Orakel gegen Zertifikat ===");
{
agree = 0; miss = 0;
for (i = 1, 6,
  pos = 0;
  for (n = 1, #Xfound, if (Vec(Xfound[n]) == chis[i] && !pos, pos = n));
  if (!pos,
    print(labs[i], "  -- vom Orakel nicht erreicht"); miss = miss + 1
  ,
    Do = lift(Mod(Dfound[pos][3], p));
    Dc = lift(Mod(Dcert[i], p));
    st = if (Do == Dc, "GLEICH",
             if (Do == lift(Mod(-Dc, p)), "GLEICH BIS AUFS VORZEICHEN", "VERSCHIEDEN"));
    if (Do == Dc, agree = agree + 1);
    print(labs[i], "   ", st);
    print("   Orakel:     ", Do);
    print("   Zertifikat: ", Dc)));
print("");
print("uebereinstimmend: ", agree, " von 6,   nicht erreicht: ", miss);
}
quit
