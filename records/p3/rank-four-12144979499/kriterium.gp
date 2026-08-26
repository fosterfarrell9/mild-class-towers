\\ Transversalitaetskriterium (p = 3, Rang d >= 3) auf den D-Matrizen des Orakels.
\\ Liest D-matrizen.txt im Arbeitsverzeichnis.  Ablauf:
\\   1. Bockstein-Matrix B[l,i] = Koeffizient von X_i^3 in rho_l = D_i[i,l];
\\      Kegel C_beta = P(ker B).  Leer => Kriterium nicht anwendbar.
\\   2. Punkte des Kegels ueber F_{3^e} durchgehen (e bis 4 je nach Kegeldim.),
\\      Rangbedingung rk D_x = d-2 pruefen.
\\   3. An Rangpunkten Theta_x: (F_q^d)/x -> Hom(ker D_x, x^perp/im D_x)
\\      aufstellen ((d-1) x 2 - Matrix nach Wahl eines Funktionals ell) und
\\      Surjektivitaet (Rang 2) testen.  Surjektiv => transversales Element
\\      => mild (Transversalitaetskriterium, d >= 3).
p = 3;
lines = readstr("D-matrizen.txt");
CH = List(); DD = List();
{
for (i = 1, #lines,
  my(s = lines[i]);
  if (#s > 0 && Vecsmall(s)[1] == 91,          \\ Zeile beginnt mit '['
    my(v = eval(s));
    listput(CH, v[1]); listput(DD, v[2])));
}
d = #CH[1];
printf("d = %d, Matrizen: %d\n", d, #CH);
E(i) = vector(d, j, (j==i));
getD(chi) = {for (n = 1, #CH, if (CH[n] == chi, return(DD[n]))); error("Matrix fehlt: ", chi)}
Di = vector(d, i, getD(E(i)));
\\ Polarisation: Del[i][j] = D_{e_i+e_j} - D_i - D_j fuer i != j, und 2 D_i fuer i = j
Del = vector(d, i, vector(d));
{
for (i = 1, d, Del[i][i] = 2*Di[i]);
for (i = 1, d-1, for (j = i+1, d,
  my(M = getD(E(i)+E(j)) - Di[i] - Di[j]);
  Del[i][j] = M; Del[j][i] = M));
}
B = matrix(d, d, l, i, Di[i][i,l]);
K3 = matker(Mod(B, 3));
c = #K3;
printf("Bockstein-Matrix: Rang %d, Kegel P(ker B) projektive Dimension %d\n", d - c, c - 1);
if (c == 0, print("VERDIKT: KEGEL LEER - Kriterium nicht anwendbar"); quit);
kerB = lift(K3);
emax = if (c == 1, 1, if (c == 2, 4, 3));
{
found = 0;
for (e = 1, emax,
  my(q = 3^e, a = ffgen(3^e, 't), one = ffgen(3^e, 't)^0);
  my(elts = vector(q, n,
    my(dg = digits(n-1, 3)); sum(i = 1, #dg, dg[#dg+1-i]*a^(i-1), 0*a)));
  my(tested = 0, nrk2 = 0, nrklow = 0, ntrans = 0);
  \\ normierte Vertreter von P(F_q^c): erste nichtverschwindende Koordinate = 1
  for (j0 = 1, c,
    forvec (rest = vector(c - j0, i, [1, q]),
      my(t = vector(c, i, 0*a));
      t[j0] = one;
      for (i = 1, c - j0, t[j0 + i] = elts[rest[i]]);
      my(x = (kerB * t~) * one);
      tested++;
      my(Dx = sum(i = 1, d, x[i]^2 * Di[i])
              + sum(i = 1, d-1, sum(j = i+1, d, x[i]*x[j] * Del[i][j])));
      Dx = Dx * one;
      my(rk = d - #matker(Dx));
      if (rk < d-2, nrklow++; next);
      if (rk > d-2, next);
      nrk2++;
      \\ Orthogonalitaet im D_x c x^perp (auf dem Kegel garantiert; Kontrolle)
      if (#select(y -> y != 0, Vec(x~ * Dx)) > 0,
        print("WARNUNG: x~ D_x != 0 bei e=", e, " t=", t); next);
      my(kerDx = matker(Dx));                       \\ d x 2
      my(leftk = matker(mattranspose(Dx)));         \\ d x 2, enthaelt x
      my(ell = 0);
      for (m = 1, #leftk,
        if (matrank(matconcat([x~, leftk[,m]~]~)) == 2 && ell == 0,
          ell = leftk[,m]));
      if (ell == 0, print("WARNUNG: kein ell gefunden bei e=", e); next);
      \\ Delta D(x, e_j) fuer alle j
      my(Mx = vector(d, j, sum(i = 1, d, x[i] * Del[i][j]) * one));
      \\ Kontrolle: Werte auf ker D_x liegen in x^perp
      my(okperp = 1);
      for (j = 1, d, for (b = 1, 2,
        if (x~ * (Mx[j] * kerDx[,b]) != 0, okperp = 0)));
      if (!okperp,
        print("WARNUNG: Theta-Werte nicht in x^perp bei e=", e, " t=", t); next);
      \\ Theta-Matrix ueber dem Komplement von x
      my(j0x = 0);
      for (j = 1, d, if (x[j] != 0 && !j0x, j0x = j));
      my(js = select(j -> j != j0x, vector(d, j, j)));
      my(T = matrix(d-1, 2, aa, bb, ell~ * (-(Mx[js[aa]] * kerDx[,bb]))));
      if (matrank(T) == 2,
        ntrans++;
        if (!found,
          found = 1;
          print("TRANSVERSAL: e = ", e, ", Kegelkoordinaten t = ", t);
          print("  x = ", x~);
          print("  rk D_x = ", rk, ",  Theta-Matrix = ", T)));
    ));
  printf("e=%d: %d Punkte, davon rk=d-2: %d, rk<d-2: %d, transversal: %d\n",
         e, tested, nrk2, nrklow, ntrans));
if (found,
  print("VERDIKT: TRANSVERSALES ELEMENT GEFUNDEN => mild (Kriterium, d >= 3)"),
  print("VERDIKT: kein transversales Element auf dem Kegel (bis e=", emax, ")"));
}
quit
