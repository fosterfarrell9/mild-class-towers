\\ Deterministic format-2 certificate builder for p=3.
\\
\\ The relative BNF is a search oracle only.  Every emitted entry has already
\\ passed exact field, Artin, AC1, AC2, J-term, and class-coordinate checks.
\\ The accepted output is [N_{L/K}(I') * J] in Cl(K)/3.
\\
\\ Run from the repository root:
\\   P3_DISC=-3640387 P3_EXPECTED_CYC='[18,3,3]' \\
\\   P3_RESULT_DIR=/tmp/p3-build P3_CERT_PATH=/tmp/certificate.gp \\
\\   P3_EXPECTED_TENSOR='[[...],[...],[...]]' \\
\\     gp -qf experiments/p3-cert/build_certificate.gp

default(parisizemax, if(type(getenv("P3_PARISIZEMAX")) == "t_STR" && getenv("P3_PARISIZEMAX") != "", eval(getenv("P3_PARISIZEMAX")), 4*10^9));

P = 3;
SURVEY_FIELDS = [[-3321607,[63,3,3]],[-3640387,[18,3,3]],[-4019207,[207,3,3]],[-4447704,[24,6,6]],[-4472360,[30,6,6]],[-4818916,[48,3,3]],[-4897363,[33,3,3]],[-5048347,[18,6,3]],[-5067967,[69,3,3]],[-5153431,[216,3,3]],[-5288968,[72,3,3]],[-5769988,[12,6,6]]];
P3_DISC_ENV = getenv("P3_DISC");
EXPECTED_DISC = if(P3_DISC_ENV == "", -3640387, eval(P3_DISC_ENV));
P3_CYC_ENV = getenv("P3_EXPECTED_CYC");
if (P3_CYC_ENV == "", FIELD_RECORD = select(row -> row[1] == EXPECTED_DISC, SURVEY_FIELDS); if (#FIELD_RECORD != 1, error("discriminant is not in the fixed survey list and P3_EXPECTED_CYC is unset")); EXPECTED_CYC = FIELD_RECORD[1][2], EXPECTED_CYC = eval(P3_CYC_ENV));
BASE_POL = if (EXPECTED_DISC % 4 == 1, y^2-y+(1-EXPECTED_DISC)/4, y^2-EXPECTED_DISC/4);
P3_RESULT_ENV = getenv("P3_RESULT_DIR");
if (P3_RESULT_ENV == "", error("P3_RESULT_DIR is required"));
RESULT_DIR = P3_RESULT_ENV;
MATRIX_TSV = concat(RESULT_DIR, "/matrices.tsv");
SUMMARY_GP = concat(RESULT_DIR, "/arithmetic-summary.gp");
CERT_PATH = getenv("P3_CERT_PATH");
if (CERT_PATH == "", error("P3_CERT_PATH is required"));
P3_TENSOR_ENV = getenv("P3_EXPECTED_TENSOR");
if (P3_TENSOR_ENV == "", error("P3_EXPECTED_TENSOR is required"));
EXPECTED_TENSOR = eval(P3_TENSOR_ENV);
if (#EXPECTED_TENSOR != 3 || sum(i = 1, 3, #EXPECTED_TENSOR[i] != 27), error("expected tensor must have shape 3 by 27"));
PARI_VERSION_VECTOR = version();
PARI_VERSION_CODE_VALUE = PARI_VERSION_VECTOR[1]*65536 + PARI_VERSION_VECTOR[2]*256 + PARI_VERSION_VECTOR[3];
CERT_HANDLE = 0;
CERT_FIRST_ENTRY = 1;

CHARACTERS = [["x1",[1,0,0]],["x2",[0,1,0]],["x3",[0,0,1]],["x1+x2+x3",[1,1,1]],["x1+x2",[1,1,0]],["x1+x3",[1,0,1]]];

modp(a) = ((a % P) + P) % P;

row_mod(M, cyc) =
{
  my(R = M);
  for (i = 1, matsize(R)[1],
    for (j = 1, matsize(R)[2], R[i,j] = R[i,j] % cyc[i]));
  R;
}

vector_mod(v, cyc) = vector(#v, i, v[i] % cyc[i]);

p_coordinates(Kb, idxK, ideal) =
{
  my(c = bnfisprincipal(Kb, ideal, 0));
  vector(#idxK, i, modp(c[idxK[i]]));
}

\\ Reduced binary powering: keeps every intermediate ideal reduced, so the
\\ construction stays within the default stack on fields with large class
\\ groups.  The result is a possibly different representative of the same
\\ class; every subsequent audit works with the representative chosen here.
idealpow_reduced(Lb, I, e) =
{
  my(ans = idealhnf(Lb, 1), base = idealred(Lb, idealhnf(Lb, I)));
  while (e,
    if (e % 2, ans = idealred(Lb, idealmul(Lb, ans, base)));
    base = idealred(Lb, idealmul(Lb, base, base));
    e \= 2);
  ans;
}

ideal_from_coords(Lb, coordinates) =
{
  my(ans = idealhnf(Lb, 1));
  for (i = 1, #coordinates,
    if (coordinates[i],
      ans = idealred(Lb,
        idealmul(Lb, ans, idealpow_reduced(Lb, Lb.gen[i], coordinates[i])))));
  idealred(Lb, ans);
}

one_minus_ideal(Lb, aut, ideal) =
{
  idealhnf(Lb, idealmul(Lb, ideal,
    idealinv(Lb, nfgaloisapply(Lb, aut, ideal))));
}

aut_power(Lb, abs_pol, aut, exponent) =
{
  if (exponent == 0, x,
    if (exponent == 1, lift(aut), lift(nfgaloisapply(Lb, aut, aut))));
}

select_relative_generator(Lb, abs_pol, y_in_L) =
{
  my(auts = nfgaloisconj(Lb), alpha = Mod(x, abs_pol), choices = List());
  for (i = 1, #auts,
    if (nfgaloisapply(Lb, auts[i], y_in_L) == y_in_L && auts[i] != alpha,
      listput(choices, auts[i])));
  if (#choices != 2, error("expected two nontrivial elements of Gal(L/K)"));
  if (nfgaloisapply(Lb, choices[1],
        nfgaloisapply(Lb, choices[1],
          nfgaloisapply(Lb, choices[1], alpha))) != alpha,
    error("selected automorphism does not have order three"));
  lift(choices[1]);
}

frob_prime(Kb, Lb, rnf, abs_pol, aut, prime_K) =
{
  my(up = rnfidealup(rnf, prime_K, 1), fac_L, rows, prime_L,
     modpr, norm_prime, answer = -1, ok, left, right, candidate);
  fac_L = idealfactor(Lb, up);
  rows = matsize(fac_L)[1];
  if (rows == 3, return(0));
  if (rows != 1, error("unexpected splitting count in cyclic cubic: ", rows));
  if (fac_L[1,2] != 1, error("ramification detected in declared unramified field"));
  prime_L = fac_L[1,1];
  modpr = nfmodprinit(Lb, prime_L);
  norm_prime = idealnorm(Kb, prime_K);
  for (a = 0, P-1,
    ok = 1;
    candidate = aut_power(Lb, abs_pol, aut, a);
    for (z = 1, #Lb.zk,
      left = nfmodpr(Lb, nfgaloisapply(Lb, candidate, Lb.zk[z]), modpr);
      right = nfmodpr(Lb, Lb.zk[z], modpr)^norm_prime;
      if (left != right, ok = 0; break));
    if (ok, answer = a; break));
  if (answer < 0, error("could not identify exact Frobenius exponent"));
  answer;
}

artin_ideal(Kb, Lb, rnf, abs_pol, aut, ideal_K) =
{
  my(fac_K = idealfactor(Kb, ideal_K), answer = 0);
  for (j = 1, matsize(fac_K)[1],
    answer = modp(answer + fac_K[j,2] *
      frob_prime(Kb, Lb, rnf, abs_pol, aut, fac_K[j,1])));
  answer;
}

compact_norm(Kb, rnf, compact) =
{
  my(n = matsize(compact)[1]);
  matrix(n, 2, i, j,
    if (j == 1,
      nfalgtobasis(Kb,
        rnfeltnorm(rnf, rnfeltabstorel(rnf, compact[i,1]))),
      compact[i,2]));
}

compact_inverse(compact) =
{
  matrix(matsize(compact)[1], 2, i, j,
    if (j == 1, compact[i,j], -compact[i,j]));
}

compact_append(compact, factor, exponent) =
{
  my(n = matsize(compact)[1]);
  matrix(n+1, 2, i, j,
    if (i <= n, compact[i,j], if (j == 1, factor, exponent)));
}

\\ Choose a deterministic odd residue characteristic at which all compact
\\ factors are units and verify that their quotient is +1 rather than -1.
compact_modular_witness(Kb, compact) =
{
  my(ell = 3, primes, prime, suitable, modpr, residue, plus_one, minus_one);
  while (1,
    if (isprime(ell),
      primes = idealprimedec(Kb, ell);
      for (q = 1, #primes,
        prime = primes[q]; suitable = 1;
        for (i = 1, matsize(compact)[1],
          if (nfeltval(Kb, compact[i,1], prime) != 0,
            suitable = 0; break));
        if (suitable,
          modpr = nfmodprinit(Kb, prime);
          residue = nfmodpr(Kb, 1, modpr);
          for (i = 1, matsize(compact)[1],
            residue *= nfmodpr(Kb, compact[i,1], modpr)^compact[i,2]);
          plus_one = nfmodpr(Kb, 1, modpr);
          minus_one = nfmodpr(Kb, -1, modpr);
          if (plus_one == minus_one,
            error("odd residue field does not distinguish signs"));
          if (residue != plus_one,
            error("exact AC2 modular sign is not +1"));
          return([ell, prime]))));
    ell += 2);
}

\\ Exact principal ideal of a compact element, without expanding the element.
\\ Factor the principal ideal of each compact factor, combine equal prime
\\ ideals and their signed exponents first, then materialize the cancelled
\\ factorization.  This is the GP-level analogue of PARI's internal
\\ famat_idealfactor route, implemented only with public GP functions.
compact_principal_ideal(nf, compact) =
{
  my(primes = List(), exponents = List(), fac, prime, exponent, found,
     n, factorization);
  for (i = 1, matsize(compact)[1],
    if (compact[i,2],
      fac = idealfactor(nf, idealhnf(nf, compact[i,1]));
      for (j = 1, matsize(fac)[1],
        prime = fac[j,1];
        exponent = compact[i,2] * fac[j,2];
        found = 0;
        for (k = 1, #primes,
          if (primes[k] == prime,
            exponents[k] += exponent; found = 1; break));
        if (!found,
          listput(primes, prime); listput(exponents, exponent)))));
  n = #primes;
  factorization = matrix(n, 2, i, j,
    if (j == 1, primes[i], exponents[i]));
  idealhnf(nf, idealfactorback(nf, factorization));
}

\\ Return [t_AC, sign correction, AC2 unit coordinates].  If t_code is the
\\ compact generator of i(J)*(1-sigma)^2 I', then t_AC=t_code^(-1), up to
\\ the unique sign correction.  Since [L:K]=3, N(-1)=-1, so both base units
\\ are exact norms and no search through extension units is needed.
normalize_t_ac(Kb, rnf, t_code, a_prime) =
{
  my(norm_code = compact_norm(Kb, rnf, t_code),
     unit = bnfisunit(Kb, compact_append(norm_code, a_prime, 1)),
     t_ac, correction = 1, norm_ac, quotient);
  if (#unit == 0, error("N(t_code)*a' is not a base unit"));
  if (#unit != 1, error("unexpected unit rank for imaginary quadratic base"));
  if (unit[1] != 0,
    correction = -1);
  t_ac = compact_inverse(t_code);
  if (correction == -1, t_ac = compact_append(t_ac, -1, 1));
  norm_ac = compact_norm(Kb, rnf, t_ac);
  quotient = bnfisunit(Kb, compact_append(norm_ac, a_prime, -1));
  if (#quotient == 0 || sum(i = 1, #quotient, quotient[i] != 0),
    error("exact norm identity N(t_AC)=a' failed"));
  [t_ac, correction, quotient];
}

solve_ac_value(Kb, idxK, Lb, Lbnr, rnf, aut, J, a_prime) =
{
  my(started = getwalltime(), cyc_L = Lb.cyc, sigma_matrix,
     one_minus, operator, iJ, rhs, solution_data, coordinates,
     I_prime, operated, principal, principal_data, t_code, normalized,
     t_ac, t_ac_ideal, ac1_hnf, base_pair_hnf, norm_I, norm_only,
     norm_plus_J, D_value, J_value, norm_t, quotient_compact, modular);

  sigma_matrix = bnrgaloismatrix(Lbnr, aut);
  one_minus = row_mod(matid(#cyc_L) - sigma_matrix, cyc_L);
  operator = row_mod(one_minus * one_minus, cyc_L);
  iJ = rnfidealup(rnf, J, 1);
  rhs = -bnfisprincipal(Lb, iJ, 0);
  solution_data = matsolvemod(operator, Col(cyc_L), rhs, 1);
  if (solution_data == 0, error("AC class equation has no solution"));
  coordinates = vector_mod(Vec(solution_data[1]), cyc_L);
  I_prime = ideal_from_coords(Lb, coordinates);

  operated = one_minus_ideal(Lb, aut, one_minus_ideal(Lb, aut, I_prime));
  principal = idealhnf(Lb, idealmul(Lb, operated, iJ));
  principal_data = bnfisprincipal(Lb, principal, 5);
  if (sum(i = 1, #principal_data[1], principal_data[1][i] != 0),
    error("AC1 ideal is not principal"));
  t_code = principal_data[2];
  normalized = normalize_t_ac(Kb, rnf, t_code, a_prime);
  t_ac = normalized[1];

  \\ Independently reconstruct the principal ideal of the stored compact
  \\ t_AC after cancelling its prime-ideal factorization.  This is an exact
  \\ HNF computation, not a class-group-coordinate check.
  t_ac_ideal = compact_principal_ideal(Lb, t_ac);
  if (t_ac_ideal != idealhnf(Lb, idealinv(Lb, principal)),
    error("compact t_AC principal-ideal audit failed"));
  ac1_hnf = idealhnf(Lb,
    idealmul(Lb, idealmul(Lb, operated, t_ac_ideal), iJ));
  if (ac1_hnf != idealhnf(Lb, 1),
    error("exact HNF audit of AC1 failed"));

  base_pair_hnf = idealhnf(Kb,
    idealmul(Kb, idealhnf(Kb, a_prime), idealpow(Kb, J, P)));
  if (base_pair_hnf != idealhnf(Kb, 1),
    error("div(a')+3J identity failed"));

  norm_I = rnfidealnormrel(rnf, rnfidealabstorel(rnf, I_prime));
  norm_plus_J = idealhnf(Kb, idealmul(Kb, norm_I, J));
  norm_only = p_coordinates(Kb, idxK, norm_I);
  D_value = p_coordinates(Kb, idxK, norm_plus_J);
  J_value = p_coordinates(Kb, idxK, J);
  if (vector(3, i, modp(D_value[i] - norm_only[i])) != J_value,
    error("J-term coordinate audit failed"));

  norm_t = compact_norm(Kb, rnf, t_ac);
  quotient_compact = compact_append(norm_t, a_prime, -1);
  modular = compact_modular_witness(Kb, quotient_compact);

  [I_prime, coordinates, t_ac, norm_I, norm_plus_J,
   norm_only, J_value, D_value, 1, 1, 1,
   getwalltime() - started, modular[1], modular[2]];
}

subgroup_for_character(cyc_K, idxK, point) =
{
  my(chi = vector(#cyc_K, i, 0));
  for (j = 1, #idxK,
    chi[idxK[j]] = ((cyc_K[idxK[j]] / P) * point[j]) % cyc_K[idxK[j]]);
  charker(cyc_K, chi);
}

write_matrix_rows(handle, label, point, doubled, entries) =
{
  for (ell = 1, #entries,
    my(e = entries[ell], D = e[8], N = e[6], Jc = e[7]);
    filewrite(handle, Str(label, "\t", point[1], "\t", point[2], "\t",
      point[3], "\t", doubled, "\t", ell, "\t",
      D[1], "\t", D[2], "\t", D[3], "\t",
      N[1], "\t", N[2], "\t", N[3], "\t",
      Jc[1], "\t", Jc[2], "\t", Jc[3], "\t", e[12])));
}

write_certificate_entries(Lb, rel_pol, abs_pol, aut_x, label, point, J_basis, a_primes, entries) =
{
  for (ell = 1, #entries,
    my(e = entries[ell], entry =
      [label, ell, Col(point), rel_pol, abs_pol, nfalgtobasis(Lb, aut_x),
       a_primes[ell], J_basis[ell], e[1], e[3], e[13], e[14],
       Col(e[8]), Lb.zk]);
    if (!CERT_FIRST_ENTRY, filewrite1(CERT_HANDLE, ",\n"));
    CERT_FIRST_ENTRY = 0;
    filewrite1(CERT_HANDLE, Str(entry)));
}

analyze_character(Kb, bnr, idxK, kappa, J_basis, a_primes, label, point, matrix_handle) =
{
  my(started = getwalltime(), H, rel_pol, equation, abs_pol, y_in_L,
     Lb, rnf, Lbnr, aut0, artin0, pivot = 0, lambda, aut_x,
     artin_x, entries, doubled_entries, double_point, aut_2x,
     outpath, out, relative_degree, disc_ok);

  print("CHARACTER ", label, " point=", point);
  H = subgroup_for_character(Kb.cyc, idxK, point);
  if (abs(matdet(H)) != P, error("character subgroup does not have index 3"));
  rel_pol = bnrclassfield(bnr, H, 1);
  equation = rnfequation(BASE_POL, rel_pol, 1);
  abs_pol = equation[1];
  y_in_L = Mod(equation[2], abs_pol);
  Lb = bnfinit(abs_pol, 1);       \\ GRH-dependent, deliberately uncertified
  rnf = rnfinit(Kb, rel_pol);
  nfinit(rnf);
  Lbnr = bnrinit(Lb, 1, 1);
  relative_degree = poldegree(abs_pol) / poldegree(BASE_POL);
  disc_ok = nfdisc(abs_pol) == EXPECTED_DISC^P;
  if (relative_degree != P || !disc_ok,
    error("class field degree/discriminant unramifiedness audit failed"));

  aut0 = select_relative_generator(Lb, abs_pol, y_in_L);
  artin0 = vector(3, i, artin_ideal(Kb, Lb, rnf, abs_pol, aut0, kappa[i]));
  for (i = 1, 3, if (point[i] != 0, pivot = i; break));
  lambda = lift(Mod(artin0[pivot], P) / Mod(point[pivot], P));
  if (lambda == 0, error("Artin character is zero"));
  if (vector(3, i, modp(artin0[i] - lambda*point[i])) != [0,0,0],
    error("Artin character is not the prescribed projective line"));
  aut_x = aut_power(Lb, abs_pol, aut0, lambda);
  artin_x = vector(3, i, artin_ideal(Kb, Lb, rnf, abs_pol, aut_x, kappa[i]));
  if (artin_x != point,
    error("normalized automorphism does not realize prescribed character"));

  entries = vector(3, ell,
    solve_ac_value(Kb, idxK, Lb, Lbnr, rnf, aut_x,
                   J_basis[ell], a_primes[ell]));

  double_point = vector(3, i, modp(2*point[i]));
  aut_2x = aut_power(Lb, abs_pol, aut_x, 2);
  if (vector(3, i,
      artin_ideal(Kb, Lb, rnf, abs_pol, aut_2x, kappa[i])) != double_point,
    error("sigma_(2x) orientation audit failed"));
  doubled_entries = vector(3, ell,
    solve_ac_value(Kb, idxK, Lb, Lbnr, rnf, aut_2x,
                   J_basis[ell], a_primes[ell]));
  for (ell = 1, 3,
    if (entries[ell][8] != doubled_entries[ell][8],
      error("arithmetic D_(2x)=D_x audit failed for ", label, ", e_", ell)));

  write_certificate_entries(Lb, rel_pol, abs_pol, aut_x, label, point,
                            J_basis, a_primes, entries);

  write_matrix_rows(matrix_handle, label, point, 0, entries);
  write_matrix_rows(matrix_handle, label, double_point, 1, doubled_entries);

  outpath = concat([RESULT_DIR, "/character-", label, ".gp"]);
  out = fileopen(outpath, "w");
  filewrite(out, Str("character_data = ",
    [2, label, point, H, rel_pol, abs_pol, Lb.no, Lb.cyc,
     nfdisc(abs_pol), aut_x, artin_x, entries,
     double_point, aut_2x, doubled_entries,
     getwalltime() - started], ";"));
  fileclose(out);
  print("  h(L)=", Lb.no, " Cl(L)=", Lb.cyc,
        " Artin=", artin_x, " elapsed_ms=", getwalltime()-started);
  [label, point, Lb.no, Lb.cyc, artin_x, entries,
   double_point, doubled_entries, getwalltime()-started];
}

{
  my(total_started = getwalltime(), Kb, certified, idxK, kappa,
     J_basis, a_primes, torsion_coordinates, bnr, matrix_handle,
     results, summary_handle, pari_version_string);

  system(concat(["mkdir -p '", RESULT_DIR, "'"]));
  matrix_handle = fileopen(MATRIX_TSV, "w");
  filewrite(matrix_handle,
    "label\tx1\tx2\tx3\tdoubled\tinput\td1\td2\td3\t"
    "norm_only1\tnorm_only2\tnorm_only3\tj1\tj2\tj3\telapsed_ms");

  Kb = bnfinit(BASE_POL, 1);
  if (nfdisc(BASE_POL) != EXPECTED_DISC, error("base discriminant mismatch"));
  if (Kb.cyc != EXPECTED_CYC || Kb.no != prod(i = 1, #EXPECTED_CYC, EXPECTED_CYC[i]),
    error("base class group mismatch"));
  certified = bnfcertify(Kb);
  if (certified != 1, error("mandatory bnfcertify(K) failed"));
  idxK = select(i -> Kb.cyc[i] % P == 0, vector(#Kb.cyc, i, i));
  if (#idxK != 3, error("expected 3-class rank three"));

  CERT_HANDLE = fileopen(CERT_PATH, "w");
  filewrite1(CERT_HANDLE, Str(
    "[2,", PARI_VERSION_CODE_VALUE, ",", P, ",", BASE_POL, ",",
    EXPECTED_DISC, ",[", Kb.cyc, ",", Kb.no, ",", Kb.gen, ",",
    Kb.tu[1], ",", Kb.tu[2], ",", Kb.zk, ",", EXPECTED_TENSOR,
    "],[\n"));

  \\ Repository convention: kappa_i is PARI's p-relevant class-group
  \\ generator, used modulo 3.  J_i=kappa_i^(cyc_i/3) is the matching
  \\ basis of Cl(K)[3].
  kappa = vector(3, i, idealhnf(Kb, Kb.gen[idxK[i]]));
  J_basis = vector(3, i,
    idealred(Kb, idealpow(Kb, kappa[i], Kb.cyc[idxK[i]] / P)));
  a_primes = vector(3, i,
    nfeltdiv(Kb, 1,
      bnfisprincipal(Kb, idealpow(Kb, J_basis[i], P), 3)[2]));
  torsion_coordinates = matrix(3, 3, i, j,
    if (i == j, 1, 0));
  for (ell = 1, 3,
    my(class_coordinates = bnfisprincipal(Kb, J_basis[ell], 0),
       expected_coordinates = vector(#Kb.cyc, i,
         if (i == idxK[ell], Kb.cyc[i] / P, 0)));
    if (vector(#Kb.cyc, i,
        (class_coordinates[i] - expected_coordinates[i]) % Kb.cyc[i]) !=
        vector(#Kb.cyc, i, 0),
      error("Cl(K)[3] basis-coordinate audit failed for e_", ell));
    if (idealhnf(Kb,
        idealmul(Kb, idealhnf(Kb, a_primes[ell]),
                       idealpow(Kb, J_basis[ell], P))) != idealhnf(Kb, 1),
      error("base AC pair failed for e_", ell)));

  bnr = bnrinit(Kb, 1, 1);
  results = vector(#CHARACTERS, q,
    analyze_character(Kb, bnr, idxK, kappa, J_basis, a_primes,
      CHARACTERS[q][1], CHARACTERS[q][2], matrix_handle));
  fileclose(matrix_handle);
  filewrite1(CERT_HANDLE, "\n]]\n");
  fileclose(CERT_HANDLE);

  pari_version_string = Str(version());
  summary_handle = fileopen(SUMMARY_GP, "w");
  filewrite(summary_handle, Str("arithmetic_summary = ",
    [2, P, BASE_POL, EXPECTED_DISC, Kb.cyc, Kb.no, certified,
     idxK, kappa, J_basis, a_primes, torsion_coordinates,
     pari_version_string, results, getwalltime()-total_started], ";"));
  fileclose(summary_handle);
  print("ARITHMETIC PILOT COMPLETE elapsed_ms=", getwalltime()-total_started);
}

quit
