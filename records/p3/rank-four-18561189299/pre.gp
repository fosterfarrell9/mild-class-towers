\\ Vorpruefung eines Kandidaten: fundamental? Klassengruppe? 3-Struktur? Polynom.
D = eval(getenv("DISC"));
{
if (!isfundamental(D), print("NICHT FUNDAMENTAL: ", D); quit);
}
cl = quadclassunit(D);
print("D = ", D, "   h = ", cl.no, "   Cl = ", cl.cyc);
t3 = select(x -> x % 3 == 0, cl.cyc);
print("3-relevante Faktoren: ", t3, "   3-Rang = ", #t3);
{
if (#t3 != 4, print("WARNUNG: 3-RANG IST NICHT 4"));
if (#select(x -> x % 9 == 0, t3) == 0,
  print("3-Anteil elementarabelsch (3,3,3,3)"),
  print("3-Anteil NICHT elementar: Faktor durch 9 teilbar"));
}
pol = if (Mod(D,4) == Mod(1,4), Str("y^2 - y + ", (1-D)/4), Str("y^2 + ", -D/4));
print("POL = ", pol);
write("pol.txt", pol);
quit
