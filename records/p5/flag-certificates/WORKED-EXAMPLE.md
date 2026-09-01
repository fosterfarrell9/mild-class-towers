# A worked example at p = 5: discriminant -781922404

Certificate `781922404.json`; coefficient field $\mathbf F_{25} = \mathbf F_5[t]/(t^2+3)$.

## 1. Relations and change of variables

The relation tensor is the entry `cubic_relation_matrix` of
`records/p5/D-781922404/result.gp`, the output of the run whose arithmetic
certificate is `certificates/p5/K-195480601-p5/certificate.gp`; its three
rows are Lie elements (no cube term occurs).  The certificate prescribes
`basis = [[[1, 0], [1, 0], [1, 0]], [[1, 2], [3, 0], [0, 0]], [[1, 1], [0, 0], [0, 0]]]` (row convention
$\varphi(X_i) = \sum_a \text{basis}[i][a]\,X_a$).  After the
substitution and Gaussian elimination (words of equal degree ordered
by $X_1 > X_2 > X_3$), the three relations become

$r_{1} = X_1 X_1 X_3 + (1 + 2t)\, X_1 X_2 X_3 + 3\, X_1 X_3 X_1 + X_1 X_3 X_2 + (3 + 3t)\, X_2 X_1 X_3 + X_2 X_2 X_3 + X_2 X_3 X_1 + 3\, X_2 X_3 X_2 + (2 + t)\, X_2 X_3 X_3 + X_3 X_1 X_1 + (3 + 3t)\, X_3 X_1 X_2 + (1 + 2t)\, X_3 X_2 X_1 + X_3 X_2 X_2 + (1 + 3t)\, X_3 X_2 X_3 + (2 + t)\, X_3 X_3 X_2$   (head word $X_1 X_1 X_3$)

$r_{2} = X_1 X_2 X_2 + (4 + t)\, X_1 X_2 X_3 + (1 + 2t)\, X_1 X_3 X_2 + 3\, X_2 X_1 X_2 + 2t\, X_2 X_1 X_3 + X_2 X_2 X_1 + 4t\, X_2 X_2 X_3 + (1 + 2t)\, X_2 X_3 X_1 + 2t\, X_2 X_3 X_2 + (1 + 3t)\, X_2 X_3 X_3 + 2t\, X_3 X_1 X_2 + (4 + t)\, X_3 X_2 X_1 + 4t\, X_3 X_2 X_2 + (3 + 4t)\, X_3 X_2 X_3 + (1 + 3t)\, X_3 X_3 X_2$   (head word $X_1 X_2 X_2$)

$r_{3} = X_1 X_3 X_3 + (2 + 3t)\, X_2 X_2 X_3 + (1 + 4t)\, X_2 X_3 X_2 + (2 + 3t)\, X_2 X_3 X_3 + 3\, X_3 X_1 X_3 + (2 + 3t)\, X_3 X_2 X_2 + (1 + 4t)\, X_3 X_2 X_3 + X_3 X_3 X_1 + (2 + 3t)\, X_3 X_3 X_2$   (head word $X_1 X_3 X_3$)

## 2. Completion: overlaps and the genuine reduction

An *overlap* (Lemma A.5 of the paper) is a pair of head words $AB$,
$BC$ with nonempty words $A, B, C$; the overlapped word $ABC$ can be
rewritten in two ways, and the difference $gC - Ag'$ of the two
corresponding multiples is reduced with respect to the current basis.
In order of degree (1 overlap in all):

| overlapped word | head words | result |
|---|---|---|
| $X_1 X_1 X_3 X_3$ | $X_1 X_1 X_3$, $X_1 X_3 X_3$ | **new element with head $X_1 X_2 X_3 X_2$** |


The overlap $X_1 X_1 X_3 X_3$ of $X_1 X_1 X_3$ and $X_1 X_3 X_3$ leaves a nonzero
remainder after reduction with respect to the basis found so far:

$r_{4} = X_1 X_2 X_3 X_2 + (4 + t)\, X_1 X_2 X_3 X_3 + 2\, X_1 X_3 X_2 X_2 + 2t\, X_1 X_3 X_2 X_3 + 4\, X_2 X_1 X_2 X_3 + 3\, X_2 X_2 X_1 X_3 + 4t\, X_2 X_2 X_2 X_3 + 3\, X_2 X_2 X_3 X_1 + (1 + t)\, X_2 X_2 X_3 X_2 + (1 + 2t)\, X_2 X_2 X_3 X_3 + (2 + t)\, X_2 X_3 X_1 X_3 + 4\, X_2 X_3 X_2 X_1 + (3 + t)\, X_2 X_3 X_2 X_2 + (3 + t)\, X_2 X_3 X_2 X_3 + (4 + 3t)\, X_2 X_3 X_3 X_1 + (1 + 3t)\, X_2 X_3 X_3 X_2 + (2 + 2t)\, X_2 X_3 X_3 X_3 + (4 + 4t)\, X_3 X_1 X_2 X_3 + (3 + 4t)\, X_3 X_1 X_3 X_2 + 3\, X_3 X_2 X_1 X_3 + 3\, X_3 X_2 X_2 X_1 + (1 + 4t)\, X_3 X_2 X_2 X_2 + 2t\, X_3 X_2 X_2 X_3 + (3 + 4t)\, X_3 X_2 X_3 X_1 + 4t\, X_3 X_2 X_3 X_2 + (2 + 3t)\, X_3 X_2 X_3 X_3 + (4 + 4t)\, X_3 X_3 X_1 X_2 + (3 + 2t)\, X_3 X_3 X_2 X_1 + 3t\, X_3 X_3 X_2 X_2 + 3t\, X_3 X_3 X_2 X_3 + (1 + 2t)\, X_3 X_3 X_3 X_2$

whose head word $X_1 X_2 X_3 X_2$ joins the basis.


Every remaining overlap reduces to $0$, so the completed basis is
finite with head words $X_1 X_1 X_3$, $X_1 X_2 X_2$, $X_1 X_3 X_3$, $X_1 X_2 X_3 X_2$
(Lemma A.5 (ii)).

## 3. Word counts

The set of proper prefixes of the head words (empty word included)
has $s = 6$ elements.  The counts $c_n$ of words with no factor
among the head words satisfy a linear recurrence of order at most
$6$; comparison with the coefficients $b_n$ of $1/(1-3z+3z^3)$ up
to degree $s + 2 = 8$ settles all degrees:

| $n$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| $c_n$ | 1 | 3 | 9 | 24 | 63 | 162 | 414 | 1053 | 2673 |
| $b_n$ | 1 | 3 | 9 | 24 | 63 | 162 | 414 | 1053 | 2673 |

Hence the relations are strongly free over the stated field, and the invariance of the Hilbert series under base change carries the conclusion down to $\mathbf F_5$ (Lemma 3.7);
the 5-class tower group of the field of discriminant $-781922404$ is mild
of cohomological dimension 2 (Proposition 3.8).
