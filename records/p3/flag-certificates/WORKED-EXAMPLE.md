# A worked example: discriminant -278439307

Certificate `278439307.json`; coefficient field $\mathbf F_3$.

## 1. Relations and change of variables

The relation tensor is the entry `tensor_3_by_27` of the census
record for this discriminant in `records/p3/results/slice2/verification-starthinker.json`.  The certificate prescribes
`basis = [[1, 0, 0], [1, 1, 0], [0, 0, 1]]` (row convention
$\varphi(X_i) = \sum_a \text{basis}[i][a]\,X_a$).  After the
substitution and Gaussian elimination (words of equal degree ordered
by $X_1 > X_2 > X_3$), the three relations become

$r_{1} = X_1 X_1 X_3 + X_1 X_2 X_3 + X_1 X_3 X_1 + 2\, X_2 X_1 X_3 + 2\, X_2 X_3 X_3 + X_3 X_1 X_1 + 2\, X_3 X_1 X_2 + X_3 X_2 X_1 + 2\, X_3 X_2 X_3 + 2\, X_3 X_3 X_2 + 2\, X_3 X_3 X_3$   (head word $X_1 X_1 X_3$)

$r_{2} = X_1 X_2 X_2 + 2\, X_1 X_3 X_2 + X_2 X_1 X_2 + X_2 X_1 X_3 + X_2 X_2 X_1 + X_2 X_2 X_3 + 2\, X_2 X_3 X_1 + X_2 X_3 X_2 + 2\, X_2 X_3 X_3 + X_3 X_1 X_2 + X_3 X_2 X_2 + 2\, X_3 X_2 X_3 + 2\, X_3 X_3 X_2$   (head word $X_1 X_2 X_2$)

$r_{3} = X_1 X_3 X_3 + 2\, X_2 X_2 X_3 + 2\, X_2 X_3 X_2 + X_2 X_3 X_3 + X_3 X_1 X_3 + 2\, X_3 X_2 X_2 + X_3 X_2 X_3 + X_3 X_3 X_1 + X_3 X_3 X_2$   (head word $X_1 X_3 X_3$)

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

$r_{4} = X_1 X_2 X_3 X_2 + X_1 X_3 X_2 X_2 + 2\, X_2 X_1 X_2 X_3 + 2\, X_2 X_2 X_1 X_3 + X_2 X_2 X_2 X_3 + 2\, X_2 X_2 X_3 X_1 + X_2 X_2 X_3 X_3 + 2\, X_2 X_3 X_2 X_1 + X_2 X_3 X_2 X_3 + X_3 X_1 X_3 X_2 + 2\, X_3 X_2 X_1 X_3 + 2\, X_3 X_2 X_2 X_1 + 2\, X_3 X_2 X_2 X_2 + 2\, X_3 X_2 X_2 X_3 + X_3 X_2 X_3 X_1 + X_3 X_2 X_3 X_2 + X_3 X_2 X_3 X_3 + 2\, X_3 X_3 X_1 X_2 + X_3 X_3 X_2 X_2 + X_3 X_3 X_2 X_3 + X_3 X_3 X_3 X_2$

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

Hence the relations are strongly free over the stated field;
the 3-class tower group of the field of discriminant $-278439307$ is mild
of cohomological dimension 2.
