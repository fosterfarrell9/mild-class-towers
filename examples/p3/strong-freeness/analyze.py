"""Evaluation-point constants, vendored for the public tree.

The values are the ones fixed by the reconstruction pipeline; see the
paper's reconstruction section.
"""

BASIS_POINTS = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
MINIMAL_POINTS = BASIS_POINTS + [(1, 1, 1)]
SIX_POINTS = MINIMAL_POINTS + [(1, 1, 0), (1, 0, 1)]
