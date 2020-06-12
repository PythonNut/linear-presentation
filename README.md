# A (sort-of) Linear-Time Algorithm for Drawing Linear Knot Diagrams
## Motivation
When doing work in [Knot
Theory](https://en.wikipedia.org/wiki/Knot_theory), it is often
helpful to use _combinatorial encodings_ such as the [_Gauss
code_](http://katlas.math.toronto.edu/wiki/Gauss_Codes) to represent
our knots (another common encoding is the [_DT
code_](http://katlas.org/wiki/DT_%28Dowker-Thistlethwaite%29_Codes),
but we won't talk about that today). Such encodings offer computers an
elegant way of representing / manipulating knots --- indeed, many
invariants can be computed directly from combinatorial
representations.

Basically, the Gauss code is a really elegant tool for applying
algebraic/algorithmic techniques to the study of knots. However, it
can be challenging to build geometric intuition for how manipulations
of the Gauss code translate to manipulations of a knot and vice versa.
The problem is not the translation process is inherently challenging,
but rather that it is tedious and time-consuming.

`linear-presentation` is designed to streamline this process. Given a
Gauss code, it creates a diagram for the associated knot such that
  (1) all of the crossings are collinear,
  (2) the over/understrands of each crossing are perpendicular at the
    crossing point, and
  (2) the first time a crossing is encountered, the strand is oriented
    pointing to the right.
Since Gauss codes are presented as linear strings, there is a sense in
which this feels like the "natural" diagram to associate to the code.
And in the other direction: In our personal experience, this
presentation makes writing down a Gauss code associated to a diagram
much faster than it would otherwise be.

## A Brief Description of the Algorithm
