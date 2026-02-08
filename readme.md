# ANIMA-MECHANISMI
An Index of Computer Science for Practitioners

---

> *Comets arrive, bright with borrowed fire,*  
> *promising their blaze will stay.*  
> *Morning steals them anyway.*
> 
> *Beneath the trailing glow,*  
> *all that may be and may never be gathers.*  
> *A river bends unseen,*  
> *tracing arcs like frozen light.*
> 
> *Pearls gleam in your hand,*  
> *but seeds you let fall*  
> *sprout into trees,*  
> *and trees into forests.*
> 
> *After glitter fades,*  
> *a silent thunder*  
> *echoes from the dawn of time.*

---

## Chapter 1: From Particularity to Universality

### Propositional Logic

A proposition is a statement that is either true or false.

**Connectives:**
- `P AND Q` - both P and Q are true
- `P OR Q` - at least one of P or Q is true
- `NOT P` - P is false
- `P IMPLIES Q` - if P is true, then Q is true
- `P IFF Q` - P and Q have the same truth value

**De Morgan's Laws:**
- `NOT (P AND Q) IFF (NOT P) OR (NOT Q)`
- `NOT (P OR Q) IFF (NOT P) AND (NOT Q)`

### First-Order Logic

Extends propositional logic with quantifiers over objects.

**Quantifiers:**
- `FOR ALL x, P(x)` - P holds for every x
- `EXISTS x, P(x)` - P holds for some x

**Variables:**
- Bound variable: quantified by FOR ALL or EXISTS
- Free variable: not quantified

### Proof Techniques

**Direct proof:**
- To prove P IMPLIES Q, assume P is true and derive Q.

**Proof by contradiction:**
- To prove P, assume NOT P and derive a contradiction.

**Proof by contrapositive:**
- To prove P IMPLIES Q, prove (NOT Q) IMPLIES (NOT P).

**Proof by cases:**
- To prove P, split into exhaustive cases and prove each.

**Counterexample:**
- To disprove FOR ALL x, P(x), find one x where P(x) is false.

---

## Chapter 2: From Collection to Connection

### Set

A set is a collection of distinct objects.

**Notation:**
- `{1, 2, 3}` - explicit set
- `x IN S` - x is an element of S
- `A SUBSET B` - every element of A is in B
- `EMPTY SET` - the set with no elements

**Operations:**
- `A UNION B` - elements in A or B
- `A INTERSECT B` - elements in both A and B
- `A MINUS B` - elements in A but not in B
- `POWER SET of A` - set of all subsets of A

### Cartesian Product

The Cartesian product of sets A and B is the set of all ordered pairs.

**Notation:**
- `A CROSS B = {(a, b) | a IN A AND b IN B}`

### Algebraic Data Types: Products and Sums

**Product type `(A CROSS B)`:**
- Combines two types into ordered pairs.

**Sum type `(A + B)`:**
- `({left} CROSS A) UNION ({right} CROSS B)`
- Disjoint union - element from A or B, tagged with origin.

### Relations

A relation R from set A to set B is a subset of `A CROSS B`.

**Notation:**
- `(a, b) IN R`, or: `a R b`

**Properties of a relation R on a set A:**
- Reflexive: `FOR ALL x, x R x`
- Symmetric: `FOR ALL x, y, (x R y) IMPLIES (y R x)`
- Transitive: `FOR ALL x, y, z, (x R y AND y R z) IMPLIES (x R z)`
- Antisymmetric: `FOR ALL x, y, (x R y AND y R x) IMPLIES (x = y)`

### Graphs

Given a relation E on a set V, we can visualize it as a graph.

A graph `G = (V, E)` consists of:
- V - set of vertices
- E - edges on V `(E SUBSET V CROSS V)`
  
If `(u, v) IN E`, we say there is an edge from u to v.

**Types:**
- Directed graph: relation E using ordered pairs `(u, v)`
- Undirected graph: relation E is symmetric

**Terminology:**
- Path: sequence `v1, v2, ..., vn` where `(vi, vi+1) IN E`
- Cycle: path where first and last vertex are the same
- Connected: `FOR ALL u, v IN V, EXISTS path from u to v`
- Tree: connected graph with no cycles

### Equivalence Relations

An equivalence relation is reflexive, symmetric, and transitive.

**Notation:**
- `a EQUALS b` (a and b are equivalent)

In programming, equality (=) is an equivalence relation.

**Equivalence class of x:**
- `[x] = {y | x EQUALS y}`

**Partition:**
- An equivalence relation divides a set into disjoint equivalence classes.

### Partial Orders

A partial order is a relation that is reflexive, antisymmetric, and transitive.

**Notation:**
- `a PRECEDES b` (a precedes b in the order)
- Or: `(S, PRECEDES)` is a partially ordered set (poset)

**Special elements:**
- Minimal: no element precedes it
- Maximal: it precedes no element
- Least: precedes all elements
- Greatest: all elements precede it

### Functions

A function is a special relation where each input maps to exactly one output.

A function `f: A -> B` is a subset of `A CROSS B` where:
- `FOR ALL a IN A, EXISTS unique b IN B, (a, b) IN f`

**Notation:**
- `f(x) = y` means `(x, y) IN f`

**Properties:**
- Injective (one-to-one): different inputs map to different outputs
- Surjective (onto): every element of B is mapped to
- Bijective: both injective and surjective

**Composition:**
- `(g * f)(x) = g(f(x))` - apply f, then apply g

### Indexed Collections

An indexed collection is a function from an index set to a set of values.

Given index set I and set S:
- `f: I -> S`
- For each `i IN I`, `f(i)` is an element of S

### Algebraic Data Types: Function Types

**Function type `(A -> B)`:**
- Type of functions mapping from A to B.
- Example: `Integer -> Boolean` represents functions that take an integer and return a boolean.

**Higher-order functions:**
- Functions that take functions as arguments or return functions.
- Example: `(A -> B) -> (B -> C) -> (A -> C)` is function composition.

---

## Chapter 3: From Enumeration to Generation

### Natural Numbers

The natural numbers `N = {0, 1, 2, 3, ...}` represent counting.

Construction via Peano axioms (details omitted):
- Zero is a natural number
- Every natural number has a successor
- Zero is not the successor of any number
- Different numbers have different successors

### Sequences

A sequence is an indexed collection using natural numbers as the index set.

**Notation:**
- `s: N -> S`
- `s(0), s(1), s(2), ...` or `s_0, s_1, s_2, ...`

**Finite sequences:**
- Functions from `{0, 1, ..., n-1}` to S
- Length n sequence

### Mathematical Induction

Principle of Mathematical Induction:

To prove `FOR ALL n IN N, P(n)`:
1. Base case: Prove `P(0)`
2. Inductive step: Prove `FOR ALL k, P(k) IMPLIES P(k+1)`

**Strong induction:**
- Inductive step uses `P(0), P(1), ..., P(k)` to prove `P(k+1)`

### Inductive Definitions

Objects defined recursively using base cases and inductive rules.

**Example - Factorial:**
- `0! = 1` (base case)
- `n! = n * (n-1)!` (inductive case)

**Example - Sum:**
- `sum(0) = 0` (base case)
- `sum(n) = n + sum(n-1)` (inductive case)

### Algebraic Data Types: Recursive Types

Recursive types are defined in terms of themselves.

**List type:**
- `List(A) = Empty | Cons(A, List(A))`
- A list is either empty or an element paired with another list

**Example lists over N:**
- Empty
- Cons(1, Empty)
- Cons(2, Cons(1, Empty))

**Tree type:**
- `Tree(A) = Leaf | Node(Tree(A), A, Tree(A))`
- A tree is either a leaf or a node with left subtree, value, right subtree

### Recursion

A function is recursive if it calls itself.

**Structure:**
- Base case: direct answer
- Recursive case: reduce to simpler problem

**Example - Fibonacci:**
- `fib(0) = 0`
- `fib(1) = 1`
- `fib(n) = fib(n-1) + fib(n-2)`

**Termination requires:**
- Each recursive call moves toward base case

### Combinatorics

**Permutations of k from n:**
- Recursive: `P(n, k) = n * P(n-1, k-1)`
- Closed form: `P(n, k) = n!/(n-k)!`

**Permutations of n objects:**
- `P(n) = n!`

**Combinations:**
- Recursive: `C(n, k) = C(n-1, k-1) + C(n-1, k)` (Pascal's triangle)
- Closed form: `C(n, k) = n!/(k! * (n-k)!)`

---

## Chapter 4: From Plurality to Uncertainty

### Construction of Number Systems

Number systems are constructed from simpler ones using equivalence relations.

**Integers `Z = {..., -2, -1, 0, 1, 2, ...}`**
- Construction: equivalence classes on `N CROSS N`

**Rationals `Q = {p/q | p, q in Z, q != 0}`**
- Construction: equivalence classes on `Z CROSS (Z - {0})`

**Reals R (all points on the number line)**
- Construction: Dedekind cuts or Cauchy sequences of rationals

### Algebraic Structures

An algebraic structure is a set equipped with one or more operations satisfying specific axioms.

**Formal definition:**
- `(S, *)` where S is a set and `* : S CROSS S -> S` is a binary operation

Different axioms define different structures: monoids, groups, rings, fields.

### Monoids

A monoid is a set with an associative binary operation and an identity element.

**Definition:**
`(M, *, e)` where:
- `* : M CROSS M -> M` (binary operation)
- `FOR ALL a, b, c IN M:`
  - `(a * b) * c = a * (b * c)` (associativity)
  - `e * a = a * e = a` (identity)

Monoids are a fundamental structure in abstract algebra, along with groups, rings, and fields.

**Examples:**
- `(N, +, 0)` - natural numbers under addition
- `(N, *, 1)` - natural numbers under multiplication
- `(Z, +, 0)` - integers under addition
- `(Q, +, 0)` - rationals under addition
- `(Q - {0}, *, 1)` - nonzero rationals under multiplication
- `(R, +, 0)` - reals under addition
- `(sequences over A, concatenation, empty sequence)`
- `(lists over A, concatenation, empty list)`

### Limits and Convergence

**Limit of a sequence:**

A sequence `(s_n)` has limit L if:
- `FOR ALL epsilon > 0, EXISTS N, FOR ALL n >= N, dist(s_n, L) < epsilon`

where `dist(a, b)` is the distance between a and b.

Notation: `LIM (n -> infinity) s_n = L`

The sequence gets arbitrarily close to L for sufficiently large n.

**Convergence:**
- A sequence converges if it has a finite limit.
- A series `SUM(a_n)` converges if its partial sums converge.

### Sample Spaces and Events

A probability space consists of:
- Sample space Omega - set of all possible outcomes
- Event - subset of Omega

**Example:**
- Coin flip: `Omega = {H, T}`
- Die roll: `Omega = {1, 2, 3, 4, 5, 6}`

### Probability Function

A probability function P assigns probabilities to events.

**Axioms:**
- `P(Omega) = 1`
- `P(A) >= 0` for all events A
- `P(A UNION B) = P(A) + P(B)` when A, B are disjoint

### Conditional Probability

**Conditional probability of A given B:**
- `P(A | B) = P(A AND B) / P(B)`

**Independence:**
- Events A and B are independent if `P(A AND B) = P(A) * P(B)`

### Random Variables

A random variable is a function from the sample space to real numbers.

**Formal definition:**
- `X: Omega -> R`

**Probability mass function for discrete random variable X:**
- `p(x) = P(X = x)` for each value x

**Properties:**
- `SUM p(x) = 1` (over all possible values x)
- `p(x) >= 0` for all x

### Expected Value

The expected value (mean) of discrete random variable X:
- `E[X] = SUM(xi * p(xi))`

where the sum is over all values xi with `p(xi) > 0`.

### Variance

The variance of discrete random variable X:
- `Var(X) = E[(X - E[X])^2] = E[X^2] - (E[X])^2`
- `= SUM((xi - E[X])^2 * p(xi))`

**Standard deviation:**
- `SD(X) = sqrt(Var(X))`

### Covariance

Covariance of discrete random variables X and Y:
- `Cov(X, Y) = E[(X - E[X])(Y - E[Y])]`
- `= SUM SUM((xi - E[X])(yj - E[Y]) * p(xi, yj))`

**Correlation (normalized covariance):**
- `Cor(X, Y) = Cov(X, Y) / (SD(X) * SD(Y))`

Ranges from -1 to 1.

### Limit Theorems

**Law of Large Numbers:**

Sample mean converges to population mean as sample size increases.

- `LIM (n -> infinity) (X1 + X2 + ... + Xn) / n = mu`

**Central Limit Theorem:**

Standardized sum of independent random variables converges to normal distribution.

### Entropy

Entropy measures uncertainty in a random variable.

For discrete random variable X with probabilities `p1, p2, ..., pn`:
- `H(X) = - SUM(pi * log_2(pi))`

Measured in bits:
- 1 bit = one binary choice
- `H(X)` = average bits needed to encode X

### Mutual Information

Mutual information measures shared information between random variables.

For random variables X and Y:
- `I(X; Y) = H(X) + H(Y) - H(X, Y)`

where `H(X, Y)` is joint entropy.

Represents how much knowing Y reduces uncertainty about X.

---

## Chapter 5: From Communication to Articulation

### Alphabets

An alphabet is a finite set of symbols.

**Notation:**
- `A = {0, 1}` - binary alphabet
- `A = {a, b, c, ..., z}` - lowercase letters

Alphabets provide the basic symbols for constructing strings.

### Strings

A string over alphabet A is a finite sequence (indexed collection) over A.

**Notation:**
- Empty string: `EPSILON` (length 0)
- Length of string s: `|s|`
- Concatenation: `s1 s2` (join strings together)

**Examples over alphabet `{0, 1}`:**
- EPSILON, 0, 1, 01, 101, 0011

### Kleene Star

The Kleene star `A*` is the set of all strings over alphabet A.

**Formal definition:**
- `A* = A(0) UNION A(1) UNION A(2) UNION A(3) UNION ...`
- where `A(n) = {strings of length n over A}`

- `A(0) = {EPSILON}`
- `A(1) = A`
- `A(2) = {ab | a, b IN A}`
- `A(n) = {s1 s2 | s1 IN A(n-1), s2 IN A}`

**Examples:**
- `{0, 1}* = {EPSILON, 0, 1, 00, 01, 10, 11, 000, ...}`
- `{a}* = {EPSILON, a, aa, aaa, aaaa, ...}`

`A+` denotes all non-empty strings `(A* - {EPSILON})`.

### Regular Languages

A regular language over alphabet A can be built from basic languages using regular operations.

**Basic languages:**
- `EMPTY SET` - {}
- `{EPSILON}` - contains only the empty string
- `{a}` - single symbol, for each a IN A

**Regular operations:**
- Union: `L1 UNION L2`
- Concatenation: `L1 L2 = {s1 s2 | s1 IN L1, s2 IN L2}`
- Kleene star: `L* = L(0) UNION L(1) UNION L(2) UNION ...`
  - where `L(n) = {strings formed by concatenating n strings from L}`

**Closure properties:**

Regular languages are closed under union, intersection, complement, concatenation, and Kleene star.

Regular expressions provide notation for describing regular languages.

### Context-Free Languages

A context-free grammar generates strings through production rules.

**Definition:**

`G = (V, T, P, S)` where:
- V - set of variables (non-terminals)
- T - set of terminals (alphabet)
- P - set of production rules
- S - start variable

Production rules in BNF notation:
- `<variable> ::= <expression>`

**Example - balanced parentheses:**
- `S ::= EPSILON`
- `S ::= (S)`
- `S ::= SS`

Generates: EPSILON, (), (()), ()(), (())(), ...

**Example - arithmetic expressions:**
- `<expr> ::= <number>`
- `<expr> ::= <expr> + <expr>`
- `<expr> ::= <expr> * <expr>`
- `<expr> ::= (<expr>)`

**Example - statements:**
- `<stmt> ::= <var> := <expr>`
- `<stmt> ::= <stmt> ; <stmt>`
- `<stmt> ::= if <expr> then <stmt> else <stmt>`
- `<stmt> ::= while <expr> do <stmt>`

These grammars define the syntax of programs.

---

## Chapter 6: From Interpretation to Validation

### State

A state represents the values of program variables.

Formally, a state is a function `s: Var -> Val`
- where Var is the set of variables
- and Val is the set of values

**Notation:**
- `s(x)` - value of variable x in state s
- `s[x := v]` - state identical to s except x maps to v

**Definition of state update:**
- `s[x := v](y) = v` if `y = x`
- `s[x := v](y) = s(y)` if `y != x`

### Denotational Semantics

Denotational semantics interprets programs as mathematical functions.

Each construct maps to a mathematical object:
- Expressions denote values
- Statements denote state transformations (functions `State -> State`)

**For expressions:**
- `[[n]](s) = n`
- `[[x]](s) = s(x)`
- `[[e1 + e2]](s) = [[e1]](s) + [[e2]](s)`
- `[[e1 * e2]](s) = [[e1]](s) * [[e2]](s)`

**For statements:**
- `[[x := e]](s) = s[x := [[e]](s)]`
  - (state where x maps to value of e)

- `[[c1; c2]](s) = [[c2]]([[c1]](s))`
  - (compose state transformations)

This provides a compositional, mathematical meaning for programs.

### Rewrite Systems

A rewrite system is a pair `(T, R)` where:
- T - set of terms
- R - set of rewrite rules (pairs of terms)

A rewrite rule `(L, R)` is written `L -> R`.

A term t rewrites to t' (written `t => t'`) if t contains L and t' replaces it with R.

**Desirable properties:**
- Confluence: if `t => t1` and `t => t2`, then `EXISTS t3` such that `t1 =>* t3` and `t2 =>* t3`
- Termination: no infinite rewrite sequences exist
- Normal form: a term that cannot be rewritten further

### Lambda Calculus

Lambda calculus is a formal system for expressing computation via functions.

**Syntax:**
- `M ::= x | LAMBDA x. M | M1 M2`

**Beta-reduction:**
- `(LAMBDA x. M) N ==> M[x := N]`
- where `M[x := N]` substitutes N for free occurrences of x in M.

### Operational Semantics: Small-Step

Small-step semantics defines single computation steps.

**Transition relation:**
- `<e, s> -> <e', s'>`
- Expression e in state s reduces to e' in state s'

**Example rules for expressions:**
```
<x, s> -> <s(x), s>

<e1, s> -> <e1', s'>
-------------------------
<e1 + e2, s> -> <e1' + e2, s'>
```

**Example rules for statements:**
```
<x := e, s> -> <x := e', s'>  if <e, s> -> <e', s'>
<x := v, s> -> <skip, s[x := v]>

<c1, s> -> <c1', s'>
-------------------------
<c1; c2, s> -> <c1'; c2, s'>

<skip; c2, s> -> <c2, s>
```

### Operational Semantics: Big-Step

Big-step semantics defines complete evaluation.

**Evaluation relation:**
- `<e, s> => v`
- Expression e in state s evaluates to value v

**Example rules for expressions:**
```
<n, s> => n

<e1, s> => n1    <e2, s> => n2
-----------------------------------
<e1 + e2, s> => n1 + n2
```

**Example rules for statements:**

`<c, s> => s'`

Statement c in state s produces state s'

```
<x := e, s> => s[x := [[e]](s)]

<c1, s> => s'    <c2, s'> => s''
---------------------------------
<c1; c2, s> => s''
```

### Hoare Logic

Hoare logic provides formal reasoning about program correctness.

**Judgment:**
- `{P} c {Q}`
- where P, Q are assertions (predicates on states), c is a command

**Inference rules:**

**Assignment:**
```
------------------------
{P[x := e]} x := e {P}
```

**Sequence:**
```
{P} c1 {Q}    {Q} c2 {R}
--------------------------
{P} c1; c2 {R}
```

**Conditional:**
```
{P AND b} c1 {Q}    {P AND NOT b} c2 {Q}
------------------------------------------
{P} if b then c1 else c2 {Q}
```

**While (loop invariant I):**
```
{I AND b} c {I}
----------------------------------
{I} while b do c {I AND NOT b}
```

**Consequence:**
```
P' IMPLIES P    {P} c {Q}    Q IMPLIES Q'
--------------------------------------------
{P'} c {Q'}
```

**Example:**
- `{x = 5} y := x + 1 {y = 6}`
- `{x > 0} while x > 0 do x := x - 1 {x = 0}`

### Temporal Logic

Temporal logic formalizes reasoning about execution sequences.

**Key operators over infinite state sequence:**
- `NEXT phi` - phi holds in next state
- `EVENTUALLY phi` - phi holds in some future state
- `ALWAYS phi` - phi holds in all future states
- `phi1 UNTIL phi2` - phi1 holds until phi2 holds

**Example:**
- `ALWAYS (request IMPLIES EVENTUALLY response)`

---

## Chapter 7: From Animation to Imitation

### Transition Systems

A transition system models state evolution through discrete transitions.

**Formal definition:**

`(S, s0, T)` where:
- S - set of states
- s0 - initial state
- `T SUBSET S CROSS S` - transition relation

A transition `(s, s') IN T` means the system can move from state s to s'.

**Execution:** sequence of states `s0, s1, s2, ...` where `(si, si+1) IN T`

### Semiautomata

A semiautomaton adds an input alphabet to a transition system.

**Formal definition:**

`(S, A, delta)` where:
- S - set of states
- A - input alphabet
- `delta: S CROSS A -> S` - transition function

Reading input string `a1 a2 ... an` produces state sequence:
- `s0, delta(s0, a1), delta(delta(s0, a1), a2), ...`

### Finite State Automata

A finite state automaton recognizes languages by accepting or rejecting strings.

**Formal definition:**

`(Q, A, delta, q0, F)` where:
- Q - finite set of states
- A - input alphabet
- `delta: Q CROSS A -> Q` - transition function (deterministic)
- q0 - initial state
- `F SUBSET Q` - set of accept states

**Deterministic FSA (DFA):**
- `delta(q, a)` gives exactly one next state

**Nondeterministic FSA (NFA):**
- `delta: Q CROSS A -> P(Q)` (power set of Q)
- `delta(q, a)` gives a set of possible next states

**Equivalence:**

For every NFA, there exists an equivalent DFA recognizing the same language.

Construction: subset construction (states of DFA = subsets of NFA states)

### Pushdown Automata

A pushdown automaton extends FSA with a stack for unbounded memory.

**Formal definition:**

`(Q, A, G, delta, q0, Z0, F)` where:
- Q - finite set of states
- A - input alphabet
- G - stack alphabet
- `delta: Q CROSS (A UNION {EPSILON}) CROSS G -> P(Q CROSS G*)`
- q0 - initial state
- Z0 - initial stack symbol
- F - accept states

**Deterministic PDA:**
- At most one transition for each `(q, a, Z)` configuration

**Nondeterministic PDA:**
- Multiple possible transitions

Note: Deterministic and nondeterministic PDAs are NOT equivalent.

Nondeterministic PDAs recognize exactly the context-free languages.

### Turing Machines

A Turing machine has an infinite tape for unbounded memory.

**Formal definition:**

`(Q, A, G, delta, q0, qaccept, qreject)` where:
- Q - finite set of states
- A - input alphabet
- G - tape alphabet `(A SUBSET G)`
- `delta: Q CROSS G -> Q CROSS G CROSS {L, R}` - transition function
- q0 - initial state
- qaccept - accept state
- qreject - reject state

**Operation:**

Read symbol on tape, write symbol, move head left (L) or right (R)

**Deterministic Turing Machine (DTM):**
- delta is a function (single transition per configuration)

**Nondeterministic Turing Machine (NTM):**
- `delta: Q CROSS G -> P(Q CROSS G CROSS {L, R})`
- Multiple possible transitions

**Equivalence:** DTM and NTM recognize the same languages.

### Universal Turing Machine

A universal Turing machine U can simulate any Turing machine M.

**Input to U:**
- `<M, w>` (encoded string representing machine M and input w)

**Operation:**
- U simulates M on input w
- U accepts if M accepts w
- U rejects if M rejects w

### Church-Turing Thesis

All reasonable models of computation are equivalent in power.

**Equivalences:**
- Turing machines EQUIV Lambda calculus

These models compute exactly the same class of functions (computable functions).

The thesis is not a mathematical theorem but a claim about the nature of computation, supported by evidence: every computational model proposed has proven equivalent to Turing machines.

---

## Chapter 8: From Universality to Impossibility

### Asymptotic Complexity

Asymptotic notation describes function growth as input approaches infinity.

**Big-O (upper bound):**

`f(n) = O(g(n))` if:
- `EXISTS c IN R+, EXISTS n0 IN N, FOR ALL n >= n0: f(n) <= c * g(n)`

**Big-Omega (lower bound):**

`f(n) = Omega(g(n))` if:
- `EXISTS c IN R+, EXISTS n0 IN N, FOR ALL n >= n0: f(n) >= c * g(n)`

**Big-Theta (tight bound):**

`f(n) = Theta(g(n))` if:
- `f(n) = O(g(n)) AND f(n) = Omega(g(n))`

### Time Complexity

Time complexity measures computational steps as a function of input size.

**Common classes:**
- `O(1)` - constant
- `O(log n)` - logarithmic
- `O(n)` - linear
- `O(n log n)` - linearithmic
- `O(n^2)` - quadratic
- `O(2^n)` - exponential

### Space Complexity

Time complexity measures computational steps as a function of input size.

**Common classes:**
- `O(1)` - constant
- `O(log n)` - logarithmic
- `O(n)` - linear
- `O(n^k)` - polynomial

### Complexity Classes: P

P is the class of decision problems solvable in polynomial time.

**Formal definition:**
- `P = {L | EXISTS deterministic Turing machine M and polynomial p(n) such that M decides L in time O(p(n))}`

### Complexity Classes: NP

NP is the class of decision problems verifiable in polynomial time.

**Formal definition:**
- `NP = {L | EXISTS nondeterministic Turing machine M and polynomial p(n) such that M decides L in time O(p(n))}`

**Equivalently:**
- `NP = {L | EXISTS polynomial-time verifier (algorithm) V such that: x IN L IFF EXISTS certificate c where V(x, c) accepts (returns YES)}`

### Polynomial-Time Reductions

A polynomial-time reduction (Karp reduction) transforms one problem to another in polynomial time.

**Formal definition:**

A KARP-REDUCES-TO B if:
- `EXISTS polynomial-time computable function f such that:`
- `FOR ALL x: x IN A IFF f(x) IN B`

**Properties:**
- If A KARP-REDUCES-TO B and B IN P, then A IN P
- If A KARP-REDUCES-TO B and A NOT IN P, then B NOT IN P

### NP-Completeness

A problem is NP-complete if it is in NP and every problem in NP reduces to it.

**Formal definition:**

Problem L is NP-complete if:
1. L IN NP
2. FOR ALL A IN NP: A KARP-REDUCES-TO L

**Properties:**
- If any NP-complete problem is in P, then P = NP
- If any NP-complete problem is not in P, then P != NP

### P vs NP Problem

The P vs NP problem asks whether P = NP.

**Question:**

Can every problem whose solution is verifiable in polynomial time also be solved in polynomial time?

**Status:** Open problem

### Turing Reductions

A Turing reduction solves one problem using another as an oracle (a hypothetical subroutine that solves a problem in one step).

**Formal definition:**

A TURING-REDUCES-TO B if:
- EXISTS Turing machine M with oracle for B such that M decides A

**Properties:**
- If B is decidable, then A is decidable
- If A is undecidable, then B is undecidable

### Halting Problem

The halting problem: Given encoding of Turing machine M and input w, does M halt on w?

**Input:** `<M, w>` where `<M>` is the encoding of a Turing machine M, w is a string

**Output:** YES if M halts on w, NO if M loops forever

**Theorem:** The halting problem is undecidable.

**Proof by contradiction:**

Assume H decides the halting problem (H accepts `<M, w>` iff M halts on w).

Construct machine D:
```
D(<M>):
  if H(<M>, <M>) accepts: loop forever
  if H(<M>, <M>) rejects: halt
```

Consider `D(<D>)`:
- If `H(<D>, <D>)` accepts: `D(<D>)` loops (contradiction)
- If `H(<D>, <D>)` rejects: `D(<D>)` halts (contradiction)

Therefore H cannot exist.

### Rice's Theorem

All non-trivial properties of Turing machine languages are undecidable.

**Theorem:**

Let P be a property of languages such that:
1. P is non-trivial
2. P depends only on the language

Then "Does Turing machine M recognize a language with property P?" is undecidable.

### Gödel's Incompleteness Theorems

In any consistent formal system capable of expressing arithmetic, there exist true statements that cannot be proven.

Incompleteness in logic and undecidability in computation both reveal fundamental limits of formal systems.