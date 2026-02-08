# Chapter 1: From Particularity to Universality

## Propositional Logic

A proposition is a statement that is either true or false.

**Connectives:**

- `P AND Q` — both P and Q are true  
- `P OR Q` — at least one of P or Q is true  
- `NOT P` — P is false  
- `P IMPLIES Q` — if P is true, then Q is true  
- `P IFF Q` — P and Q have the same truth value  

**De Morgan's Laws:**

- `NOT (P AND Q) IFF (NOT P) OR (NOT Q)`  
- `NOT (P OR Q) IFF (NOT P) AND (NOT Q)`  

## First-Order Logic

Extends propositional logic with quantifiers over objects.

**Quantifiers:**

- `FOR ALL x, P(x)` — P holds for every x  
- `EXISTS x, P(x)` — P holds for some x  

**Variables:**

- **Bound variable:** quantified by `FOR ALL` or `EXISTS`  
- **Free variable:** not quantified  

## Proof Techniques

**Direct proof:**  
To prove `P IMPLIES Q`, assume P is true and derive Q.

**Proof by contradiction:**  
To prove P, assume `NOT P` and derive a contradiction.

**Proof by contrapositive:**  
To prove `P IMPLIES Q`, prove `(NOT Q) IMPLIES (NOT P)`.

**Proof by cases:**  
To prove P, split into exha
