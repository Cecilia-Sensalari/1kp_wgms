# The B/C/D selection rule, explained out loud

Notes for explaining the outgroup-term constraint in conversation. The measurements are in
[`ortholog_ks_path_distance_calibration.md`](ortholog_ks_path_distance_calibration.md) (5.5 and
5.6) and the implementation in `make_branch_subtrees.py`; this is only the order in which the
argument lands when spoken, with the numbers worth having ready.


## The thirty-second version

> We were checking one of three numbers. The rate correction needs three ortholog Ks
> measurements, and one of them is subtracted. The rule only ever looked at the first. It turns
> out the missing check costs nothing to add, because the thing it constrains is a fixed offset
> from the thing we were already constraining. What it mostly buys us is knowing which branches
> are shaky — it rarely finds a better species to pick. And where it fires, the date really does
> move: by tens of percent.


## The build-up, in four steps

**1. The correction uses three numbers, not one.**

To correct the focal species A for rate differences, using a sister B and an outgroup O:

```
                    ┌── A        Ks_adjusted = ( Ks(A,B) + Ks(A,O) − Ks(B,O) ) / 2
            ┌───────┤                                                ^^^^^^^
            │       └── B                                            subtracted
      ──────┤
            └────────── O
```

Three ortholog distributions, each of which can be too deep to have a measurable peak. The third
enters with a minus sign and the same weight as the other two.

**2. We were only checking the first one.**

The old budget asked: is the sister close enough to A? It never asked what the sister's distance
to the *outgroup* was — even though that moves whenever the sister moves. A sister further from A
is usually also further from O.

The number to quote here, from the sister-rate test: on branch 667 every candidate passed the
budget, three of them had an outgroup term beyond Ks 2.7, and the six sisters disagreed by **28%**
about the answer. On the other five branches tested, that spread was 3–14%. Sorting all 32 sisters
by their outgroup term instead of by their distance to A:

| Ks(B,O) | n | median deviation |
|---|---|---|
| below 2.0 | 24 | 1.8% |
| 2.0–2.7 | 5 | 2.9% |
| **above 2.7** | 3 | **9.1%** |

**3. The fix is almost free, and this is the part worth taking slowly.**

Distances add along the tree. Call `a`, `b`, `c` the branches from the split down to A, B and O:

```
path(A,B) = a + b        path(B,O) = b + c        so their difference is c − a
```

The sister's own branch, `b`, **cancels**. The gap between "how far B is from A" and "how far B is
from the outgroup" does not depend on which sister you pick — it is a property of the focal
species and the outgroup alone. So requiring the subtracted term to be measurable is just the
condition we already had, with a tighter ceiling.

The sentence to say:

> **Both conditions limit the same thing — how far the sister may sit from A.** The first comes
> from wanting a close sister, the second from needing the outgroup term to exist. A candidate has
> to satisfy both, so the budget becomes the smaller of the two.

In practice this means picking the clades **oldest first** — D, then C, then B — so that every
representative is chosen already knowing the outgroups it will be subtracted against. The
measurable range is a path of 1.0, an expected Ks of 2.7, from the calibration.

**4. Nothing is thrown away.**

Where no candidate clears the ceiling, we take the least bad one and **flag** the branch rather
than dropping it. Every branch still gets a dataset. The flags are `ok`, `saturated:D`, or
`saturated:D;partly:C` when the picks don't all fail the same way.


## Be honest about what it bought

This is the part not to oversell, and it is better coming from you than extracted by a question.

- It changed the species of only **45 of 2354 branches**. In the rest, no better candidate exists.
- **375 branches (16%)** have a sister whose subtracted term cannot be made measurable by *any*
  choice; **347 (15%)** the same for C. Every one of them was being reported as *reliable* before
  this check existed.
- So the reselection is close to free and rarely binding. **The product is the flags**: we now
  know which 16% of the tree is resting on a subtraction we cannot measure, and we did not know
  that before.

Say "rarely binding", not "ineffective" - the two are different, and the test below is what
makes the difference. Where the constraint can act, it is averting a shift of tens of percent;
where it cannot, the flag marks one.


## What the fern run added (6–10 November)

Ten branches, half flagged, half clean. Two results:

- **The flags are truthful.** Checked against the peaks actually in the Ks database: 27 of 30
  branch/role/outgroup combinations agree. The three that disagree are one group of four pairs
  sitting within 6% of the limit — predicted Ks 2.53–2.65, measured 2.57–2.90 — so the only
  errors are at the boundary, decided by the calibration's own ~8% scatter. Away from the
  boundary it is exact.
- **The flag predicts disagreement.** At matched depth, flagged branches show about twice the
  disagreement between outgroups (1.6% → 3.3–3.9%) and up to three times the bootstrap
  uncertainty (1.5% → 4.0%).

## And what the purpose-built test added (12 November)

The fern run tested the flag, not the rule — only 2 of those 10 branches got different species at
all, and recomputing those two gave a confounded answer both times, because replacing a species
changes the measurement as well as the saturation.

The clean test holds the species fixed and varies only the outgroup, within a single run. Two
datasets were built for it among the gymnosperms, each dating one divergence with several
outgroups — some whose subtracted term is comfortably measurable, some beyond the limit.

**The flagged outgroup lands 22–61% away from where the others agree, in 6 of 6 comparisons.**
The outgroups with a measurable subtracted term agree among themselves to 0.2–2.6%, so the
departure is twenty to a hundred times the method's own scatter. Nine of the ten flagged
estimates are inflated, one deflated — the damage varies in sign, so it cannot be corrected for.

Say the limitations before they are asked for: saturation is confounded with the outgroup's
lineage (the saturated ones are always flowering plants), it rests on six comparisons in one
clade, and it still says nothing about whether the substituted sister is better than the rejected
one — only that subtracting an unmeasurable peak damages the result.
