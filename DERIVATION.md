# Replacing repeated restricted sums with cardinality-weighted queries

Date: 2026-10-30. Research prototype; no claim of novelty, proof certification, or competitiveness with production #SAT solvers.

## Scope and modularity

The OpenAI manuscript has a rational-interpolation route in Sections 2-6 and a different, sharper recursive tensor-contraction route in Appendix D. This experiment modifies the Boolean-differences reconstruction of the former. It does not implement or improve the latter's stated bound.

The replaceable interface is: given valid point-specific groups and measurement rows, evaluate restricted sums of products over Boolean assignments in a coefficient ring. The upstream toy accepts either integers or truncated rational power series. We extend the values carried through the group sums to include one additive cardinality statistic.

The upstream baseline is Chelsea A. McElveen's MIT-licensed executable exposition, commit 7a73e39b11c71d2b3c7209074e04e84208612b58:
https://github.com/chelsinchaos/measurement-indexed-counting

The original mathematical identity is Lemma 3.3 (label `lem:boolean-differences` in the downloaded source; use the label if numbering changes) and the Taylor recovery in Section 3.4:
https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/main_full_version.tex

## Derivation

Let A be the set of coordinates on which any correction Delta P(t,z) depends modulo z^(q+1), and let a=|A|. For each assignment x, write p_x(X) for the truncated polynomial G_x(X,z), regarded as a polynomial in the independent Boolean inputs X over Q[z]/(z^(q+1)). Its multilinear reduction has degree at most q and uses only variables in A.

The paper's Boolean-differences identity, restricted to its actual active coordinates, gives

    sum_x p_x(x)
      = sum_{E subset A, |E|<=q} sum_{U subset E}
          (-1)^(|E|-|U|) sum_{x: x_i=1 for i in E} p_x(1_U).

Fix U and x with x_i=1 for all i in U. Let w be the number of one-bits of x in A minus U. A superset E of U compatible with x is obtained by choosing j of those w positions. Reordering finite sums gives

    sum_x p_x(x)
      = sum_{U subset A, |U|<=q}
          sum_{x: x_i=1 for i in U} p_x(1_U) B_(q-|U|)(w),

where

    B_h(w) = sum_{j=0}^{min(h,w)} (-1)^j binom(w,j).

Its closed form is 1 if w=0, zero if 1<=w<=h, and (-1)^h binom(w-1,h) if w>h. No approximation or sampling is used.

For a fixed U, the weight w is additive across disjoint variable groups. Therefore a group sum can store a histogram indexed by its contribution to w; multiplication of group sums becomes histogram convolution. The final weighted sum applies B_h to the histogram. All dependence of the original transformed factors on other groups still passes through the same measurements; the added cardinality statistic does not invalidate their separation.

## Query count and limits

The original number of restricted calls per point is

    C(n,q) = sum_{j=0}^{min(n,q)} binom(n,j) 2^j.

Restricting to active coordinates changes n to a. Regrouping then uses

    S(a,q) = sum_{j=0}^{min(a,q)} binom(a,j)

weighted calls. These weighted calls are more expensive than scalar calls. A histogram has at most a+1 entries, each containing a truncated z-series. Convolution and coefficient bit lengths add polynomial overhead. The reduction is not a claim of a C/S wall-clock speedup.

When q>=a the query counts are 3^a versus 2^a. When q is much smaller than a, the leading query-count ratio is roughly 2^q. The dominant q log(a/q) part of the exponent remains; this does not improve the paper's advertised n/(log n)^(1/2+o(1)) scale.

The prototype caches local assignment tables, like the upstream toy. Those caches can be exponential. The paper's polynomial-space guarantee must not be ascribed to this implementation. The identity itself also admits streamed local enumeration with polynomial-size histogram accumulators, but that version has not been benchmarked here.

## A smaller alternative to full histograms

Only B_h(w) is needed, and

    B_h(w) = [s^h] (1-s)^w / (1-s).

For each U, instead of keeping all w values, carry the coefficients of (1-s)^w through order h=q-|U|. Multiplication across groups is convolution truncated at order h. Sum those coefficients at the end. This requires h+1 entries rather than a+1. At h=0 it is just the ordinary restricted sum. This moment representation is implemented as `moment_reconstruction`.

High h can make moments more expensive than sparse histograms; neither representation should be chosen unconditionally. Any adaptive rule needs to count arithmetic and memory costs, not just calls.

## Exact cofactor pruning when q covers the active support

For a query U, every possible weight satisfies 0 <= w <= a-|U|. If q>=a,
then h=q-|U| >= a-|U|, and the binomial identity above gives B_h(w)=0
for every w>0. Only assignments with x_i=0 for i in A minus U can
contribute. Together with the existing x_i=1 constraints on U, this fixes
the active coordinates to the indicator of U.

`cofactor_reconstruction` applies those zero constraints while constructing
local measurement buckets, before any tuple or clause evaluation. Nonnegative
group weights sum to zero exactly when each is zero, so this local filtering
preserves the global sum. If q<a, it uses the original histogram unchanged.
The identity works in the rational truncated ring and in the valid finite
fields used by the modular backend.

This removes work that would later cancel. It does not reduce the number of
U queries. If a=n, it partitions all 2^n assignments and retains polynomial
evaluation/interpolation overhead. If a is nearly n, as in these benchmarks,
the active-subset enumeration is still nearly exhaustive. A practical advantage
would require genuinely cheaper treatment of remaining variables or a different
regime; making n larger alone does not supply that advantage.

## Verification in this package

`python verify.py` performs 120 comparisons of the original, active-support, histogram, and moment reconstructions against an independent direct expansion of H(t+z), and 94 modular arithmetic checks. Cases include different seeds/designs, small Taylor orders, empty formulas, unused variables, and contradictions. Finite checks accompany the finite-sum derivation; they do not replace it.

It additionally performs 193 cofactor checks: 60 independent coefficient
cases in rational arithmetic and two finite fields (180 comparisons), nine
empty/constant/contradictory edge cases, and four full modular counts. The
coefficient cases exercise both the pruning branch (38) and the unchanged
fallback (22). `python cofactor_checks.py` also audits all 88 modular Taylor
coefficients across the five benchmark cases and records the evaluation
counters. See `cofactor-audit.json` and `cofactor-results.json`.

The portable benchmark compares exact full counts against brute force and, when installed, dd.autoref. Every timed paper-based run includes fresh design construction, point preparation, reconstruction, and final recovery. Three repeats rotate method order by default. Optional memory instrumentation runs separately and measures Python allocations only. See README.md and sample-results.json for the actual measured scope. The unmodified upstream module's attribution is retained in third_party/LICENSE.txt. No full Lean check or production-solver superiority is claimed.

## Components worth replacing next

| Component | Candidate replacement | Requirement / failure mode |
|---|---|---|
| Repeated restricted Taylor sums | Active-coordinate pruning; cardinality histograms; truncated binomial moments | Tested here; weighted-query overhead can outweigh fewer calls |
| Rational arithmetic | Finite-field evaluation and CRT recovery | Exclude primes dividing coefficient denominators or interpolation determinants; sufficient combined modulus for final count |
| Universal schedules and parameters | Search small point-specific designs against measured cost | Preserve P(0), degree, and row-space conditions; avoid optimizing only local groups |
| Local group enumeration | Component caching or tree-decomposition / tensor contraction | Signed rational factors cannot automatically be passed as Boolean clauses to a standard #SAT solver |
| Final Hermite reconstruction | Precompute weights for H(0), avoiding recovery of every coefficient of H | Polynomial overhead improvement; does not remove exponential Taylor work |
| Full deformation on structured input | Count directly using a small boundary / separator | Useful established special-case algorithm; not evidence for a new general SAT bound |

Algebraic model counting provides the established framework for enriching scalar counts with polynomial values: https://arxiv.org/abs/1211.4475 . Tensor-network counting and contraction planning are also established: https://mccompetition.org/assets/files/2020/MCW_2020_paper_1.pdf . The brief search performed here does not establish that this specific regrouping application is new.

The baseline comparison is essential: structure that makes deformation cheap may also make a conventional method cheaper. Improvements inside this reconstruction should not be presented as general improvements to counting without that comparison.
