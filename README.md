# Is the randomized Tukey HSD really low-powered?

This repository investigates the statistical power of the randomized Tukey
HSD, as implemented by
[RandomizedTukeyHSD.jl](https://github.com/jerome-white/RandomizedTukeyHSD.jl),
and of the alternatives to it. It is a living document: it states what we
are trying to find out, what we currently believe and how well we know it,
and what we will do next. It is written to guide contributors, human and
agent, and also to teach the background to its maintainer. If a section
here disagrees with a result, the section is wrong and gets edited.

How to contribute is in [CONTRIBUTING.md](CONTRIBUTING.md). The plan is a
set of issues grouped into families; see [The plan](#7-the-plan).

## How to read the status labels

Every claim below carries a status, defined in `CONTRIBUTING.md`:

- **[verified]**: checked against a source that was read, or against code
  that was read, or against a recorded run.
- **[reasoned]**: follows from an argument, but nobody has checked it
  against data or a source.
- **[unverified]**: from memory, conversation or hearsay.

When in doubt, a claim takes the lower label. A claim moves up only with
evidence, cited where the claim is made.

## Contents

1. [The question](#1-the-question)
2. [Primer: power and error rates](#2-primer-power-and-error-rates)
3. [The test and its FWER control](#3-the-test-and-its-fwer-control)
4. [Why it may lose power](#4-why-it-may-lose-power)
5. [Reservations about the comparison](#5-reservations-about-the-comparison)
6. [Candidate remedies](#6-candidate-remedies)
7. [The plan](#7-the-plan)
8. [Experimental design](#8-experimental-design)
9. [Open questions](#9-open-questions)
10. [Decision log](#10-decision-log)
11. [Glossary](#11-glossary)
12. [References](#12-references)

## 1. The question

When evaluating responses to treatments over a collection of blocks it
is natural to ask how treatments differ, and whether those differences
are the result of real effects. In the case of information retrieval
(IR; the focus of this work) that means studying the differences
between scores of systems across a range of topics. Researchers often
want to know whether differences in scoring are statistically
significant.

With `k` systems there are `k(k-1)/2` system pairs. Testing every pair
at 5% makes at least one false alarm very likely. The probability of
that is the family-wise error rate (FWER). The randomized Tukey HSD
(RTHSD) is designed to keep it at 5% across all pairs at once. It is
thus a technique that has been suggested as a modern approach to IR
evaluation ([Carterette 2012](https://doi.org/10.1145/2094072.2094076)
and [Sakai 2018](https://doi.org/10.1007/978-981-13-1199-4_4)).

Otero, Parapar and Barreiro
([ECIR 2025](https://doi.org/10.1007/978-3-031-88711-6_27); we read
[arXiv v1](https://arxiv.org/abs/2501.03930)) compared RTHSD with other
statistical tests in two experiments. The first is a simulation built
from TREC-8 runs: 10, 30 or 50 topics, 3, 5 or 10 systems, and either
every hypothesis true or every hypothesis false. Its measure is complete
power, the probability of rejecting all of the false hypotheses. The
second uses Million Query 2009 data and measures average power, the
share of truly different pairs found, on random subsets of 10 to 100
topics. They report that while RTHSD does well with respect to FWER, it
is underpowered. This excerpt closes the second experiment:

> We also highlight that the randomised TukeyHSD, which yields
> remarkably good results in terms of expected FWER, performs poorly
> in terms of power.  For sample sizes larger than 50, is indeed the
> worst test. In light of these results, we recommend researchers to
> avoid employing this test, since it likely leave out interesting and
> significant results.

The simulation reaches the same conclusion: they describe RTHSD there
as having "very low power". They ultimately recommend researchers avoid
RTHSD in favour of the Wilcoxon signed-rank test with Benjamini-Hochberg
correction. [verified: arXiv v1 read in full, text and figures; the
published version was not compared.]

In light of this finding, we want to know three things:

1. **Does the low power hold more broadly?** The paper measured complete
   power with every hypothesis either true or false, and tested up to 100
   topics. Does the result hold under other definitions of power, with a
   mix of true and false hypotheses, with more topics, and against
   baselines that make the same kind of error guarantee (section 8)?
2. **If so, why?** Is it intrinsic to controlling the error rate this
   strictly, or does it come from choices in this particular procedure?
3. **Can it be fixed?** Is there a variant that keeps the guarantee and
   recovers power, and what does it cost?

## 2. Primer: power and error rates

### Two ways to be wrong

A test looks at data and says "different" or "no evidence of a
difference". Either answer can be wrong:

|                           | Truth: no difference | Truth: real difference |
|---------------------------|----------------------|------------------------|
| Test says "different"     | Type I error (alpha) | Correct (1 - beta)     |
| Test says "no difference" | Correct (1 - alpha)  | Type II error (beta)   |

Each cell is the probability of that outcome given the column's truth.

- **alpha** is the probability of a false alarm when there is no
  difference. We choose it in advance, usually 0.05.
- **beta** is the probability of missing a difference that is real.
- **Power** is `1 - beta`: the probability of detecting a real
  difference.

Each column sums to one. alpha and beta do *not* sum to one, and there is
no fixed relation between them. What does hold is a trade-off: for a given
test, effect size, amount of noise and number of topics, demanding a lower
alpha makes the test stricter and lowers its power. Power also grows with
the true size of the difference, with more topics, and with less noise.
[verified: standard; any statistics text.]

### The ROC view

Change the threshold of a single test and record two numbers: the false
alarm rate (alpha) and the power (the rate of detecting real differences).
Plotted as alpha against power, the points form the ROC curve of the test.
A better test has a curve that lies above another's. [verified: standard.]

### Many tests at once

With 45 pairs, if each is tested at 0.05 and the pairs were independent,
the chance of at least one false alarm when nothing differs is
`1 - 0.95^45`, about 0.90. Pairs are not independent, so the real figure
differs, but the problem is real. [reasoned: arithmetic, with the
independence caveat.] There are two main ways to define "false alarm"
across many tests:

- **FWER, the family-wise error rate**: the probability of *at least one*
  false alarm in the whole family of tests. Controlling it at 0.05 means
  that in 95% of experiments with nothing to find, you find nothing.
  Tukey HSD, Bonferroni and Holm control it.
- **FDR, the false discovery rate**: the *expected fraction* of your
  discoveries that are false. Controlling it at 0.05 says that, on
  average, 5% of what you report is wrong, not that you will likely report
  nothing wrong. Benjamini-Hochberg controls it.

FDR is the weaker promise, so procedures that control it are allowed more
power. FDR is, in expectation, `1 - precision` of the list of discoveries,
so FDR against power resembles the precision-recall view. FWER is a
probability of "any error" and has no counterpart in precision and recall.
[reasoned for the FDR and precision statement; standard for the rest.]

Comparing a FWER method with an FDR method on power alone is therefore not
like for like: the FDR method is allowed to be bolder. This is one of the
reservations in section 5.

### Plotting error against power

The ROC idea extends to procedures for many tests: for each procedure,
vary its nominal level, measure the *achieved* FWER (or FDR) and the power
in simulation, and plot achieved error on the x-axis against power on the
y-axis. Each procedure gives a curve, and a better procedure sits above the
others. This lets us compare procedures at the same achieved error rate,
which is the like-for-like way to compare them. [reasoned: a design for
the experiments, not yet done.]

### Ways to count power across many pairs

"Power" is not one number when some pairs differ and some do not:

- **Per-pair (average) power**: the fraction of truly different pairs
  that the procedure finds.
- **Any-pair power**: the probability of finding at least one of the
  truly different pairs.
- **Complete power**: the probability of finding *all* of them.

These can rank procedures differently, and complete power is very
demanding when some pairs differ only slightly. All three are reported
(section 8).

## 3. The test and its FWER control

### What the library computes

Scores form a matrix of topics (blocks) by systems (treatments). The
procedure [verified: read in `src/RandomizedTukeyHSD.jl` and
`src/solve.jl` of RandomizedTukeyHSD.jl]:

1. Compute the mean score of each system and, for each pair, the absolute
   difference of means.
2. Repeat many times (the trials): shuffle each topic's scores across the
   systems, so that each topic keeps its set of scores but their
   assignment to systems is random; compute the system means; record the
   **range** of the means, `maximum(means) - minimum(means)`.
3. A pair's p-value is the fraction of trials whose range reached the
   pair's observed absolute difference. There is no "+1" in the count, so
   a p-value can be zero; the library reports its Monte Carlo uncertainty
   instead of changing the convention.

Shuffling within a topic keeps the topic's difficulty, which is the same
for every system, so the topic effect cancels in the differences.
[reasoned.]

### Why this controls FWER

Suppose all systems are interchangeable (exchangeable within each topic).
Then the observed range of means `R` has the same distribution as the
range in a shuffled trial. Define the global p-value as the fraction of
trials with range at least `R`. A pair is rejected at level alpha only if
its p-value is at most alpha, and every pair's observed difference is at
most `R`, so every pair's p-value is at least the global p-value. So *any*
rejection implies that the global p-value is at most alpha, which happens
with probability at most alpha. Hence the chance of any false alarm is at
most alpha. [reasoned: the argument is the standard one for max-statistic
permutation tests; it has not been checked here by simulation.]

Two consequences:

- The p-values are **already adjusted for multiplicity**. They must not be
  fed to Holm, Bonferroni or Benjamini-Hochberg a second time. [reasoned.]
- The shuffle takes the threshold from the largest difference among *all*
  systems. The next sections build on that fact.

### What the argument does not cover

The argument above assumes *every* system is interchangeable. It is called
**weak control**. **Strong control** is the requirement that the guarantee
holds when some systems really differ and others do not, which is the
usual situation. Classical Tukey HSD has strong control under its
normal-data assumptions. For the randomized version the shuffle also mixes
in the systems that truly differ, so the shuffled distribution is not the
distribution of the null systems alone, and whether the guarantee still
holds is not obvious. In the multiple-testing literature this is the
question of "subset pivotality". The published comparison also did not
test mixtures of true and false nulls. [unverified: open question; a
simulation with mixed nulls answers it for a given setting.]

A related concern: the shuffle assumes systems are exchangeable *in
distribution*, not only in mean. Systems with different variances or
different score shapes (for example bounded metrics near 0 or 1) violate
that even when their means are equal, and the p-values then answer a
slightly different question. [reasoned.]

## 4. Why it may lose power

Each item is a hypothesis to test, not a finding.

1. **The guarantee is the strictest one.** Controlling the chance of any
   false alarm across all pairs costs power compared with controlling the
   fraction of false discoveries. Part of the gap in the published
   comparison may be just this price. It can explain the distance to
   Benjamini-Hochberg, but not the distance to Holm and Bonferroni, which
   pay the same price (section 5). The cost grows with the number of
   systems (45 pairs at 10 systems, 1,225 at 50). [reasoned.]
2. **The procedure is single-step.** A *step* is one round of decision.
   Single-step means a single threshold, taken from the range over all
   systems, is applied to every pair at once. One weak or noisy system
   widens the range of the shuffled means for everyone and so costs power
   on pairs that are clearly different. [reasoned.] Step-down procedures
   work in rounds: they test the most extreme difference first, and if it
   is rejected they test the remaining hypotheses against a threshold
   computed from a smaller set. At the same FWER they can be more
   powerful, never less. In the Newman-Keuls and REGW style, sorted
   systems are compared against the range of the group they span; in the
   Westfall-Young style, thresholds are recomputed from the hypotheses not
   yet rejected. Which form fits a permutation test on this design is a
   design question. [unverified: from memory.]
3. **The statistic is not studentized.** The library's statistic is the
   range of raw means. A *studentized* statistic divides each difference by
   an estimate of its standard error (a paired difference divided by its
   standard deviation over topics, scaled by the square root of the number
   of topics), so noisy pairs and quiet pairs are put on one scale. This
   is not z-scoring a system's scores, which would remove the very means
   being compared. With unequal variances across systems, the raw range
   is dominated by the noisiest systems and loses sensitivity on the rest.
   [reasoned.]
4. **Complete power is harsh.** It counts a run as a success only if every
   truly different pair is found, so one hard pair decides the outcome, and
   it penalizes FWER methods the most. [reasoned.]
5. **The statistic uses means.** In the paper's Figure 2 the
   Wilcoxon-based tests beat the t-test-based ones at every correction, by
   a wide margin. If per-topic differences in scores such as average
   precision are skewed, means lose power relative to ranks, and a
   procedure built on the range of means inherits that. [reasoned: from
   the figure and general knowledge; the mechanism is untested.]

Whether any of these is the main cause is unknown. The experiments are
built to separate them.

## 5. Reservations about the comparison

These are objections to be tested, not conclusions. The details of the
paper come from an earlier reading and were not rechecked here
[unverified], and none of the points has been shown to change the ranking.

1. **FWER against FDR.** The winner they recommend controls FDR, a weaker
   promise (section 2), and the unadjusted Wilcoxon test controls neither.
   Fairer baselines include FWER-controlling methods such as Holm and
   Bonferroni, and the comparison should be made at equal achieved error.
2. **Complete power** is the measure in their simulations (section 2).
3. **No mixed nulls.** Their simulations had either all hypotheses true
   or all false, as recorded. The mixed case is the one where FWER and FDR
   differ most.
4. **A truth label that favours liberal tests.** In their large check
   (Million Query 2009), a pair counted as truly different if its mean
   differed by more than a very small margin, which made nearly all pairs
   "different". Rejecting more then looks like being right.
5. **Scale.** Their simulations used 10, 30 or 50 topics and 3, 5 or 10
   systems, built from TREC-8 (50 topics, 129 runs), with 100,000
   permutations. These sizes are normal for classic lab experiments and
   small for modern collections. The only large check, 687 topics, was
   reported qualitatively. The number of permutations is not a concern:
   the Monte Carlo error of a p-value near 0.05 at 100,000 permutations is
   about 0.0007. [reasoned: arithmetic.]

A note on topics. More topics give every sound test more power in absolute
terms. The paper's claim is about *ranking*: that the randomized Tukey HSD
falls behind the others beyond 50 topics. These do not contradict each
other, but we cannot tell from the text which holds at large sizes. One of
the experiments measures it directly. [unverified.]

## 6. Candidate remedies

Each is a candidate to implement and measure, with the cost it brings.

- **Step-down variant** (section 4, item 2): same guarantee, no less power.
  Cost: new procedure, and a check that it keeps strong control in this
  design.
- **Studentized statistic** (item 3): expected to help most when system
  variances differ. Cost: the null distribution of the studentized range
  must be generated by the same shuffle, which is more work per trial.
- **Restrict the family**: if only comparisons against a baseline matter,
  test only those. This is a choice about the question asked, and it
  reduces the multiplicity penalty at no cost in validity.
- **An FDR variant**: compute a *per-pair* permutation p-value from that
  pair's own difference (not the adjusted p-value of the range procedure)
  and apply Benjamini-Hochberg. Cost: it is close to Wilcoxon plus
  Benjamini-Hochberg with a mean-difference statistic, and arguably no
  longer a Tukey procedure. The library does not produce these p-values
  today, because it stores only counts of the range.
- **More topics**: a property of the data, not of the method.

Whatever else is found, a FWER-controlling method will not match an
FDR-controlling method in power at equal nominal level. [reasoned.] The
useful result is how much of the gap is the guarantee and how much is the
procedure.

## 7. The plan

The plan is a set of issues, grouped into families. A family is a parent
issue; its sub-issues carry it out. Families, in order:

- **A. Protocol**: the questions, hypotheses, metric definitions and the
  rule for deciding what counts as a win, written down before experiments.
- **B. Simulated data**: a generative model (topic effect, system effect,
  bounded and skewed noise, unequal variances), null structures (global,
  mixed, exactly-equal groups), seeds shared across methods, and checks of
  the generator itself.
- **C. Methods**: each method behind one common interface. Single-step
  randomized Tukey HSD (through the library), Wilcoxon with
  Benjamini-Hochberg, paired t-test with Holm, per-pair permutation with
  Benjamini-Hochberg, step-down maxT, studentized randomized Tukey HSD.
- **D. Experiments**: first a replication of the published setup; then
  calibration of error rates under global and mixed nulls; then power
  against effect size, number of topics and number of systems; then the
  comparison at equal achieved error.
- **E. Visualization and reporting**: error-power curves, tables, and a
  written summary of results.
- **F. Infrastructure**: Julia project layout, reproducibility, and run
  scripts for long jobs.
- **G. Upstream**: proposals back to RandomizedTukeyHSD.jl for anything
  that earns its place.

Issue numbers will be linked here once the issues exist.

Replication comes first. If the published numbers cannot be reproduced,
nothing downstream can be trusted.

## 8. Experimental design

Principles that keep a comparison like for like, which every experiment is
held to:

1. **Known truth.** Data are generated from a model in which we know which
   systems are equal and which differ.
2. **Mixed nulls.** Include scenarios where only some systems differ, and
   groups of exactly equal systems.
3. **Realistic scores.** Include bounded, skewed scores and unequal
   variances across systems, not only normal noise.
4. **Common random numbers.** Every method sees the same simulated data for
   the same seed.
5. **Error first, then power.** Measure achieved FWER and FDR for every
   method before comparing power. Compare power at equal achieved error and
   also against each method's own claimed guarantee.
6. **Several power measures**: per-pair, any-pair and complete.
7. **Vary the design**: topics from 10 to several hundred, systems from 3
   to 10 and beyond, small and large effects.
8. **Enough replications** that the Monte Carlo error of each estimate is
   small and reported.
9. **The plan is fixed before the results are seen.** Changes after the
   fact are recorded with their reason.

Practical notes:

- A "truly different" pair needs a definition. With continuous system
  effects every nonzero difference is real but possibly undetectable.
  Prefer exactly-equal systems for nulls and clearly separated systems for
  alternatives, and treat the in-between as its own scenario.
- Each simulated dataset costs thousands of shuffles per method. Large
  grids exceed what an agent may run (see `CONTRIBUTING.md`); they are
  delivered as a script with the command, expected duration and output path.
- The library gives per-pair p-values from the single-step procedure, so
  that method needs no change to the library. Step-down, studentized and
  per-pair variants are prototyped here first, and move upstream only if
  they earn it.

## 9. Open questions

- Does the randomized Tukey HSD keep strong control (mixed nulls)?
- How much of the published power gap remains once the comparison is made at
  equal achieved error?
- How does power behave beyond 50 topics, up to several hundred?
- Does a step-down variant recover power? Which form is valid for
  permutation tests in this design?
- Does studentizing help, and when?
- Which data should the replication use? The published simulations were
  built from TREC-8 runs; whether the runs can be obtained and used is
  unchecked.

## 10. Decision log

| Date       | Decision                                                      |
|------------|---------------------------------------------------------------|
| 2026-11-13 | The implementation language is Julia.                         |
| 2026-11-13 | Replicate the published setup before anything else.           |
| 2026-11-13 | Adopt the conventions of RandomizedTukeyHSD.jl (#1).          |
| 2026-11-13 | Prototype new procedures here; upstream only if useful.       |

## 11. Glossary

- **Block / topic**: one unit scored by every system.
- **Treatment / system**: one of the things being compared.
- **Range**: the largest system mean minus the smallest.
- **Exchangeable**: interchangeable, so that shuffling them does not change
  the distribution.
- **FWER**: probability of at least one false alarm in the family of tests.
- **FDR**: expected fraction of discoveries that are false.
- **Weak control**: the error guarantee holds only when all hypotheses are
  true. **Strong control**: it holds whichever hypotheses are true.
- **Single-step / step-down**: one threshold for all hypotheses / thresholds
  recomputed in rounds.
- **Studentized**: divided by an estimate of its own standard error.
- **Common random numbers**: giving every method the same random data.

## 12. References

All cited from memory or from an earlier session; none re-read for this
document. [unverified]

- Otero, Parapar, Barreiro. ECIR 2025. https://arxiv.org/abs/2501.03930
- Carterette. Multiple testing in statistical analysis of systems-based
  information retrieval experiments. ACM TOIS, 2012.
- Sakai. Laboratory Experiments in Information Retrieval: Sample Sizes,
  Effect Sizes, and Statistical Power. Springer, 2018.
- Tukey. The problem of multiple comparisons, 1953.
- Benjamini, Hochberg. Controlling the false discovery rate. JRSS B, 1995.
- Holm. A simple sequentially rejective multiple test procedure, 1979.
- Westfall, Young. Resampling-Based Multiple Testing, 1993.
