# How the matched-null test works

## The problem: clustering always answers

Give a clustering method a single smooth cloud of points and it will
still return clusters. Model-selection criteria behave the same way:
BIC, the bootstrap likelihood-ratio test, the gap statistic all pick
*some* number of groups, and on data with skew or correlation they often
pick several. So when an analysis reports that a dataset contains, say,
four personality types, the number may describe the population, or it
may describe only the shape of the data under that method. Nothing in
the output tells you which.

The matched-null test separates the two by asking a single question:
**would a dataset with the same margins and the same correlations, but
no types, have produced the same answer?**

## A null with the structure kept and the types removed

We want a reference dataset that matches the real one in every respect
that is *not* in question, and is empty of the one thing that is. Two
features are not in question and must be preserved:

- every variable’s marginal distribution (its skew, its bounds, its
  lumps at round numbers on a Likert scale);
- the correlations among variables.

The feature in question, latent groups, must be absent by construction.

A Gaussian copula does exactly this. Write the data as columns
$`X_1, \dots, X_p`$ with empirical marginal distributions
$`F_1, \dots, F_p`$ and correlation matrix $`R`$. The null twin
$`\tilde X`$ is built in two steps:

1.  draw $`Z \sim \mathcal{N}(0, \tilde R)`$, where $`\tilde R`$ is
    chosen so that the induced correlations match $`R`$;
2.  push each column back through its own empirical quantile function,
    $`\tilde X_j = F_j^{-1}\!\big(\Phi(Z_j)\big)`$, which reuses the
    observed values of variable $`j`$.

By Sklar’s theorem the joint law of $`\tilde X`$ is the Gaussian copula
with margins $`F_j`$. Three things follow. Each $`\tilde X_j`$ has
*exactly* the empirical distribution of $`X_j`$: the twin is a
reshuffling of the real values, so no marginal test can tell them apart.
The correlation matrix of $`\tilde X`$ equals $`R`$ to within sampling
error. And the joint density is unimodal, with no islands or gaps, so
**the twin contains no cluster structure by construction**. Any clusters
your pipeline reports on the twin are therefore manufactured by the
pipeline, not carried by the data.

``` r

library(matchednull)

set.seed(1)
x <- matrix(rnorm(400 * 3), 400, 3) %*%
  chol(matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3, 3))

twin <- copula_null(x)
all(sort(twin[, 1]) == sort(x[, 1]))   # margins: identical
round(cor(x) - cor(twin), 2)           # correlations: close
```

## The test

Let $`T(\cdot)`$ be your pipeline, wrapped as a function that returns
one number, typically the selected number of clusters, but it can be any
scalar measure of clustering strength. Compute the real statistic
$`t_0 = T(X)`$, then generate $`R`$ twins and compute
$`t_1, \dots, t_R`$ with $`t_r = T(\tilde X^{(r)})`$. The one-sided
$`p`$-value is

``` math
p \;=\; \frac{1 + \#\{\, r : t_r \ge t_0 \,\}}{R + 1}.
```

The verdict reads *null-like* when $`t_0`$ sits inside the interval the
twins produce, and *exceeds the null* when the real data stand out. The
null is defined at the level of the data, not the pipeline, so the same
twins work for a $`k`$-means heuristic, a published typology’s exact
workflow, or a mixture model.

``` r

matched_null_test(x, cluster_fn = function(d) mclust::Mclust(d, verbose = FALSE)$G, R = 200)
```

## Does the test have power?

A false-positive control that never fires is worthless. The test must
stay quiet on typeless data *and* fire when types are really there. The
two panels below plant real types of two kinds and increase the signal.

![Positive controls. Red points mark datasets where the test detects the
planted types (real exceeds the null); grey points are read as
null-like. Left: types that differ only in their within-component
correlation, with identical margins. Right: types separated in their
means.](figures/positive_control.png)

Positive controls. Red points mark datasets where the test detects the
planted types (real exceeds the null); grey points are read as
null-like. Left: types that differ only in their within-component
correlation, with identical margins. Right: types separated in their
means.

In the left panel the types share identical margins and differ only in
the *orientation* of their within-group correlations. A single Gaussian
copula has one correlation matrix and cannot reproduce two, so as the
within-group correlation grows the test fires. This is the regime that
ordinary mean-comparison and marginal tests miss entirely: the groups
are invisible one variable at a time and appear only in the joint
dependence. In the right panel the types are separated in their means in
the familiar way, and the test again fires as the separation grows.
Where there is structure the matched null cannot carry, the test finds
it.

## Is the test calibrated on real data?

The motivating application is the claim that personality inventories
contain a small number of latent “types.” Running the test across all
fourteen `mclust` covariance parameterizations, on several public Big
Five and HEXACO datasets, gives the specification curve below.

![Selected number of types across all fourteen mclust covariance models
(steps), against the null band the twins produce (grey). The median
verdict is null-like in every dataset.](figures/spec_curve.png)

Selected number of types across all fourteen mclust covariance models
(steps), against the null band the twins produce (grey). The median
verdict is null-like in every dataset.

The number of types the pipeline selects (the step function) sits inside
the band the twins produce (grey) for the great majority of
specifications, and the median verdict is null-like in every dataset.
The handful of specifications that climb above the band are the most
flexible covariance models, which read the data’s skew and correlation
as extra components; the twins, which share that skew and correlation,
climb with them. The apparent “types” are the shape of the data under
`mclust`, not latent kinds of people.

## A second opinion: taxometrics

The mixture-model result lines up with an older tradition. The
Comparison Curve Fit Index (CCFI) from taxometric analysis scores
whether a construct is better described as categorical (taxonic) or
continuous (dimensional), with values below the midpoint favouring a
dimensional structure.

![Comparison Curve Fit Index by trait for NEO-120 and HEXACO. Values
fall on the dimensional side of the midpoint
throughout.](figures/ccfi.png)

Comparison Curve Fit Index by trait for NEO-120 and HEXACO. Values fall
on the dimensional side of the midpoint throughout.

Every trait, in both inventories, falls on the dimensional side. Two
methodological traditions that rarely meet, mixture-model cluster
counting and taxometrics, return the same verdict, which is reassuring
precisely because their assumptions differ.

## Heavy tails: the t-copula stress test

A verdict of *exceeds the Gaussian null* licenses only the conclusion
that the data carry structure beyond their margins and correlations.
That structure need not be types. Heavy-tailed dependence, where extreme
values across variables tend to arrive together, also exceeds a Gaussian
null, and real questionnaire and clinical data are heavy-tailed more
often than not. To keep the two readings apart, rerun the test with
t-copula twins: same margins, same correlations, but tails that co-move.

``` r

matched_null_test(x, cluster_fn, R = 200, copula = "t", df = 8)  # moderate tails
matched_null_test(x, cluster_fn, R = 200, copula = "t", df = 3)  # heavy tails
```

The ladder reads simply. A result that exceeds the Gaussian twins *and*
the t twins is hard to attribute to tails. A result the t twins
reproduce was tail dependence, not types. Either way the verdict is
sharper than a single null could give.

## Where the test can be fooled

The one regime the null cannot flag is types so separated that the
separation shows up in the margins themselves, turning a variable
visibly bimodal. There the twin, which copies that bimodal margin,
inherits the split and reports the same count, so the test stays quiet
on a structure that is in fact real. This is a feature at the boundary
rather than a defect: a variable that is plainly two-humped needs no
null to reveal its groups. Inspect the margins for pronounced
multimodality before relying on the count test, and tie-break granular
Likert-type scales before formal unimodality tests.

## Reference

The method, its positive controls, and its false-positive calibration
are described in the accompanying paper, *Types Without Taxa: A
Covariance-Matched-Null Multiverse Test of Categorical versus Continuous
Personality Structure* (Meng, 2026; preregistration:
<https://doi.org/10.17605/OSF.IO/2EKCG>).
