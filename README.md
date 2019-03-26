# Pvalue_distribution

Brief description: Fit a mixture distribution to a dataset of pvalues with R

# Distribution mixture for a set of p-values

Given a set of p-value, this is an example about how to fit a mixture distribution to this data

Specifically, next mixture distributions were fitted:

1. Uniform distribution + 1 beta distribution
2. Uniform distribution + 2 beta distributions
3. Uniform distribution + 2 triangular distributions
4. Uniform distribution + 2 exponential distributions

Using a sample of p-values, the aim of this example is to estimate the proportion of p-values comming from the null hypothesis (H<sub>0</sub>) of several tests (regardless of your type) with distribution U(0,1)

The methodology is quite similar to that one explained in this article:

1. Pounds, S. & Morris, S. W. [Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values](https://academic.oup.com/bioinformatics/article/19/10/1236/184434). Bioinformatics 19, 1236–42 (2003)

# Files

- __*pvalues.txt*__: Dataset with pvalues coming from different statistical tests.
- __*pvalue_distribution_mixture.R*__: Main script to determine the proportion of p-values coming from the null hypothesis (H<sub>0</sub>)
- __*pvalue_distribution_mixture_functions.R*__: functions to run the main script. Mainly, the implementation of the likelihood functions.
