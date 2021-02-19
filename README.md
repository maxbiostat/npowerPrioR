# Normalised power prior examples
Tutorial code for sampling from an approximate normalised power posterior using Stan.

This repository holds code and instructions to help you use the methods in [On the normalised power prior](https://arxiv.org/abs/2004.14912) for your own analyses.

You too can use the methods developed here to do your own (normalised) power prior analysis. 
Let's use the Bernoulli experiment in [Neuenschwander et al. (2009)](https://www.ncbi.nlm.nih.gov/pubmed/19735071) as an example.
The steps are:
1. Write your model in Stan, like [this](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/stan/simple_Bernoulli_prior.stan);
2. Then, use [grid_builder.r](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/grid_builder.r) to estimate `c(a0)` at a number of points, for instance as done [here](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/simple_Bernoulli_estimate_c(a0).r);
3. Now you can modify your Stan program to take a dictionary of approximate values for `c(a0)`, like [so](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/stan/simple_Bernoulli_posterior_normalised_approximate.stan);
4. Finally, run your posterior analysis, as exemplified [here](https://github.com/maxbiostat/propriety_power_priors/blob/master/code/simple_Bernoulli_posterior.r);

Any doubts, shoot me a message at `lmax` dot `fgv` at gmail.

Many thanks to Chris Koenig and [Ben Jones](https://www.plymouth.ac.uk/staff/ben-jones) for testing early version of the code.
