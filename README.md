# PhD study

This repository contains some toolboxes and code examples developed during my PhD study at Aalborg University. Its main purpose is to support my PhD thesis, but it also contains some tutorials developed for teaching on related topics. The links below display 
the tutorials via [nbviewer](https://nbviewer.jupyter.org/) to ensure a proper rendering of formulas.

## Bayesian networks

- [Bayesian networks - structure learning from complete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/sLearn_fullyObs.ipynb). This tutorial demonstrates how to learn the graph structure and optimal discretization policy of a Bayesian network (BN) representation from complete / fully observed data, using my toolbox [sLearningAndDiscretizationTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R). The toolbox is a wrapper for the [bnlearn](https://www.bnlearn.com/) package, which implements the underlaying score-based routines for structure learning.

## Linear regression

- [Linear regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/LinearRegression.ipynb). This tutorial introduces linear regression; first, from a maximum likelihood estimation (MLE) perspective, and second, from a Bayesian perspective. In both cases, the tutorial implements a selection of different learning algorithms. 

- [Bayesian linear regression with Stan](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/exStan_BayesianLinearRegression.ipynb). This tutorial show how to implement Bayesian linear regression models using the probabilistic programming language [Stan](https://mc-stan.org/).

## Gaussian processes

- [Gaussian process regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-processes/GaussianProcessRegression.ipynb). This tutorial introduces Gaussian process regression; first, in a single-output setting, and second, in a multi-output setting. For the single-output case, the tutorial implements a selection of different learning algorithms, and some of the capabilities of the open source software [GPy](https://sheffieldml.github.io/GPy/) shown for both cases.

## Neural networks
...

## Gaussian mixture models

- [Model-based cluster analysis](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-mixtures/GaussianMixtures.ipynb). This tutorial considers how Gaussian mixture models may be used for cluster analysis; it implements the expectation maximization (EM) learning algorithm, and introduces the evidence lower bound, as well as the Baysian information criterion (BIC) and the integrated complete-data likelihood (ICL), for model selection.

## Sensitivity analysis
...
