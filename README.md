# PhD study

This repository contains some toolboxes and code examples developed during my PhD study at Aalborg University. Its main purpose is to support my PhD thesis, but it also contains some tutorials developed for teaching on related topics. The links below display 
the tutorials via [nbviewer](https://nbviewer.jupyter.org/) to ensure a proper rendering of formulas.

## Bayesian networks

- [Structure learning and dynamic discretization toolbox - sLearningAndDiscretizationTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R). The toolbox is a wrapper for the [bnlearn](https://www.bnlearn.com/) package, which implements the underlaying score-based routines for structure learning.

- [Parameter learning toolbox - pLearningTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). The toolbox is a wrapper for the [bnlearn](https://www.bnlearn.com/) package.

- [Bayesian networks - structure learning and automated discretization form complete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_sLearn_fullyObs.ipynb). This tutorial demonstrates how to learn the graph structure and optimal discretization policy of a Bayesian network (BN) representation from complete / fully observed data, using my toolbox [sLearningAndDiscretizationTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R). 

- [Bayesian networks - structure learning and automated discretization from incomplete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_sLearn_partiallyObs.ipynb). This tutorial demonstrates how to learn the graph structure and optimal discretization policy of a Bayesian network (BN) representation from incomplete / partially observed data, using my toolbox [sLearningAndDiscretizationTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R).

- [Bayesian networks - parameter learning from complete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_pLearn_fullyObs.ipynb). This tutorial demonstrates how to learn the parameters of a Bayesian network (BN) representation from complete / fully observed data, using my toolbox [pLearningTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). 

- [Bayesian networks - inference for discrete Bayesian networks](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_inference.ipynb). This tutorial demonstrates how to make inferences using general bn.fit objects; these may be learned using the [bnlearn](https://www.bnlearn.com/) package alone, or in combination with my toolboxes [sLearningAndDiscretizationTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R) and [pLearningTools](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R), which are wrappers for the bnlearn package. For this tutorial we will use the inference functionalities of the [bnlearn](https://www.bnlearn.com/) package, as well as the [`gRain`](http://people.math.aau.dk/~sorenh/software/gR/) package to make maximum a-posteriori inferences, as well as posterior inferences that account for parameter uncertainties. 

## Linear regression

- [Linear regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/LinReg.ipynb). This tutorial introduces linear regression; first, from a maximum likelihood estimation (MLE) perspective, and second, from a Bayesian perspective. In both cases, the tutorial implements a selection of different learning algorithms. 

- [Bayesian linear regression with Stan](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/exStan_BayesLinReg.ipynb). This tutorial show how to implement Bayesian linear regression models using the probabilistic programming language [Stan](https://mc-stan.org/).

- [EM for Bayesian linear regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/exEM_BayesLinReg.ipynb). This tutorial considers how the expectation maximization (EM) algorithm may be used to learn a parameter setting for a Bayesian linear regression model.

## Bayesian hierarchical models
...

## Gaussian processes

- [Gaussian process regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-processes/GPR.ipynb). This tutorial introduces Gaussian process regression; first, in a single-output setting, and second, in a multi-output setting. For the single-output case, the tutorial implements a selection of different learning algorithms, and some of the capabilities of the open source software [GPy](https://sheffieldml.github.io/GPy/) are demonstrated for both cases.

## Neural networks
....

## Gaussian mixture models

- [EM for Gaussian mixtures](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-mixtures/exEM_GMMs.ipynb). This tutorial considers how Gaussian mixture models may be used for cluster analysis; it implements the expectation maximization (EM) learning algorithm, and introduces the evidence lower bound, as well as the Bayesian information criterion (BIC) and the integrated complete-data likelihood (ICL), for model selection.

## Sensitivity analysis
...
