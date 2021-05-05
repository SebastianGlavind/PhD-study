# PhD study

This repository contains a set of toolboxes and code examples developed during my PhD study at Aalborg University. Its main purpose is to support my PhD thesis, but it also contains some tutorials developed for teaching on related topics. The links below display 
the tutorials via [nbviewer](https://nbviewer.jupyter.org/) to ensure a proper rendering of the formulas.

## Bayesian networks

The following toolboxes and tutorials are implemented in R;

### Toolboxes

The toolboxes are among others used in Glavind and Faber (2018), and Glavind and Faber (2020).

- [Structure learning and dynamic discretization toolbox - `sLearningAndDiscretizationTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R). The toolbox is a wrapper for the [`bnlearn`](https://www.bnlearn.com/) package, which implements the underlaying score-based routines for structure learning.

- [Parameter learning toolbox - `pLearningTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). The toolbox is a wrapper for the [`bnlearn`](https://www.bnlearn.com/) package.

### Structure learning

- [Bayesian networks - structure learning and automated discretization form complete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_sLearn_fullyObs.ipynb). This tutorial demonstrates how to learn the graph structure and optimal discretization policy of a Bayesian network (BN) representation from complete / fully observed data using my toolbox [`sLearningAndDiscretizationTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R). 

- [Bayesian networks - structure learning and automated discretization from incomplete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_sLearn_partiallyObs.ipynb). This tutorial demonstrates how to learn the graph structure and optimal discretization policy of a Bayesian network (BN) representation from incomplete / partially observed data using my toolbox [`sLearningAndDiscretizationTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R).

### Parameter learning

- [Bayesian networks - parameter learning from complete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_pLearn_fullyObs.ipynb). This tutorial demonstrates how to learn the parameters of a Bayesian network (BN) representation from complete / fully observed data using my toolbox [`pLearningTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). 

- [Bayesian networks - EM for parameter learning from incomplete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_pLearn_EM_partiallyObs.ipynb). This tutorial demonstrates how to learn the parameters of a Bayesian network (BN) representation from incomplete / partially observed data using the expectation-maximization (EM) algorithm. An implementation of the EM algorithm is found in my toolbox [`pLearningTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). 

- [Bayesian networks - Gibbs sampling for parameter learning from incomplete data](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_pLearn_Gibbs_partiallyObs.ipynb). This tutorial demonstrates how to learn the parameters of a Bayesian network (BN) representation from incomplete / partially observed data using Gibbs sampling. The implementation makes use of some functionalities from my toolbox [`pLearningTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R). 

### Inference

- [Bayesian networks - inference for discrete Bayesian networks](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/BNs_inference.ipynb). This tutorial demonstrates how to make inferences using general bn.fit objects; these objects may have be learned using the [`bnlearn`](https://www.bnlearn.com/) package alone or in combination with my toolboxes [`sLearningAndDiscretizationTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/sLearningAndDiscretizationTools.R) and [ `pLearningTools`](https://github.com/SebastianGlavind/PhD-study/blob/master/Bayesian-networks/Toolboxes/pLearningTools.R), which are wrappers for the `bnlearn` package. For this tutorial we will use the inference functionalities of the `bnlearn` package, as well as the [`gRain`](http://people.math.aau.dk/~sorenh/software/gR/) package to make maximum a-posteriori inferences, as well as posterior inferences that account for parameter uncertainties. 

## Linear regression

The following tutorials are implemented in Python;

- [Linear regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/LinReg.ipynb). This tutorial introduces linear regression; first, from a maximum likelihood estimation (MLE) perspective, and second, from a Bayesian perspective. In both cases, the tutorial implements a selection of different learning algorithms. 

- [Linear regression - assumptions and interpretations](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/LinReg_assumptionsEtc.ipynb). This notebook considers and assesses the underlaying assumptions of linear regression in detail and discusses the interpretation of these models. 

- [Bayesian linear regression with Stan](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/exStan_BayesLinReg.ipynb). This tutorial shows how to implement Bayesian linear regression models using the probabilistic programming language [Stan](https://mc-stan.org/). 

- [EM for Bayesian linear regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Linear-regression/exEM_BayesLinReg.ipynb). This tutorial considers how the expectation maximization (EM) algorithm may be used to learn a parameter setting for a Bayesian linear regression model.

## Bayesian hierarchical models

The following tutorials are implemented in Python;

- [Bayesian hierarchical models with Stan](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Hierarchical-models/HierModel_OMAE2020.ipynb). This tutorial introduces how to implement Bayesian hierarchical regression models using the probabilistic programming language [Stan](https://mc-stan.org/) by studying the fatigue data set in Glavind et al. (2020). Moreover, the concept of Bayesian model averaging is introduced as a means for making inferences for new out-of-sample fatigue sensitive details. 

## Gaussian processes

The following tutorials are implemented in Python;

- [An intuitive introduction to Gaussian processes - Gaussian processes as the generalization of Gaussian probability distributions](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-processes/GPintro_intuitive.ipynb). In this tutorial, we will consider how to generalize the Gaussian probability distribution to infinite dimensions to approach the Gaussian process. The tutorial is intuitive in the sense that we will consider graphics instead of proofs to motivate the transition from the finite dimensional Gaussian probability distribution to the Gaussian process.

- [Gaussian process regression](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-processes/GPR.ipynb). This tutorial introduces Gaussian process regression; first, in a single-output setting, and second, in a multi-output setting. For the single-output case, the tutorial implements a selection of different learning algorithms, and some of the capabilities of the open source software package [`GPy`](https://sheffieldml.github.io/GPy/) are demonstrated for both cases.

- Bayesian optimization using a Gaussian process prior for hyperparameter tuning, see the tutorials described below on [*Gradient boosting regression using XGBoost*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMReg_BostonHousing.ipynb) and [*Gradient boosting classification using XGBoost*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMClas_Wine.ipynb), respectively.

## Neural networks

The following tutorials are implemented in Python;

- [Neural network regression using keras and tensorflow](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Neural-networks/NNReg_BostonHousing.ipynb). This tutorial introduces neural network regression with keras and tensorflow by considering the Boston housing data set; first, in a single-output setting; and second, in a multi-output setting. Finally, the tutorial considers hyperparameters tuning in general models using random search cross-validation. 

- [Neural network classification using keras and tensorflow](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Neural-networks/NNClas_Wine.ipynb). This tutorial introduces neural network classification with keras and tensorflow by considering the Wine recognition data set. The tutorial first study how a neural network is implemented for classification tasks and then considers how to tune hyperparameters in general models using random search cross-validation. 

## Tree-based learners

The following tutorials are implemented in Python;

- [Gradient boosting regression using XGBoost](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMReg_BostonHousing.ipynb). This tutorial introduces gradient boosting regression with [`XGBoost`](https://xgboost.readthedocs.io/en/latest/) by considering the Boston housing data set. The tutorial first study how gradient boosting is implemented in a single-output setting as well as the effect of different data pre-processing steps. Then, it is shown how gradient boosting may be extended to a multi-output setting. Finally, the tutorial considers hyperparameters tuning in general models using Bayesian optimization with a Gaussian process prior based [`GPyOpt`](https://github.com/SheffieldML/GPyOpt). 
 
- [Gradient boosting classification using XGBoost](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMClas_Wine.ipynb). This tutorial introduces gradient boosting classification with [`XGBoost`](https://xgboost.readthedocs.io/en/latest/) by considering the Wine recognition data set. The tutorial first study how gradient boosting is implemented as well as the effect of different data pre-processing steps. Then, the tutorial considers hyperparameters tuning in general models using Bayesian optimization with a Gaussian process prior based [`GPyOpt`](https://github.com/SheffieldML/GPyOpt). Finally, the tutorial elaborates on the feature importance functionalities of `XGBoost`.

## Gaussian mixture models

The following tutorials are implemented in Python;

- [EM for Gaussian mixtures](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gaussian-mixtures/exEM_GMMs.ipynb). This tutorial considers how Gaussian mixture models may be used for cluster analysis; it implements the expectation maximization (EM) learning algorithm, and introduces the evidence lower bound, as well as the Bayesian information criterion (BIC) and the integrated complete-data likelihood (ICL), for model selection.

## Sensitivity analysis and feature selection

The following tutorials are implemented in Python;

- [Variance-based sensitivity analysis for independent inputs](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Sensitivity-analysis/SA_varianceBased_independentInputs.ipynb). This tutorial implements a set of methods, which are applicable when the inputs are independent. First, a surrogate-based method is considered that decomposes the variance based on linear regression considerations. Second, two simulation-based methods are introduced; the first method performs conditional sampling by binning the input space, and the second method performs efficient conditional sampling.  

- [Variance-based sensitivity analysis for correlated inputs](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Sensitivity-analysis/SA_varianceBased_correlatedInputs.ipynb). This tutorial implements a set of methods, which are applicable when the inputs are correlated. First, two surrogate-based methods are considered; the first method decomposed the variance based on (linear) regression considerations, and the second method decomposes the variance based on a polynomial chaos expansion. Second, two simulation-based methods are introduced; the first method performs conditional sampling by binning the input space, and the second method performs conditional sampling for randomly sampled input realizations.  

## Hyperparameter tuning, model selection and automated machine learning (AutoML)

- Random search for hyperparameter tuning, see the tutorials described above on [*Neural network regression using keras and tensorflow*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Neural-networks/NNReg_BostonHousing.ipynb) and [*Neural network classification using keras and tensorflow*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Neural-networks/NNClas_Wine.ipynb), respectively.

- Bayesian optimization using a Gaussian process prior for hyperparameter tuning, see the tutorials described above on [*Gradient boosting regression using XGBoost*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMReg_BostonHousing.ipynb) and [*Gradient boosting classification using XGBoost*](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Gradient-boosting/GBMClas_Wine.ipynb), respectively.

## Algorithms for optimization

The following tutorials are implemented in Python;

- [Deterministic algorithms for unconstrained, continuous-valued optimization](https://nbviewer.jupyter.org/github/SebastianGlavind/PhD-study/blob/master/Optimization/Optimization_con_uncon.ipynb). This tutorial considers set of local derivative-based optimization algorithms for unconstrained, continuous-valued optimization. The algorithms covered are first-order methods, i.e., gradient decent and its variations (e.g., conjugate gradient decent and Adam), and second-order methods, i.e. Newton's method and quasi-Newton methods (DFP and BFGS).

## References
***
Sebastian T. Glavind and Michael H. Faber, “A framework for offshore load environment modeling”, in proceedings of the ASME 2018 37th International Conference on Ocean, Offshore and Arctic Engineering (OMAE2018), OMAE2018-77674, 2018.

Sebastian T. Glavind and Michael H. Faber, “A framework for offshore load environment modeling”, Journal of Offshore Mechanics and Arctic Engineering, vol. 142, no. 2,
pp. 021702, OMAE-19-1059, 2020.

Sebastian T. Glavind, Henning Brüske and Michael H. Faber, “On normalized fatigue crack growth modeling”, in proceedings of the ASME 2020 39th International Conference on Ocean, Offshore and Arctic Engineering (OMAE2020), OMAE2020-18613, 2020.
***
