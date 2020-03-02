# LRMoE

- [Overview](#overview)
- [Model Description](#model-description)
- [Package Installation and Usage](#package-installation-and-usage)
- [Supported Distributions](#supported-distributions)
- [Issue Report and Suggestions](#issue-report-and-suggestions)
- [Planned Development](#planned-development)

# Overview 

**LRMoE** is an R package tailor-made for actuarial applications which allows actuarial researchers and practitioners to model and analyze insurance loss frequencies and severities using the Logit-weighted Reduced Mixture-of-Experts (LRMoE) model. The flexibility of LRMoE models is theoretically justified in [Fung et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0167668719303956), and an application of LRMoE for modelling correlated insurance claim frequencies is in [Fung et al. (2019)](https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/class-of-mixture-of-experts-models-for-general-insurance-application-to-correlated-claim-frequencies/E9FCCAD03E68C3908008448B806BAF8E).

The package **LRMoE** offers several new distinctive features which are motivated by various actuarial applications and mostly cannot be achieved using existing packages for mixture models. Key features include:
* A wider coverage on frequency and severity distributions and their zero inflation;
* The flexibility to vary classes of distributions across components;
* Parameter estimation under data censoring and truncation;
* A collection of insurance rate making and reserving functions; and
* Model selection and visualization tools.

# Model Description

Let ![](https://latex.codecogs.com/svg.latex?x_{i}) denote the ![](https://latex.codecogs.com/svg.latex?(P+1))-dimensional covariate vector (age, sex, policy region, etc.) of policyholder ![](https://latex.codecogs.com/svg.latex?i). Based on the covariates, the policyholder is classified into one of ![](https://latex.codecogs.com/svg.latex?g) latent risk classes by a logit *gating function* such that the probability of class ![](https://latex.codecogs.com/svg.latex?j) is

![](https://latex.codecogs.com/svg.latex?\pi_{j}(x_{i};&space;\alpha)&space;=&space;\frac{\exp(\alpha_j^T&space;x_i)}{\sum_{j^\prime&space;=&space;i}^{g}\exp(\alpha_{j^\prime}^T&space;x_i)}.)

Let ![](https://latex.codecogs.com/svg.latex?y_{i}) denote the ![](https://latex.codecogs.com/svg.latex?D)-dimensional response vector (claim number, claim amount, etc.) of policyholder ![](https://latex.codecogs.com/svg.latex?i). Conditional on being assigned to latent risk class ![](https://latex.codecogs.com/svg.latex?j), each dimension of ![](https://latex.codecogs.com/svg.latex?y_{i}) are conditionally independent with the ![](https://latex.codecogs.com/svg.latex?d)-th marginal density given by the *expert function*

![](https://latex.codecogs.com/svg.latex?g_{jd}&space;=&space;\delta_{jd}I_{\{y_{id}&space;=&space;0\}}&space;&plus;&space;(1-\delta_{jd})f_{jd}(y_{id};&space;\psi_{jd}))

where ![](https://latex.codecogs.com/svg.latex?\delta_{jd}) represents a zero-inflation probability mass, ![](https://latex.codecogs.com/svg.latex?I_{y_{jd}=0}) is the indicator function, and ![](https://latex.codecogs.com/svg.latex?f_{jd}) is a commonly-used parametric distribution for actuarial loss modelling with parameter ![](https://latex.codecogs.com/svg.latex?\psi_{jd}). For policyholder ![](https://latex.codecogs.com/svg.latex?i), the density function of ![](https://latex.codecogs.com/svg.latex?y_{i}) is therefore

![](https://latex.codecogs.com/svg.latex?f(y_{i};&space;x_{i},&space;\alpha,&space;\delta,&space;\Psi)&space;=&space;\sum_{j=1}^{g}&space;\pi_{j}(x_{i};&space;\alpha)\prod_{d=1}^{D}g_{jd}(y_{jd};&space;\delta_{jd},&space;\psi_{jd}).)

The parameters ![](https://latex.codecogs.com/svg.latex?(\alpha,&space;\delta,&space;\Psi)) can be estimated using the Expectation-Conditional-Maximization algorithm, which is implemented in this package.

# Package Installation and Usage

While we prepare for submitting the package to Comprehensive R Archive Network ([CRAN](https://cran.r-project.org/)), the package is available in this github repository. The package manual documenting all functions in the package can also be found in [the current repository](https://github.com/sparktseung/LRMoE/blob/master/LRMoE_0.1.0.pdf). 

In R, the **LRMoE** package can be downloaded and installed by running the following code.

```R
library(devtools)
install_github("sparktseung/LRMoE")
```

For a detailed demonstration of using the package, we have set up a a separate repository LRMoE-Paper-Demo, which also accompanies a paper submitted to the Annals of Actuarial Science (as of Feb 29, 2020). In that paper, we have provided two illustrative examples (one on a simulated dataset and another on a real dataset), in order to demonstrate the basic procedures of model fitting, selection, visualization and application.

# Supported Distributions



# Issue Report and Suggestions


# Planned Development
