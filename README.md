# Hamiltonian Monte Carlo Parameter Estimation for Location Scale Regression Models

## Introduction

Welcome to our collaborative project! In this repository, Sönke Hänel, Jan Parlesak, and I present an R package that enables parameter estimation for location scale regression models using Hamiltonian Monte Carlo (HMC), which was pasrt of the "Advanced Statistical Progarmming" course at the university of Göttingen. The original repository can be found under: https://gitlab.gwdg.de/paul.jarschke/hamiltonian-mcmc

## Key Features

- **Hamiltonian Monte Carlo (HMC)**: Our R package leverages the power of HMC to estimate parameters for location scale regression models.

- **Dual Averaging Algorithm**: With our package, you don't need to manually fine-tune hyperparameters. Our Dual Averaging algorithm does the job for you, ensuring convergence and independent sampling from the target distributions.

- **Methods** We included methods for summarizing model output and plotting target distributions, trace plots and more.

- **Documentation**: To help you get started, we've included comprehensive documentation within the package, that can be looked up using R Studio.

- **Simulation File**: The repository includes a simulation file to compare our estimates to the estimates of the lmls package.

## Acknowledgments

We extend our appreciation to Hannes Riebl for the foundational lmls package that made our work possible.
