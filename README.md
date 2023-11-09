# Hamiltonian Monte Carlo Parameter Estimation for Location Scale Regression Models

## Introduction

Welcome to our collaborative project! In this repository, Sönke Hänel, Jan Parlesak, and Paul Jarschke present an R package that enables parameter estimation for location scale regression models using Hamiltonian Monte Carlo (HMC). We've also implemented a Dual Averaging algorithm that automates the tuning of hyperparameters (step size and step length of the leapfrog integrator) until convergence criteria are met. This results in efficient and accurate sampling, followed by the generation of independent samples from the target distributions.

## Key Features

- **Hamiltonian Monte Carlo (HMC)**: Our R package leverages the power of HMC to estimate parameters for location scale regression models.

- **Dual Averaging Algorithm**: With our package, you don't need to manually fine-tune hyperparameters. Our Dual Averaging algorithm does the job for you, ensuring that the sampling process is optimized.

- **Independent Samples**: Once the hyperparameters are tuned and convergence is achieved, our package facilitates the generation of independent samples from the target distributions. This is essential for reliable and statistically valid results.

- **Documentation**: To help you get started, we've included comprehensive documentation within the. The documentation covers installation, usage, and customization of our package, providing clear and concise explanations.

- **Simulation File**: The repository includes a simulation file to compare our estimates to the estimates of the lmls package.

## Acknowledgments

We extend our appreciation to Hannes Riebl for the foundational lmls package that made our work possible.
