# Forecasting tipping times from data

## 🚀 Description 

In dynamical systems and bifurcation theory, one topic of research that has been continuosly developed in the past two decades has been the investigation of critical transitions, often called _tipping points_ and their prediction through early-warning signals.
Despite the promises made by early works regarding the predictive power of these methods for the forecast and prevention of catastrophic events in ecological and climate systems, the mathematical and scientific communities largely worked independenlty in the progress of the theory.
Mathematicians have refined more and more the constraints and modelling of the theoretical aspects underpinning tipping events and their prediction, while scientists and sought to derive statistical indicators to apply to very specific simulated systems of natural collapse.
The gap that exists, at this point, between the mathematical theory and the practical applications is substantial and very little efforts have tried to address it in full generality.

### ✏️📋 Outline

__Collaborators__: [Lauren Smith](https://profiles.auckland.ac.nz/lauren-smith) (U. of Auckland, New Zealand) and [Graham Donovan](https://profiles.auckland.ac.nz/g-donovan) (U. of Auckland, New Zealand).

We address the above issues head-on by deriving a practical, end-to-end pipeline for scientists to forecast the time at which their systems are supposed to tip in the near or far future.
This pipeline is an inference-extrapolation meshup combining several different results from the literature.
On top of this, we derive optimal conditions for the collection and analysis of timeseries data that scientists can use to ensure robust detection of the tipping time using the proposed methodology.

### 🔍📜 How to cite

_Work in progress_

## 📦 Structure of the repo

This repository is organised as follows

```bash
root/
├── doc/              # LaTeX source code of the paper 
└── sim/              # Collection of simulations in the paper 
```

### ⚙️🛠️ Organisation of the experiments

In `sim` you'll find many subdirs, each associated with an important experiment and/or a figure in the paper. The code in each subdir is organized as follows

```bash
experiment/
├── data/            # If present it holds data used in the scripts 
├── src/             
│   ├── plotting     # Customise the layouts of the figures 
│   ├── processing   # Postprocess the results of a simulation 
│   └── setup        # Defines and sets-up the simulation 
└── main             # Runs the experiment and plots the figures
```

### 💡🔬 What does each experiment do?

What follows is a brief description of those experiments so that you can gather their functionalities without necessarly interpreting the code in it:

- [01_tipping_time_distribution](./sim/01_tipping_time_distribution/): simulations of three different tipping regimes for the [saddle-node normal form](https://en.wikipedia.org/wiki/Saddle-node_bifurcation);
- [02_admissible_time_horizon](./sim/02_admissible_time_horizon/): computation of the return rate of the drifting steady equilibrium and the [pullback attractor](https://en.wikipedia.org/wiki/Pullback_attractor) of the nonautonomous system;
- [03_bias_estimators](./sim/03_bias_estimators/): accuracy and [uncertantiy quantification](https://en.wikipedia.org/wiki/Uncertainty_quantification#Bias_correction_only) of three different statistical estimators of the return rate from the sample;
- [04_observation_window](./sim/04_observation_window/): plot of the bounded observation window and the effects on the samples;
- [05_dynamic_saddle_node](./sim/05_dynamic_saddle_node/): application of the forecasting pipeline to a simulated ensemble of sample paths of a dynamic saddle-node normal form;
- [06_amoc_collapse](./sim/06_amoc_collapse/): application of the forecasting pipeline to the freshwater flux timeseries of a [Community Earth System Model](https://www.cesm.ucar.edu/) of the multistable [Atlantic Meridional Overturning Circulation](https://en.wikipedia.org/wiki/Multiple_equilibria_in_the_Atlantic_meridional_overturning_circulation) collapse.
