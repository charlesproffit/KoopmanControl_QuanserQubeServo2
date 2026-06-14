# Data-driven Control Methods for an Inverted Pendulum using Koopman Operator

This repository contains the MATLAB implementation developed for the semester project:

**"Data-driven Control Methods for an Inverted Pendulum using Koopman Operator"**

## Overview

The goal of this project is to control the Quanser QUBE-Servo 2, a rotary pendulum system, using data-driven methods based on the Koopman operator framework.

Instead of relying on a first-principles nonlinear model, the approach uses data collected from the hardware to:

1. Identify the system dynamics using Extended Dynamic Mode Decomposition (EDMD).
2. Construct a lifted bilinear representation of the nonlinear system.
3. Perform data-driven feedback linearization.
4. Design controllers on the resulting linearized model.

The objective is to stabilize the pendulum around its downward equilibrium while regulating the rotary arm angle to zero.

## Requirements

* MATLAB R2024a or later
* MOSEK

## Repository Structure

```text
.
├── README.md
├── LabView/                    % LabVIEW projects for hardware experiments
│
└── Matlab/
    ├── data/                   % Experimental data (.tdms files)
    ├── datadriven/             % Data-driven robust control toolbox
    ├── figures/                % Generated figures and results
    ├── src/
    │   ├── Controllers/
    │   │   ├── FLController.m
    │   │   ├── RobustController.m
    │   │   └── SFController.m
    │   │
    │   ├── Models/
    │   │   ├── DDFLModel.m
    │   │   ├── EDMDModel.m
    │   │   ├── MBFLModel.m
    │   │   └── NLModel.m
    │   │
    │   ├── Helpers/
    │   │   ├── compare_EDMD.m
    │   │   ├── compare_FL.m
    │   │   ├── plots_EDMD.m
    │   │   ├── QuanserQubeDynamicsCompactForm.m
    │   │   └── structure_data.m
    │   │
    │   └── Simulators/
    │       ├── OLSimulator.m
    │       └── CLSimulator.m
    │
    └── main.m
```

## Usage

1. Open MATLAB and navigate to the `Matlab` folder.
3. Run the main script:

```matlab
main
```

The script loads the experimental data, identifies the Koopman-based model, performs feedback linearization, and evaluates the different control strategies.

## Main Components

### Models

* **NLModel**: Nonlinear model of the Quanser QUBE system.
* **EDMDModel**: Koopman-based system identification using EDMD.
* **DDFLModel**: Data-driven feedback linearization model.
* **MBFLModel**: Model-based feedback linearization model.

### Controllers

* **SFController**: State-feedback controller.
* **FLController**: Feedback-linearization controller.
* **RobustController**: Robust controller.

### Simulators

* **OLSimulator**: Open-loop simulations.
* **CLSimulator**: Closed-loop simulations.

## Usage

The project workflow consists of defining the system dynamics, selecting lifting functions, collecting data, identifying Koopman-based models, and designing controllers.

### 1. Setup

Open MATLAB, navigate to the `Matlab` folder, and add all project directories to the path:

```matlab
clc;
clear;
close all;

addpath(genpath('src'));
addpath(genpath('data'));
addpath(genpath('datadriven'));

rng(43);
```

### 2. Define System Parameters

Specify the simulation and identification parameters:

```matlab
Ts = 0.008;                 % Sampling time
nstates = 4;                % Number of states
ninputs = 1;                % Number of inputs

ntrajs_training = 160;      % Training trajectories
ntrajs_testing = 40;        % Testing trajectories
nsteps = 500;               % Samples per trajectory
numax = 7;                  % Maximum input amplitude
```

The state vector is defined as:

```text
q(1) : Rotary arm angle
q(2) : Pendulum angle
q(3) : Rotary arm angular velocity
q(4) : Pendulum angular velocity
```

### 3. Define the Nonlinear Dynamics

The Quanser QUBE dynamics are expressed in compact form:

```matlab
[Mq,H,Phi,B] = QuanserQubeDynamicsCompactForm();

func_discard = @(x,x0) ...
    abs(x(1)) >= pi/2 || ...
    abs(x(2)-x0(2)) >= 2*pi || ...
    abs(x(2)) >= pi;

Model_NL = NLModel(nstates,ninputs,Ts,Mq,H,Phi,B);
```

The discard function removes trajectories that leave the desired operating region during data collection.

### 4. Define Lifting Functions

A lifting function maps the nonlinear state into a higher-dimensional observable space used by EDMD.

```matlab
K_LQR = Model_NL.designLQR( ...
    [0;0;0;0], ...
    0, ...
    diag([0.41 0 0 0]), ...
    0.016);

func_lifting = @(q)[
    1;
    q;
    Model_NL.f_nr(q) - Model_NL.g_nr(q)*K_LQR*q;
    Model_NL.g_nr(q);
];
```

A linear lifting can also be defined:

```matlab
func_linear = @(q)(q);
```

### 5. Generate Initial Conditions and Inputs

Initial conditions are sampled from different pendulum angle regions:

```matlab
func_initialStates = @() [
    0;
    -pi + angle_region * (randi([2,nregions-1]) + rand() - 1);
    5*pi*(2*rand()-1);
    5*pi*(2*rand()-1)
];
```

Random persistently exciting inputs are generated using:

```matlab
func_inputs = @() ...
    2*numax*rand(ninputs,nsteps) - numax;
```

### 6. Collect Data

An open-loop simulator is used to generate trajectories:

```matlab
Sim_NL = OLSimulator(Model_NL,func_discard);
```

Synthetic datasets can be generated directly from the model:

```matlab
dataset = Sim_NL.generateDataset( ...
    ntrajs_training,...
    ntrajs_testing,...
    nsteps,...
    func_initialStates,...
    func_inputs,...
    'euler');
```

Alternatively, experimental trajectories collected from the Quanser hardware can be loaded from TDMS files:

```matlab
[dataset,nsteps,ntrajs_training,ntrajs_testing] = ...
    structure_data( ...
        'data\data.tdms', ...
        nstates, ...
        ninputs, ...
        0.8, ...
        50, ...
        true);
```

### 7. Run the Main Pipeline

After defining the dynamics, lifting functions, and dataset, execute:

```matlab
main
```

The script performs:

1. EDMD-based system identification.
2. Construction of Koopman lifted models.
3. Data-driven feedback linearization.
4. Controller synthesis.
5. Open-loop and closed-loop validation.
6. Performance comparison and visualization.

## Author

Charles Proffit
