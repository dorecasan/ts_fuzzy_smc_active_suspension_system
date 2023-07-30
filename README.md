# Adaptive Sliding Mode Takagi-Sugeno Fuzzy Controller for Half-Car Active Suspension System

This repository contains my code implementation of an adaptive sliding mode Takagi-Sugeno (TS) fuzzy controller for a half-car active suspension system. 

## Introduction

The goal of this project is to provide an efficient and effective control strategy for the active suspension system, which aims to improve the ride comfort and handling performance of the vehicle. The controller is designed based on the Takagi-Sugeno fuzzy inference system and incorporates an adaptive sliding mode control approach.


## Instructions for Running the Code

To run the code and simulate the active suspension system, follow these steps:

1. Adjust the global parameters defined in the `globalParams.m` file according to your requirements.
2. Run the `check_lmi.m` script in MATLAB to obtain the control gains for the TS model. This script will solve the LMI optimization problem and output the control gains.
3. Open the `vssSimulation.slx` Simulink model file.
4. Click on the "Run" button in Simulink to start the simulation.

## Dependencies

- MATLAB
- Simulink
