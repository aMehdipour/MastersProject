clear; close all; clc

addpath('/home/amehdipour/repos/MastersProject/Sim/Data');
addpath('/home/amehdipour/repos/MastersProject/Sim/QCAT/qcat');
addpath('/home/amehdipour/repos/MastersProject/Sim/src');
addpath('/home/amehdipour/repos/MastersProject/Sim/');

%% Load Data
load CFD; load std_atmos; load att_profile
set(0,'DefaultFigureWindowStyle','docked')

rng(69); % Set the RNG seed for reproducability

CFD.aeroScalingFactor = 1; % Scaling factor for aero variables to add scaling error

%% Initialize
[constants, state, sim] = initialize();


