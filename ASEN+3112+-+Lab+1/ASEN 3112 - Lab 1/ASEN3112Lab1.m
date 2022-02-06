%% ASEN 3112 - Lab 1 - Torsion
% 
% Authors:
%   1) Ross White
%   2) Jacob Wilson
%   3) Tristan Workman
%   4) Alicia Wu
%
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specimen Parameters
% Parameters of both specimens
diameterExterior = 0.75; % [in]
radiusExterior = diameterExterior/2; % [in]
thickness = 1/16; % [in]
length = 10; % [in]
perimeter = 2*pi*radiusExterior; % [in]
crossSectionArea = perimeter*thickness; % [in^2]
enclosedArea = pi*(radiusExterior-(thickness/2))^2; % [in]
shearModulus = 3.75*10^6; % [psi]

% Parameters of Closed tube (CTW)
torqueNaughtCTW = 0; % [lbs-in]
torqueMaxCTW = 400; % [lbs-in]
shearStressMaxCTW = 8620; % [psi]

% Parameters of Slitted tube (OTW)
torqueNaughtOTW = 0; % [lbs-in]
torqueMaxOTW = 20; % [lbs-in]
shearStressMaxOTW = 7800; % [psi]

%% Loading Experimental Data
% Closed tube Experimental Data
dataCTW = table2array(readtable('Torsion400inlbf_142.txt', 'VariableNamingRule', 'preserve'));
timeCTW = dataCTW(:,1);  % [sec] time data for CTW specimen
angleCTW = dataCTW(:,2).*(pi/180);  % [rad] twist angle on CTW specimen
shearStrainCTW = dataCTW(:,3).*(pi/180); % [rad] shear strain on CTW specimen
torqueCTW = dataCTW(:,4); % [lbf*in] torsional torque on CTW specimen
axialStrainCTW = dataCTW(:,5); % [] axial strain on CTW specimen

% Slitted tube Experimental Data
dataOTW = table2array(readtable('Torsion20inlbf_143.txt', 'VariableNamingRule', 'preserve'));
timeOTW = dataOTW(:,1);  % [sec] time data for OTW specimen
angleOTW = dataOTW(:,2);  % [deg] twist angle on OTW specimen
shearStrainOTW = dataOTW(:,3); % [deg] shear strain on OTW specimen
torqueOTW = dataOTW(:,4); % [lbf*in] torsional torque on OTW specimen
axialStrainOTW = dataOTW(:,5); % [] axial strain on OTW specimen

%% Analysis of the Closed Thin Wall Specimen
% • Plot the torque vs.  shear strain provided by the extensometer, as well 
%   as the torque vs.  shear strain calculated using the total rotation 
%   angle imposed by the testing machine.

% Calculation of Torque and Shear Strain using Twist Angle
J_CTW = 4*crossSectionArea^2/(perimeter/thickness); % [in^4]
GJ_CTW = shearModulus*J_CTW; % [lbs-in^2]
torqueCalculatedCTW = GJ_CTW.*angleCTW; % [lbs-in^2]
shearStressCalculatedCTW = torqueCalculatedCTW./(2.*thickness.*enclosedArea); % [psi]
shearStrainCalculatedCTW = shearStressCalculatedCTW./shearModulus; % [psi]
shearStrainCalculatedCTW = shearStrainCalculatedCTW - shearStrainCalculatedCTW(1);

% Plotting
figure
subplot(1,2,1)
hold on
title('Torque vs. Shear Strain - Closed Specimen')
xlabel('Shear Strain [radians]');
ylabel('Applied Torque [lbf*in]');
plot(shearStrainCTW,torqueCTW);
plot(shearStrainCalculatedCTW,torqueCTW);
legend('Shear Strain Experimental Data','Shear Strain Calculated using Thin Wall Theory');
hold off

% • Use least squares fitting to calculate the torsional rigidity, GJ, for 
%   the two ways to obtain the shear strain.  Provide the value of the 
%   associated uncertainty.

% • Compare the value of GJ obtained through the experiments with the 
%   theoretical predictions obtained using exact theory and thin wall theory.  
%   Discuss the differences.

%% Analysis of the Open Thin Wall Specimen
% • Plot the torque vs.  shear strain provided by the extensometer, as well 
%   as the torque vs.  shear strain calculated using the total rotation 
%   angle imposed by the testing machine.

% Calculation of Torque and Shear Strain using Twist Angle
% Assume alpha and beta are both 1/3 because p/t >> 5
% This means that J_alpha = J_beta = J_OTW
J_OTW = (1/3).*perimeter.*thickness.^3; % [in^4]
GJ_OTW = shearModulus*J_OTW; % [lbs-in^2]
torqueCalculatedOTW = GJ_OTW.*angleOTW; % [lbs-in^2]
shearStressCalculatedOTW = torqueCalculatedOTW.*thickness./J_OTW; % [psi]
shearStrainCalculatedOTW = shearStressCalculatedOTW./shearModulus; % [psi]
shearStrainCalculatedOTW = shearStrainCalculatedOTW - shearStrainCalculatedOTW(1);

% Plotting
subplot(1,2,2)
hold on
title('Torque vs. Shear Strain - Open Specimen')
xlabel('Shear Strain [radians]');
ylabel('Applied Torque [lbf*in]');
plot(shearStrainOTW,torqueOTW);
plot(shearStrainCalculatedOTW,torqueOTW);
legend('Shear Strain Experimental Data','Shear Strain Calculated using Thin Wall Theory');
hold off

% • Use least squares fitting to calculate the torsional rigidity, GJ, for 
%   the two ways to obtain the shear strain.  Provide the value of the 
%   associated uncertainty.

% • Compare the value of GJ obtained through the experiments with the 
%   theoretical predictions obtained using exact theory and thin wall theory.  
%   Discuss the differences.


