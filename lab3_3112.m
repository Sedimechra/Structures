%% ASEN 3112 - Lab 3 - FEM Analysis
% Section 011 - Group 3
% 
% Authors: 
%     1. Luca Bonarrigo
%     2. Tristan Workman
%     3. Abby Durell
%     4. Reina Krumvieda
%     5. Cordelia Kohuth
%     6. Kevin Pipich
%     7. Andrew Sapuppo
%
% Created: 11/17/2021 
% Last edited: 11/28/2021 
%

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Given Constants
L = 12;     % cantilever span S'B' [in]
L_E = 4.5;  % elevator span [in]
L_R = 5;    % rudder span [in]
w = 1;      % width of all members [in]
h = 1/8;    % thickness of fuselage [in]
h_E = 1/4;  % thickness of elevator [in]
h_R = 0.04; % thickness of rudder [in]
E = 10175e3;   % material elastic modulus [psi]
rho = 0.0002505;    % material mass density [lb-sec^2/in^4]
M_T = 1.131*rho;    % tail assembly mass [llb-sec^2/in]
S_T = 0.5655*rho;    % tail assembly first mass moment wrt B' [lb-sec^2]
I_T = 23.124*rho;    % tail assembly second mass moment wrt B' [lb-sec^2-in]

I_zz = w * h^3 / 12; % moment of inertia [in^4]
A = w * h;           % area of fuselage[in^2]
c_K2 = 4*E*I_zz/L^3;
c_M2 = rho*A*L/100800;

M2 = c_M2.*[19272 1458*L 5928 -642*L 0 0;...
            1458*L 172*L^2 642*L -73*L^2 0 0;...
            5928 642*L 38544 0 5928 -642*L;...
            -642*L -73*L^2 0 344*L^2 642*L -73*L^2;...
            0 0 5928 642*L 19272 -1458*L;...
            0 0 -642*L -73*L^2 -1458*L 172*L^2] +...
            [0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 M_T S_T;...
            0 0 0 0 S_T I_T];
K2 = c_K2.*[24 6*L -24 6*L 0 0;...
            6*L 2*L^2 -6*L L^2 0 0;...
            -24 -6*L 48 0 -24 6*L;...
            6*L L^2 0 4*L^2 -6*L L^2;...
            0 0 -24 -6*L 24 -6*L;...
            0 0 6*L L^2 -6*L 2*L^2];
% Delete first 2 columns and rows of both matricies (reduced stiffness/mass
% matrix)
M2(1:2,:) = [];
M2(:,1:2) = [];
K2(1:2,:) = [];
K2(:,1:2) = [];

% syms w
% eqn = det(K2 - w^2*M2) == 0;
% S = vpasolve(eqn, w);
% omega_2 = sqrt(S(find(S>0)));
% freq_2 = omega_2(1:3)./(2*pi)

[eigvec2,eigval2] = eig(K2,M2); % solves for eigenvalue/eigenvectors -> eigenvalues net us omega^2
omega_2 = sqrt(diag(eigval2)); % [rad/s]
freq_2 = omega_2(1:3)/(2*pi)


%% Matrix Definition [4 Element Problem]
c_M4 = rho * A * L / 806400;
c_K4 = 8 * E * I_zz / L^3;

M4 = c_M4.*[77088 2916*L 23712 -1284*L 0 0 0 0 0 0;...
        2916*L 172*L^2 1284*L -73*L^2 0 0 0 0 0 0;...
        23712 1284*L 154176 0 23712 -1284*L 0 0 0 0;...
        -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0 0 0;...
        0 0 23712 1284*L 154176 0 23712 -1284*L 0 0;...
        0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0;...
        0 0 0 0 23712 1284*L 154176 0 23712 -1284*L;...
        0 0 0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2;...
        0 0 0 0 0 0 23712 1284*L 77088 -2916*L;...
        0 0 0 0 0 0 -1284*L -73*L^2 -2916*L 172*L^2] +...
        [0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 M_T S_T;...
        0 0 0 0 0 0 0 0 S_T I_T];
K4 = c_K4.*[96 12*L -96 12*L 0 0 0 0 0 0;...
        12*L 2*L^2 -12*L L^2 0 0 0 0 0 0;...
        -96 -12*L 192 0 -96 12*L 0 0 0 0;...
        12*L L^2 0 4*L^2 -12*L L^2 0 0 0 0;...
        0 0 -96 -12*L 192 0 -96 12*L 0 0;...
        0 0 12*L L^2 0 4*L^2 -12*L L^2 0 0;...
        0 0 0 0 -96 -12*L 192 0 -96 12*L;...
        0 0 0 0 12*L L^2 0 4*L^2 -12*L L^2;...
        0 0 0 0 0 0 -96 -12*L 96 -12*L;...
        0 0 0 0 0 0 12*L L^2 -12*L 2*L^2];
% Delete first 2 columns and rows of both matricies (reduced stiffness/mass
% matrix)
M4(1:2,:) = [];
M4(:,1:2) = [];
K4(1:2,:) = [];
K4(:,1:2) = [];

[eigvec4,eigval4] = eig(K4,M4); % solves for eigenvalue/eigenvectors -> eigenvalues net us omega^2
omega_4 = sqrt(diag(eigval4)); % [rad/s]
freq_4 = omega_4(1:3)./(2*pi)

%% Load in data
% all_2min = table2array(readtable('test_2min_all_3'));
% all_5min = table2array(readtable('test_5min_all_1'));
% nose = table2array(readtable('test_2min_nose_2'));
% tail = table2array(readtable('test_2min_tail_3'));
% wing = table2array(readtable('test_2min_wing_2'));


%% Plotting Modes
eigval2 = diag(eigval2);
eigval4 = diag(eigval4);

% calculate new eigen vectors for each frequency
% 2 finite elements:
for i=1:3
    ev = [0;0;eigvec2(:,i)]; % add zero values to start of matrix
    [xtemp,vtemp] = ploteigenvector(L,ev,2,12,1);
    % plot
    figure()
    plot(xtemp,vtemp);
    grid on
    plottitle = strcat('2 Element Mode shape at ' + string(freq_2(i)) + ' Hz');
    title(plottitle);
    xlabel("Length along bar [in]")
    ylabel("Displacement")
end

% 4 finite elements:
for i=1:3
    ev = [0;0;eigvec4(:,i)]; % add zero values to start of matrix
    [xtemp,vtemp] = ploteigenvector(L,ev,4,12,1);
    % plot 
    figure()
    plot(xtemp,vtemp);
    grid on
    plottitle = strcat('4 Element Mode shape at ' + string(freq_4(i)) + ' Hz');
    title(plottitle);
    xlabel("Length along bar [in]")
    ylabel("Displacement")
end
    
%% FUNCTIONS
function [x,v] = ploteigenvector(L,ev,ne,nsub,scale)
    %declare local variables here if required by language
    nv=ne*nsub+1;
    Le=L/ne; 
    dx=Le/nsub; 
    k=1;
    x=zeros(nv,1);
    v=zeros(nv,1); %declare and set to zero plot arrays

    for e=1:ne  %loop over elements
        xi=Le*(e-1);
        vi=ev(2*e-1);
        qi=ev(2*e);
        vj=ev(2*e+1);
        qj=ev(2*e+2);
        for n=1:nsub %loop over subdivisions
            xk=xi+dx*n; 
            num=(2*n-nsub)/nsub; %isoP coordinate
            vk=scale*(0.125*(4*(vi+vj)+2*(vi-vj)*(num^2-3)*num+Le*(num^2-1)*(qj-qi+(qi+qj)*num))); %Hermitian interpolant
            k = k+1;
            x(k)=xk;
            v(k)=vk; %build plot functions
    %plot v (vertical) vs x (horizontal) -- language dependent
        end
    end
end

