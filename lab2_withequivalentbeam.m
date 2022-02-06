%% ASEN 3112 Lab 2
% Fall 2021

% Luca Bonarrigo
% Kevin Cook
% Mikaela Felix
% Alexander Hubben
% Bennett Spengler 
% Tristan Workman

% Housekeeping
clear ALL; clc; close ALL;

%==========================================================================
%% Load Experimental Data
data = table2array(readtable('group_1-9_data'));

load_forces = data(:,1);
F0_rxn = data(:,2);
F1_rxn = data(:,3);
F2_rxn = data(:,4);
rxn_forces = data(:,2) + data(:,3) + data(:,4);
internal_rxn = data(:,5) + abs(min(data(:,5)));
displacement = data(:,6);
% 
% figure()
% scatter(load_forces,rxn_forces);
% xlabel('Load force (lb)');
% ylabel('Rxn forces (lbf)');
% figure()
% 
% scatter(load_forces,internal_forces);
% xlabel('Load force (lb)');
% ylabel('Internal forces (lbf)');
% 
% figure()
% scatter(load_forces,displacement);
% xlabel('Load force (lb)');
% ylabel('Displacement (in)');


% figure(1)
% hold on
% scatter(load_forces,F0_rxn);
% scatter(load_forces,F1_rxn);
% scatter(load_forces,F2_rxn);
% scatter(load_forces,internal_rxn);
% xlabel('Load force (lbf)');
% ylabel('Reaction forces (lbf)');
% title('Truss Reaction Forces vs Load Forces');
% 
% figure(2)
% hold on
% scatter(load_forces,displacement);
% xlabel('Load force (lbf)');
% ylabel('Displacement (in)');
% title('Truss Displacement vs Load Forces');

%% Anysys data
ANSYS = table2array(readtable('reaction_solution'));
ANSYSrxn = ANSYS([1:4],3);
maxRxn = sum(ANSYSrxn)/4.448;
slopeA = maxRxn/50;
rxnVec = slopeA * [0:1:50];
rxnPoints = slopeA * [0:10:50];
plot([0:10:50],rxnPoints,"ko","LineWidth",1.5);
hold on
ANSYSline = plot([0:1:50],rxnVec,"k","LineWidth",1.5);
plot(load_forces,F0_rxn,'o','Color','#0072BD');
plot(load_forces,F1_rxn,'o','Color','#D95319');
plot(load_forces,F2_rxn,'o','Color','#EDB120');
y0line = plot(load_forces,y0,'--','Color','#0072BD');
y1line = plot(load_forces,y1,'--','Color','#D95319');
y2line = plot(load_forces,y2,'--','Color','#EDB120');
legend([ANSYSline(1) y0line(1) y1line(1) y2line(1)],["ANSYS","F0","F1","F2"],"Location","NorthWest")
grid on
ylim([0,50])
title("ANSYS Total Reaction Forces vs Experimental Reactions")
ylabel("Reaction Force [lbf]")
xlabel("Load Force [lbf]")

%% Linear regressions
p0 = polyfit(load_forces,F0_rxn,1);
p1 = polyfit(load_forces,F1_rxn,1);
p2 = polyfit(load_forces,F2_rxn,1);
p3 = polyfit(load_forces,internal_rxn,1);

y0 = polyval(p0,load_forces); 
y1 = polyval(p1,load_forces);
y2 = polyval(p2,load_forces);
y3 = polyval(p3,load_forces);

yresid0 = y0 - F0_rxn;
yresid1 = y1 - F1_rxn;
yresid2 = y2 - F2_rxn;
yresid3 = y3 - internal_rxn;

SSresid0 = sum(yresid0.^2);
SStotal0 = (length(F0_rxn)-1) * var(F0_rxn);
SSresid1 = sum(yresid1.^2);
SStotal1 = (length(F1_rxn)-1) * var(F1_rxn);
SSresid2 = sum(yresid2.^2);
SStotal2 = (length(F2_rxn)-1) * var(F2_rxn);
SSresid3 = sum(yresid3.^2);
SStotal3 = (length(internal_rxn)-1) * var(internal_rxn);

r2_0 = 1 - SSresid0/SStotal0;
r2_1 = 1 - SSresid0/SStotal0;
r2_2 = 1 - SSresid0/SStotal0;
r2_3 = 1 - SSresid0/SStotal0;

% figure(1)
% hold on
% plot(load_forces,y0,'--','Color','#0072BD');
% plot(load_forces,y1,'--','Color','#D95319');
% plot(load_forces,y2,'--','Color','#EDB120');
% plot(load_forces,y3,'--','Color','#7E2F8E');
% hold off

pd = polyfit(load_forces,displacement,1);
yd = polyval(pd,load_forces);
yresidd = yd - displacement;
SSresidd = sum(yresidd.^2);
SStotald = (length(displacement)-1) * var(displacement);
r2_d = 1 - SSresidd/SStotald;

% legend('F0','F1','F2','Internal forces (F3D)','Linear Regression - F0',...
%     'Linear Regression - F1','Linear Regression - F2','Linear Regression - F3D');

% figure(2)
% plot(load_forces,yd,'--','Color','#0072BD');
% legend('Displacement','Linear Regression');

% equivalent beam model
E = 69e9; % elastic modulus for aluminum alloy ranges from 69 to 70 GPA
Lb = 0.250; %[m], length between nodes along longeron
L = Lb * 16;
A = pi * (3/16)^2 - pi * (3/16 - 1/16)^2;
A = A * 0.00064516;
d = 1/2 * Lb; % distance between longeron and neutral axis
Ic_outer = pi * (3/16 * 0.0254)^4 /4;
Ic_inner = pi * ((3/16 - 1/16) * 0.0254)^4 /4;
Ic = Ic_outer - Ic_inner;
Io = 4*(Ic + A * d^2);
x = 0:Lb:L;
P = 222.4;
v1 = (P/12 * x.^3 - P*L^2/16 * x) / (E * Io);
v2 = (P*L / 4 * x.^2 - P / 12 * x.^3 - 3/16*P*L^2*x - P * L^3 * (1/16 - 1/12)) / (E*Io);
v = [v1(1:9) v2(10:17)];

% figure(3)
% plot(x(1:9),v1(1:9))
% hold on
% plot(x(9:17),v2(9:17))
% hold off

% figure(3)
% plot(x,v*1000)
% title('Equivalent Beam Model Deflections','FontSize',16)
% xlabel('Distance Throughout Beam [m]','FontSize',16)
% ylabel('Deflection [mm]','FontSize',16)
