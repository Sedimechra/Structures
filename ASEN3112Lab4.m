%% ASEN 3112 Lab 4

% Tristan Workman
% Joe Moser
% Graham Kersey
% Maklen Estrada
% Skylar Edwards


%% Housekeeping
clc;
clear all;

%% Calling dataAnaylze function, producing known deflections and associated force
[d_thin, F_thin,thinF] = dataAnalyze('ThinBar_FA2020.mat',"Thinbar_FA2020");
[d_square,F_square,squareF] = dataAnalyze('Squarerod_FA2020.mat',"Squarebar_FA2020");

%% Calling predict function to get predicted forces for each displacement
% first predicting for theoretical Pcr
[d_thin, F_thinP] = predict(139.0054,10.8);
[d_square,F_squareP] = predict(286.6831,11.85);
% second predicting for experimentally determined Pcr
[dsquare,F_thinP_E] = predict(F_thin(2),10.8);
[d_square,F_squareP_E] = predict(F_square(2),11.85);

%% Finding deflection for post buckling
thin_D = findDeflec(10.8,35000,10000000,0.125);
square_D = findDeflec(11.85,35000,10000000,0.25);

%% plotting displacement verssus force
figure()
plot(d_thin,F_thin,"Linewidth",2)
hold on
plot(d_thin,F_thinP,"Linewidth",2)
plot(d_thin,F_thinP_E,"--","Linewidth",2)
xline(thin_D,"--","Linewidth",2)
% ylim([0 210])
yticks(linspace(0,150,7))
grid on
title("Predicted versus Experimental for Thin Specimen")
xlabel("Deflection [in]")
ylabel("Force [lbf]")
legend(["Experimental" "Theoretical Pcr" "Experimental Pcr" "Plastic Post Buckling"],"Location","Best")

figure()
plot(d_square,F_square,"Linewidth",2)
hold on
plot(d_square,F_squareP,"Linewidth",2)
plot(d_square,F_squareP_E,"--","Linewidth",2)
xline(square_D,"--","Linewidth",2)
% ylim([0 210])
grid on
title("Predicted versus Experimental for Square Specimen")
xlabel("Deflection [in]")
ylabel("Force [lbf]")
legend(["Experimental" "Theoretical Pcr" "Experimental Pcr" "Plastic Post Buckling"],"Location","Best")

%% Function find deflection
% finds lateral deflection at plastic post buckling given mat properties
function deflect = findDeflec(L,sigma,E,w)
    epsilon = sigma/E;
    deflect = (2 * epsilon * L^2)/(pi^2 * w);
end

%% Function predict force
% takes in critical load and length of specemine and outputs predicted
% force relative to deflection
function [deflecVec, predForce] = predict(Pcr, L)
    % create a vector of known deflections
    deflecVec = [linspace(0,1,17) linspace(1,2,9)];
    % use Pcr and L to get predicted forces for each deflection
    predForce = Pcr .* (1 + (pi^2/(8 * L^2)) .* deflecVec.^2);
end

%% Function dataAnalyze
% takes in matfile and the name of the matfile (for indexing the structure)
% and outputs a vector of the known deflections and their associated force
% converted from voltage
function [deflecVec, adjForce, finalF] = dataAnalyze(matFile,name)
    % Load file
    file = load(matFile);
    % pull voltage/deflection vectors from the file
    voltage = file.(name)(:,2);
    deflection = file.(name)(:,3);
    % create a vector of known deflections
    deflecVec = [linspace(0,1,17) linspace(1,2,9)];
    % preallocate voltage average holder;
    voltAv = [];
    voltF = [];
    for j = 1:26
        % reset holder and count to zero
        holder = 0;
        count = 0;
        for i = 1:size(voltage)
           if deflecVec(j) == deflection(i)
               % find each same deflection and average associated voltage
               holder = holder + voltage(i);
               final = voltage(i);
               count = count + 1;
           end
        end
        voltAv(end+1) = holder./count;
        voltF(end+1) = final;
    end
    % convert voltage average to force
    adjForce = ((voltAv .* 1000)./2.37);
    finalF = ((voltF .* 1000)./2.37);
end