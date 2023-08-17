clc;
clear;
V1=1;
delta_1=0;
V2=1;
delta_2=0;
S_base = 100;

% Define line data for the 5-bus system
% Format: [From Bus, To Bus, Resistance (pu), Reactance (pu), Shunt Admittance (pu)]
lineData = [
    1 2 0.12 0.16 0;
];

% Calculate the number of buses
numBuses = max(max(lineData(:, 1)), max(lineData(:, 2)));

% Initialize the Y-Bus matrix
YBus = zeros(numBuses, numBuses);

% Populate the Y-Bus matrix
for i = 1:size(lineData, 1)
    fromBus = lineData(i, 1);
    toBus = lineData(i, 2);
    resistance = lineData(i, 3);
    reactance = lineData(i, 4);
    shuntAdmittance = lineData(i, 5);
    
    % Off-diagonal elements
    YBus(fromBus, toBus) = -1 / (resistance + 1i * reactance);
    YBus(toBus, fromBus) = YBus(fromBus, toBus);
    shuntAdmittance = 1i * shuntAdmittance;
    
    % Diagonal elements
     YBus(fromBus, fromBus) = YBus(fromBus, fromBus) - YBus(fromBus, toBus) + shuntAdmittance;
     YBus(toBus, toBus) = YBus(toBus, toBus) - YBus(fromBus, toBus) + shuntAdmittance;   
end

disp('   it,         V2,      delta_degrees');

for iterations=1:10
% Power generated/consumed at each bus
P_MW = -100';  % negative for power consumed
Q_MVAr = -50;  % negative for reactive power consumed

% Convert power to per unit
P = P_MW / S_base;  % Real power in pu
Q = Q_MVAr / S_base;  % Reactive power in pu

Pcalc = V1 * V2 * abs(YBus(1,2)) * cos (angle(YBus(1,2))- delta_2 + delta_1) + V2* V2 * abs(YBus(2,2)) *cos (angle(YBus(2,2)));
Qcalc = - V1 * V2 * abs(YBus(1,2)) * sin(angle(YBus(1,2))- delta_2 + delta_1) - V2* V2 * abs(YBus(2,2))*sin(angle(YBus(2,2)));


Jacobian(1,1) = V1 * V2 * abs(YBus(1,2)) * sin (angle(YBus(1,2))- delta_2 + delta_1);
Jacobian(1,2) = V1 * abs(YBus(1,2)) * cos (angle(YBus(1,2))- delta_2 + delta_1) + 2* V2 * abs(YBus(2,2)) *cos (angle(YBus(2,2)));
Jacobian(2,1) = V1 * V2 * abs(YBus(1,2)) * cos(angle(YBus(1,2))- delta_2 + delta_1);
Jacobian(2,2) = - V1 * abs(YBus(1,2)) * sin (angle(YBus(1,2))- delta_2 + delta_1) - 2* V2 * abs(YBus(2,2))*sin(angle(YBus(2,2)));


Change_X = inv(Jacobian) * [P - Pcalc;
    Q - Qcalc];

X = [delta_2;
    V2];

Ans = Change_X + X;

V2= V2 + Change_X(2,1);
delta_2_1 = rad2deg(delta_2 + Change_X(1,1));
delta_2 = delta_2 + Change_X(1,1);


disp([iterations, V2, delta_2_1]);

end






