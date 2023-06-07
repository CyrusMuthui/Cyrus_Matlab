%{
CAT 1; Group Members
1. Brenda Cheptoo - ENE211-0004/2018
2. Cyrus Muthui - ENE211-0010/2018
3. Olivia Paige Muva - ENE211-0017/2018
%}
%Line-to-line impedances of transmission lines (in ohms)
Z12 = 0.02 + 0.1i;
Z15 = 0.05 + 0.25i;
Z23 = 0.04 + 0.2i;
Z25 = 0.05 + 0.25i;
Z34 = 0.05 + 0.25i;
Z35 = 0.08 + 0.4i;
Z45 = 0.1 + 0.5i;

% Shunt admittances of buses (in Siemens)
Y12 = 0.03i;
Y15 = 0.02i;
Y23 = 0.025i;
Y25 = 0.02i;
Y34 = 0.02i;
Y35 = 0.01i;
Y45 = 0.075i;
Y13 = 0;
Y14 = 0;
Y24 = 0;
Y11 = 1 / Z12 + 1 / Z15 + Y12 + Y15;
Y22 = 1 / Z12 + 1 / Z23 + 1 / Z25 + Y12 + Y23 + Y25;
Y33 = 1 / Z23 + 1 / Z34 + 1 / Z35 + Y23 + Y34 + Y35;
Y44 = 1 / Z34 + 1 / Z45 + Y34 + Y45;
Y55 = 1 / Z15 + 1 / Z25 + 1 / Z35 + 1 / Z45 + Y15 + Y25 + Y35 + Y45;


% Number of buses
n = 5;

% Initialize the Y bus matrix
Ybus = zeros(n);

% Add diagonal elements to Y bus matrix
Ybus(1, 1) = 1 / Z12 + 1 / Z15 + Y12 + Y15;
Ybus(2, 2) = 1 / Z12 + 1 / Z23 + 1 / Z25 + Y12 + Y23 + Y25;
Ybus(3, 3) = 1 / Z23 + 1 / Z34 + 1 / Z35 + Y23 + Y34 + Y35;
Ybus(4, 4) = 1 / Z34 + 1 / Z45 + Y34 + Y45;
Ybus(5, 5) = 1 / Z15 + 1 / Z25 + 1 / Z35 + 1 / Z45 + Y15 + Y25 + Y35 + Y45;

% Add off-diagonal elements to Y bus matrix
Ybus(1, 2) = -1 / Z12;
Ybus(1, 5) = -1 / Z15;
Ybus(2, 1) = -1 / Z12;
Ybus(2, 3) = -1 / Z23;
Ybus(2, 5) = -1 / Z25;
Ybus(3, 2) = -1 / Z23;
Ybus(3, 4) = -1 / Z34;
Ybus(3, 5) = -1 / Z35;
Ybus(4, 3) = -1 / Z34;
Ybus(4, 5) = -1 / Z45;
Ybus(5, 1) = -1 / Z15;
Ybus(5, 2) = -1 / Z25;
Ybus(5, 3) = -1 / Z35;
Ybus(5, 4) = -1 / Z45;

% Display the Y bus matrix
disp('Y bus matrix:');
disp(Ybus);
% Known values
V1 = 1.05;           % Voltage magnitude at bus 1
delta1 = 0;         % Voltage angle at bus 1
Pgen2 = 0;          % Active power generated at bus 2
Qgen2 = 0;          % Reactive power generated at bus 2
Pgen3 = 0;          % Active power generated at bus 3
Qgen3 = 0;          % Reactive power generated at bus 3
Pgen4 = 0;          % Active power generated at bus 4
Qgen4 = 0;          % Reactive power generated at bus 4
Pgen5 = 0.48;        % Active power generated at bus 5
Pload1 = 0;         % Active power load at bus 1
Qload1 = 0;         % Reactive power load at bus 1
Pload2 = 0.96;       % Active power load at bus 2
Qload2 = 0.62;       % Reactive power load at bus 2
Pload3 = 0.35;       % Active power load at bus 3
Qload3 = 0.14;       % Reactive power load at bus 3
Pload4 = 0.16;       % Active power load at bus 4
Qload4 = 0.08;       % Reactive power load at bus 4
V5 = 1.02;          % Voltage magnitude at bus 5 (specified)
Qgen5 = 0;          % Reactive power generated at bus 5 (unknown)
Pload5 = 0.24;       % Active power load at bus 5
Qload5 = 0.11;       % Reactive power load at bus 5

% Create an initial guess for voltage magnitudes and angles
V2 = 1.0 + 0i;           % Voltage magnitude at bus 2 (initial guess)
delta2 = 0;         % Voltage angle at bus 2 (initial guess)
V3 = 1.0 + 0i;           % Voltage magnitude at bus 3 (initial guess)
delta3 = 0;         % Voltage angle at bus 3 (initial guess)
V4 = 1.0 + 0i;           % Voltage magnitude at bus 4 (initial guess)
delta4 = 0;         % Voltage angle at bus 4 (initial guess)
delta5 = 0;         % Voltage angle at bus 5 (initial guess)

% Define convergence criteria
max_iter = 100;     % Maximum number of iterations
tolerance = 1e-6;   % Tolerance for convergence

% Initialize iteration variables
iter = 0;
converged = false;

while ~converged && iter < max_iter
    % Store previous values for convergence check
    V2_prev = V2*cos(delta2 * 180/pi) + 1i * V2*sin(delta2 * 180/pi);
    delta2_prev = delta2;
    V3_prev = V3*cos(delta3 * 180/pi) + 1i * V3*sin(delta3 * 180/pi);
    delta3_prev = delta3;
    V4_prev = V4*cos(delta4 * 180/pi) + 1i * V4*sin(delta4 * 180/pi);
    delta4_prev = delta4;
    Qgen5_prev = Qgen5;  % Store previous value of Qgen5
    
    % Update voltage magnitudes and angles using Gauss-Seidel iteration method
    
    % Bus 2
    V2 = abs((((Pgen2 - Pload2)- (1i * (Qgen2 - Qload2))/conj(V2_prev)) - Y12*V1 - Y23*V3_prev - Y24*V4_prev - Y25*V5)/Y22);
    delta2 = angle((((Pgen2 - Pload2)- 1i * (Qgen2 - Qload2)/conj(V2_prev)) - Y12*V1 - Y23*V3_prev - Y24*V4_prev - Y25*V5)/Y22);
    
    % Bus 3
    V3 = abs((((Pgen3 - Pload3)- (1i * (Qgen3 - Qload3))/conj(V3_prev)) - Y13*V1 - Y23*V2_prev - Y34*V4_prev - Y35*V5)/Y33);
    delta3 = angle((((Pgen3 - Pload3)- 1i * (Qgen3 - Qload3)/conj(V3_prev)) - Y13*V1 - Y23*V2_prev - Y34*V4_prev - Y35*V5)/Y33);
    
    % Bus 4
    V4 = abs((((Pgen4 - Pload4)- (1i * (Qgen4 - Qload4))/conj(V4_prev)) - Y14*V1 - Y24*V2_prev - Y34*V3_prev - Y45*V5)/Y44);
    delta4 = angle((((Pgen4 - Pload4)- (1i * (Qgen4 - Qload4))/conj(V4_prev)) - Y14*V1 - Y24*V2_prev - Y34*V3_prev - Y45*V5)/Y44);
    
    % Bus 5
    Qgen5 = imag(V5 * (Y25 * V2 * exp(1i * delta2) + Y35 * V3 * exp(1i * delta3) + Y45 * V4 * exp(1i * delta4)));
    delta5 = angle((Pgen5 + 1i * Qgen5) / V5);

    % Bus 1
    S1 = V1 * conj(Y12) * exp(1i * delta1) + V2 * conj(Y12) * exp(1i * delta2);
    Pgen1 = real(S1);
    Qgen1 = imag(S1);
    
    % Check convergence
    if abs(V2 - V2_prev) < tolerance && abs(delta2 - delta2_prev) < tolerance && abs(V3 - V3_prev) < tolerance && abs(delta3 - delta3_prev) < tolerance && abs(V4 - V4_prev) < tolerance && abs(delta4 - delta4_prev) < tolerance && abs(Qgen5 - Qgen5_prev) < tolerance
        converged = true;
    end
    
    % Increment iteration count
    iter = iter + 1;
end

% Display results
disp('Convergence status:');
disp(converged);
disp('Number of iterations:');
disp(iter);
disp('Voltage magnitudes:');
disp(['V2 = ' num2str(V2)]);
disp(['V3 = ' num2str(V3)]);
disp(['V4 = ' num2str(V4)]);
disp('Voltage angles:');
disp(['delta2 = ' num2str(rad2deg(delta2)) ' degrees']);
disp(['delta3 = ' num2str(rad2deg(delta3)) ' degrees']);
disp(['delta4 = ' num2str(rad2deg(delta4)) ' degrees']);
disp(['delta5 = ' num2str(rad2deg(delta5)) ' degrees']);
disp('Reactive power:');
disp(['Qgen5 = ' num2str(Qgen5)]);
disp('Real and reactive power generated at bus 1:');
disp(['Pgen1 = ' num2str(Pgen1)]);
disp(['Qgen1 = ' num2str(Qgen1)]);