%% Control
close all
clear
%% Reference

Ampdes=2.2 ;
Freqdes= 2*pi/6; 
Phasedes= 3*pi/2;
Biasdes=3.3;
%% Params
%Clearing Rates (Negative FB)
k=0.1;
a=k^-3;
%Coupling Effects
% k1=1;
% k2=1;
%Non Linearity function
%alpha=1;
%K=0.1;
n=9;  %Order of equation (8 stable, 9 unstable)

%% Trajecotry and equilibria
x_0 = [1; 7; 7];
% Time span for the simulation: from t0 to tf
%tspan = [0 250];  %doing this way u let ode decide the time step size
timestep=0.01;
tspan = 0:timestep: 30;  
[t, x] = ode45(@(t, x) adimensional_goodwin_oscillator(x,a,n), tspan, x_0);
handle = @(zi) adimensional_goodwin_oscillator(zi,a,n);
x_eq = fsolve(handle, x_0);
J=jacobian(x,a,n);
eig(J)
%% Linearized System

%Derivative of the field evaluated in equilibrium:
A = [ -1,       0,    -a*n*x_eq(3)^(n-1)/(1 + x_eq(3)^n)^2;
       1,      -1,                             0;
       0,       1,                            -1];

%% Feedback Linearization
%Inner Stability
% Define the bisector line (x2dot = -x2)
x2 = linspace(-10, 10, 100); % Range for x2
x2dot = -x2; % Bisector equation: x2dot = -x2

% Plot the bisector
figure;
plot(x2, x2dot, 'b-', 'LineWidth', 2); % Blue continuous line
hold on;

% Add axes for clarity
plot([0, 0], [-10, 10], 'k-', 'LineWidth', 0.5); % Vertical axis (x2dot)
plot([-10, 10], [0, 0], 'k-', 'LineWidth', 0.5); % Horizontal axis (x2)
grid on;
delta = 0.5; % Arrow length

quiver(3, 0, -delta, 0, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 100);

quiver(-3, 0, delta, 0, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 100);
% 
% quiver(0, 3, 0, -delta, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
% 
% quiver(0, -3, 0, delta, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
% Labels and formatting
plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
xlabel('$x_2$', 'Interpreter', 'latex');
ylabel('$\dot{x}_2$', 'Interpreter', 'latex');
title('$\dot{x}_2 = -x_2$', 'Interpreter', 'latex');
legend('$\dot{x}_2 = -x_2$', 'Location', 'northeast', 'Interpreter', 'latex');
axis([-5 5 -5 5]);
axis equal;
%% Sliding 

%Sliding Surface:

% sigma = p1(z1-zd), gz=(1,0,0)
%Lie_gz(sigma)= [p1, 0, 0] * [a; 0; 0;] = p1 =/0 so i can take 1D sigma.

k=1;
p1=1;


%% DATA ELABORATION:

%% 1 Settling time of the error
% Extract the error signal and time vector from Simulink output.
% If out.error is a timeseries object, extract its Data and Time fields.
if isa(out.error, 'timeseries')
    error_signal = out.error.Data;
    time_vector = out.error.Time;
else
    % Otherwise assume out.error is already a numeric vector.
    error_signal = out.error;
    time_vector = out.tout; % adjust if your time vector is stored differently
end
error_signal=squeeze(error_signal);
% Compute the amplitude A as the maximum absolute error value.
var  = max(abs(error_signal));

% Define the threshold as 5% of the amplitude.
threshold = 0.05 * var;
% threshold = 0.05 * 3;

% Find the last time index where the absolute error exceeds or equals the threshold.
last_idx_outside = find(abs(error_signal) >= threshold, 1, 'last');

% Determine the settling time: the time immediately after the last violation.
if isempty(last_idx_outside)
    % The signal is always within bounds.
    settling_time = time_vector(1);
else
    if last_idx_outside < length(time_vector)
        settling_time = time_vector(last_idx_outside + 1);
    else
        settling_time = time_vector(last_idx_outside);
    end
end

% Display the settling time
fprintf('Settling time: %.3f seconds\n', settling_time);

% Plot the absolute error with the threshold and settling time indicated.
figure;
plot(time_vector, abs(error_signal), 'b', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 1.2, 'Label', 'Tolerance Bound');
xline(settling_time, 'g--', 'LineWidth', 1.2, 'Label', sprintf('t_s = %.3f s', settling_time));
xlabel('Time (s)');
ylabel('|Error|');
title('Error','FontSize', 14);
grid on;

%% PLOT

%Comparison
reference = Ampdes * sin(Freqdes*t + Phasedes) + Biasdes;
figure;
plot(t, reference, 'r--', 'LineWidth', 1.5);  % Plot sinusoidal signal in blue
hold on;
plot(t, x(:,1), 'b-', 'LineWidth', 1.5);         % Plot ODE solution in red dashed line
xlabel('Time (s)');
ylabel('Signal');
legend('Reference Signal', 'x1 trajectory');
title('Comparison of Sinusoidal Signal and ODE Solution');
grid on;

%Plot the input signal 
% u    = out.u;  
% figure;
% plot(u.Time, u.Data, 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Input u');
% title('Input Signal', 'FontSize', 14);
% grid on;

% Plot the state vector (Separated plot)
if isa(out.x, 'timeseries')
    X = out.x.Data;
    Time = out.x.Time;
else
X   = out.x;
Time=out.x.Time;
end
%Correct data format
X=squeeze(X);
if size(X,2) ~= 3
    X = X.';
end

ref=out.ref;
figure;

subplot(3, 1, 1);
plot(Time, X(:,1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x_1');
title('State x_1', 'FontSize', 14);
grid on;

subplot(3, 1, 2);
plot(Time, X(:,2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x_2');
title('State x_2', 'FontSize', 14);
grid on;

subplot(3, 1, 3);
plot(Time, X(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x_3');
title('State x_3', 'FontSize', 14);
grid on;

%Unic plot:

figure;
plot(ref.Time, ref.Data(:), 'c', 'LineWidth', 1.0); 
hold on;
plot(Time, X(:,1), 'r', 'LineWidth', 1.5);  % Plot x1 in re
plot(Time, X(:,2), 'g', 'LineWidth', 1.5);  % Plot x2 in green
plot(Time, X(:,3), 'b', 'LineWidth', 1.5);  % Plot x3 in blue
 
xlabel('Time (s)');
ylabel('States');
title('States vs. Time', 'FontSize', 14);
legend('Ref','x_1', 'x_2', 'x_3');
grid on;

%Linearized Case:
%Uncomment to see linearized plot
u_lin= out.u_lin;
figure;
plot(u_lin.Time, u_lin.Data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Input u_lin');
title('Input Linearized Signal', 'FontSize', 14);
grid on;

error_lin=out.error_lin;
figure;
plot(error_lin.Time, error_lin.Data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Error');
title('Error Linearized', 'FontSize', 14);
grid on;


x_lin    = out.x_lin;
figure;

subplot(3, 1, 1);
plot(x_lin.Time, x_lin.Data(:,1), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x_1');
title('State Linearized x_1', 'FontSize', 14);
grid on;

subplot(3, 1, 2);
plot(x_lin.Time, x_lin.Data(:,2), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x_2');
title('State Linearized x_2', 'FontSize', 14);
grid on;

subplot(3, 1, 3);
plot(x_lin.Time, x_lin.Data(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel(' x_3');
title('State  Linearized x_3', 'FontSize', 14);
grid on;


%% Robust Check
% Nominal parameter values
% n_nominal = 9;      % Original constant n
% b_nominal = 0.1;    % Original constant b
% num_runs = 10; % Number of simulations
% % Generate random perturbations (±30% of nominal values)
% delta_n = 0.3 * (2*rand() - 1);  % Random value between -0.3 and +0.3
% delta_b = 0.3 * (2*rand() - 1);
% 
% % Calculate perturbed parameters
% n_perturbed = n_nominal * (1 + delta_n);  % Perturbed n
% b_perturbed = b_nominal * (1 + delta_b);  % Perturbed b
% 
% % Preallocate storage for results
% simResults = struct('time', cell(num_runs,1), 'x1', cell(num_runs,1), ...
%                    'x2', cell(num_runs,1), 'x3', cell(num_runs,1));
% 
% % Configure simulations
% for i = 1:num_runs
%     % Generate parameter uncertainties (±40%)
%     delta_n = 0.5*(2*rand() - 1);
%     delta_b = 0.5*(2*rand() - 1);
%     % Calculate perturbed parameters
% 
%     % Create simulation input object
%     simIn(i) = Simulink.SimulationInput('smc_trackingRobust'); %Insert the correct simulink 
%     simIn(i) = simIn(i).setVariable('n_perturbed', n_nominal*(1 + delta_n));
%     simIn(i) = simIn(i).setVariable('b_perturbed', b_nominal*(1 + delta_b));
% end
% 
% % Run simulations
% simOut = sim(simIn, 'ShowProgress', 'on'); % Works without Parallel Toolbox
% 
% % Extract data from timeseries
% for i = 1:num_runs
%     % Get the logged timeseries - replace 'x' with your actual logged variable name
%     x_ts = simOut(i).get('x'); % Assumes state is logged as variable 'x'
% 
%     % Store results
%     simResults(i).time = x_ts.Time;
%     simResults(i).x1 = x_ts.Data(:,1);
%     simResults(i).x2 = x_ts.Data(:,2); 
%     simResults(i).x3 = x_ts.Data(:,3);
% end
% 
% % Plot results
% figure('Color', 'white', 'Position', [100 100 1200 800]);
% cmap = jet(num_runs); % Colormap for distinct trajectories
% 
% % Subplot 1: x1
% subplot(3,1,1);
% hold on;
% for i = 1:num_runs
%     plot(simResults(i).time, simResults(i).x1, 'Color', cmap(i,:), 'LineWidth', 1.5);
% end
% title('State x_1 Evolution');
% xlabel('Time');
% ylabel('x_1');
% grid on;
% 
% % Subplot 2: x2
% subplot(3,1,2);
% hold on;
% for i = 1:num_runs
%     plot(simResults(i).time, simResults(i).x2, 'Color', cmap(i,:), 'LineWidth', 1.5);
% end
% title('State x_2 Evolution');
% xlabel('Time');
% ylabel('x_2');
% grid on;
% 
% % Subplot 3: x3
% subplot(3,1,3);
% hold on;
% for i = 1:num_runs
%     plot(simResults(i).time, simResults(i).x3, 'Color', cmap(i,:), 'LineWidth', 1.5);
% end
% title('State x_3 Evolution');
% xlabel('Time');
% ylabel('x_3');
% grid on;
% 
% % Add colorbar for parameter variation reference
% colormap(jet(num_runs));
% colorbar('Ticks', linspace(0,1,num_runs), 'TickLabels', 1:num_runs, ...
%          'Position', [0.93 0.15 0.02 0.7]);

%% Function 
function lie=lie_derivative(grad, fx)
lie= grad * fx;
end

%% Function
function x_dot = adimensional_goodwin_oscillator(x,a,n)
u= 0;
x_dot = [-x(1)+a/(1+x(3)^n) + a*u;...
        -x(2)+x(1);...
        -x(3)+ x(2);];

end

function J = jacobian(x,a,n)

J = [ -1,       0,    -a*n*x(3)^(n-1)/(1 + x(3)^n)^2;
       1,      -1,                             0;
       0,       1,                            -1];
end
