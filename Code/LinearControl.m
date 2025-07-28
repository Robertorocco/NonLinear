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
out=sim('LinearControl_sim.slx');


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
if isa(out.x, 'timeseries')
    X = out.x.Data;
    Time = out.x.Time;
else
X   = out.x;
Time=out.x.Time;
end

% Check and extract u_lin consistently
if isa(out.u_lin, 'timeseries')
    u_lin_data = out.u_lin.Data;
    u_lin_time = out.u_lin.Time;
else
    u_lin_data = out.u_lin;
    u_lin_time = out.time; % or out.tout depending on your simulation output structure
end



%Correct data format
X=squeeze(X);
if size(X,2) ~= 3
    X = X.';
end

ref=out.ref;

%% PLOT
% Display the settling time
fprintf('Settling time: %.3f seconds\n', settling_time);

% Plot the absolute error with the threshold and settling time indicated.
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

% Create an inset for the zoom on the second half of the simulation
half_time_index = ceil(length(time_vector)/2);
zoom_time = time_vector(half_time_index:end);
zoom_error = abs(error_signal(half_time_index:end));

% Create inset axes in the upper right corner
ax_inset = axes('Position', [0.65, 0.65, 0.25, 0.25]);
plot(zoom_time, zoom_error, 'b', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 1);
if settling_time >= zoom_time(1)
    xline(settling_time, 'g--', 'LineWidth', 1);
end
grid on;
title('Zoom View', 'FontSize', 10);
% Adjust axes limits if needed
% xlim([zoom_time(1), zoom_time(end)]);
% Optional: Adjust y-limit to focus on relevant error range
% ylim([0, max(zoom_error)*1.1]);
%exportgraphics(gcf, 'DLinfig3.pdf', 'ContentType', 'vector');

%Comparison
% reference = Ampdes * sin(Freqdes*t + Phasedes) + Biasdes;
% figure;
% plot(t, reference, 'r--', 'LineWidth', 1.5);  % Plot sinusoidal signal in blue
% hold on;
% plot(t, x(:,1), 'b-', 'LineWidth', 1.5);         % Plot ODE solution in red dashed line
% xlabel('Time (s)');
% ylabel('Signal');
% legend('Reference Signal', 'x1 trajectory');
% %title('Comparison of Sinusoidal Signal and ODE Solution');
% grid on;
% %exportgraphics(gcf, 'Lin_fig1.pdf', 'ContentType', 'vector');
%Plot the input signal 
% u    = out.u;  
% figure;
% plot(u.Time, u.Data, 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Input u');
% title('Input Signal', 'FontSize', 14);
% grid on;

% Plot the state vector (Separated plot)

% Create a wider figure
figure('Position', [100, 100, 900, 600]);

% First subplot with reference signal
subplot(3, 1, 1);
plot(Time, X(:,1), 'b', 'LineWidth', 1.5);  % Blue line for state
hold on;
plot(ref.Time, ref.Data(:), 'r--', 'LineWidth', 1.0);  % Red dashed line for reference
xlabel('Time (s)');
ylabel('x_1');
grid on;
legend('State', 'Reference', 'Location', 'southeast');  % Moved to bottom right


% Second subplot
subplot(3, 1, 2);
plot(Time, X(:,2), 'b', 'LineWidth', 1.5);  % Blue color
xlabel('Time (s)');
ylabel('x_2');
grid on;

% Third subplot
subplot(3, 1, 3);
plot(Time, X(:,3), 'b', 'LineWidth', 1.5);  % Blue color
xlabel('Time (s)');
ylabel('x_3');
grid on;
sgtitle('Real system state', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');

%exportgraphics(gcf, 'DLinfig2.pdf', 'ContentType', 'vector');
%Linearized Case:
% Now plot with the extracted data
figure;
plot(u_lin_time, u_lin_data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('u_{lin}');
title('Input Linearized Signal', 'FontSize', 14);
grid on;


error_lin=out.error_lin;
figure;
plot(error_lin.Time, error_lin.Data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Error');
title('Error Linearized', 'FontSize', 14);
grid on;

% Check and extract u consistently
if isa(out.u, 'timeseries')
    u_data = out.u.Data;
    u_time = out.u.Time;
else
    u_data = out.u;
    u_time = out.time; % or out.tout depending on your simulation output structure
end

% Now plot with the extracted data - using blue color
figure;
plot(u_time, u_data, 'b', 'LineWidth', 1.5);  % Blue color for plot
xlabel('Time (s)');
ylabel('u');
title('Input Signal', 'FontSize', 14);
grid on;
%exportgraphics(gcf, 'DLinfig4.pdf', 'ContentType', 'vector');
% If multiple lines appear, you might want to add a legend
if size(u_data, 2) > 1
    legend('Component 1', 'Component 2', 'Component 3');
end

x_lin    = out.x_lin;

% Create a wider figure
figure('Position', [100, 100, 900, 600]);

% First subplot
subplot(3, 1, 1);
plot(x_lin.Time, x_lin.Data(:,1), 'b', 'LineWidth', 1.5);  % Blue line for state
hold on;
if exist('ref', 'var')  % Check if reference exists
    plot(ref.Time, ref.Data(:), 'r--', 'LineWidth', 1.0);  % Red dashed line for reference
    legend('State', 'Reference', 'Location', 'southeast');  % Bottom right legend
end
xlabel('Time (s)');
ylabel('x_{1lin}');
grid on;

% Second subplot
subplot(3, 1, 2);
plot(x_lin.Time, x_lin.Data(:,2), 'b', 'LineWidth', 1.5);  % Blue color
xlabel('Time (s)');
ylabel('x_{2lin}');
grid on;

% Third subplot
subplot(3, 1, 3);
plot(x_lin.Time, x_lin.Data(:,3), 'b', 'LineWidth', 1.5);  % Blue color
xlabel('Time (s)');
ylabel('x_{3lin}');
grid on;
sgtitle('Linear system state', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
%exportgraphics(gcf, 'DLinFig1.pdf', 'ContentType', 'vector');


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