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
out=sim('FBL_sim.slx');
%% Feedback Linearization
% %Inner Stability
% % Define the bisector line (x2dot = -x2)
% x2 = linspace(-10, 10, 100); % Range for x2
% x2dot = -x2; % Bisector equation: x2dot = -x2
% 
% % Plot the bisector
% figure;
% plot(x2, x2dot, 'b-', 'LineWidth', 2); % Blue continuous line
% hold on;
% 
% % Add axes for clarity
% plot([0, 0], [-10, 10], 'k-', 'LineWidth', 0.5); % Vertical axis (x2dot)
% plot([-10, 10], [0, 0], 'k-', 'LineWidth', 0.5); % Horizontal axis (x2)
% grid on;
% delta = 0.5; % Arrow length
% 
% quiver(3, 0, -delta, 0, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 100);
% 
% quiver(-3, 0, delta, 0, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 100);
% % 
% % quiver(0, 3, 0, -delta, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
% % 
% % quiver(0, -3, 0, delta, 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
% % Labels and formatting
% plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
% xlabel('$x_2$', 'Interpreter', 'latex');
% ylabel('$\dot{x}_2$', 'Interpreter', 'latex');
% title('$\dot{x}_2 = -x_2$', 'Interpreter', 'latex');
% legend('$\dot{x}_2 = -x_2$', 'Location', 'northeast', 'Interpreter', 'latex');
% axis([-5 5 -5 5]);
% axis equal;


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

%% PLOT

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

% Create an inset for the zoom on the second half of the simulation
half_time_index = ceil(length(time_vector)/2);
zoom_time = time_vector(half_time_index:end);
zoom_error = abs(error_signal(half_time_index:end));

% Create inset axes in the upper right corner
ax_inset = axes('Position', [0.65, 0.65, 0.25, 0.25]);
plot(zoom_time, zoom_error, 'b', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 1.0);
if settling_time >= zoom_time(1)
    xline(settling_time, 'g--', 'LineWidth', 1.0);
end
grid on;
title('Zoom View', 'FontSize', 10);

% Optional: Add a rectangle to the main plot to indicate the zoom region
axes(gca); % Switch back to main axes
rectangle('Position', [zoom_time(1), min(zoom_error), zoom_time(end)-zoom_time(1), max(zoom_error)-min(zoom_error)], ...
          'EdgeColor', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1);
exportgraphics(gcf, 'DFBLFig1.pdf', 'ContentType', 'vector');
%Plot the input signal 
u    = out.u;  
figure;
plot(u.Time, u.Data, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Input u');
title('Input Signal', 'FontSize', 14);
grid on;

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