%% This script contains: Open loop trajectories, Poincarè plane
clear 
close all


%% parameters. All >0

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
 

%% Adimensional Model

x_0 = [1; 7; 7];

x_01 = [2;2;2]; %Check that trajectory converges on LC even from inside

x_02= [2;2;5];

%% Trajectory
% Time span for the simulation: from t0 to tf
%tspan = [0 250];  %doing this way u let ode decide the time step size
timestep=0.01;
tspan = 0:timestep: 30;  
[t, x] = ode45(@(t, x) adimensional_goodwin_oscillator(x,a,n), tspan, x_0);
[t, x1] = ode45(@(t, x1) adimensional_goodwin_oscillator(x1,a,n), tspan, x_01);
[t, x2] = ode45(@(t, x2) adimensional_goodwin_oscillator(x2,a,n), tspan, x_02);
%% Equilibria 
handle = @(zi) adimensional_goodwin_oscillator(zi,a,n);
x_eq = fsolve(handle, x_0);
J=jacobian(x,a,n);
eig(J)

%% Poincaré Section at x3 = 2.3
% Define the plane x3 = C
C = 2.3; % The plane where x3 crosses this value
tolerance = 0.01; % Small tolerance for numerical accuracy

% Find indices where x3 crosses C (from negative to positive direction only)
crossing_indices = [];
for i = 1:length(t)-1
    % Look for crossings from below to above (positive direction)
    if (x(i,3) < C) && (x(i+1,3) > C)
        crossing_indices = [crossing_indices; i];
    end
end

% Interpolate to get precise crossing points
x_poincare = [];
y_poincare = [];

for i = 1:length(crossing_indices)
    idx = crossing_indices(i);
    
    % Linear interpolation for more precise crossing
    alpha = (C - x(idx,3)) / (x(idx+1,3) - x(idx,3));
    x_cross = x(idx,1) + alpha * (x(idx+1,1) - x(idx,1));
    y_cross = x(idx,2) + alpha * (x(idx+1,2) - x(idx,2));
    
    x_poincare = [x_poincare; x_cross];
    y_poincare = [y_poincare; y_cross];
end
% Create a grid for the plane x3 = C (where C = 2.3)
x1_limits = [min(x(:,1)), max(x(:,1))];
x2_limits = [min(x(:,2)), max(x(:,2))];
[X1, X2] = meshgrid(linspace(x1_limits(1), x1_limits(2), 10), ...
                    linspace(x2_limits(1), x2_limits(2), 10));
% Create a constant value matrix for x3 = 2.3
X3 = C * ones(size(X1));
%% Plot 3D Phase Portrait 
figure1 = figure('Position', [100, 100, 800, 600]);  % Reduced from maximized to a more moderate size
axes1 = axes('Parent', figure1,...
    'Position', [0.13, 0.11, 0.775, 0.815]);  % Standard axes position
hold(axes1, 'on');

% Reduced linewidth to 1.8 (from 2.5)
plot3(x(:,1), x(:,2), x(:,3), 'k-', 'LineWidth', 1, 'Color', 'b', 'DisplayName', 'Trajectory rooted in x_0'); 
hold on

% Reduced linewidth to 1.8 for the second trajectory as well
plot3(x2(:,1), x2(:,2), x2(:,3), 'k-', 'LineWidth', 1, 'Color', 'r', 'DisplayName', 'Trajectory rooted in x_1'); 
hold on

% Slightly reduced marker size to 8 (from 10)
plot3(x_eq(1), x_eq(2), x_eq(3), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'Saddle-Focus');

grid on;
%title('3D Phase Portrait of Goodwin Oscillator', 'FontSize', 16);  % Reduced font size
view(axes1, [-25.2940514884443 15.3558595531236]);

% Moderate legend size
legend1 = legend(axes1, 'show');
set(legend1, 'FontSize', 12, 'Location', 'northeast');  % Simplified legend positioning

% Adjusted label sizes
zlabel('Inhibitor Concentration', 'FontSize', 14);
ylabel('Protein Concentration', 'FontSize', 14, 'Position', [-1.5, 3, 0.3]);
xlabel('mRNA Concentration', 'FontSize', 14);

% Moderate axis properties
set(gca, 'LineWidth', 1.2);  % Reduced from 1.5
set(gca, 'FontSize', 11);    % Reduced from 12
%exportgraphics(gcf, 'fig4.pdf', 'ContentType', 'vector');
%% Time Series Plot with Enhanced Visuals
figure('Position', [100, 100, 900, 700]);  % Larger figure size

% First subplot: x1
subplot(3, 1, 1);
plot(t, x(:,1), 'r-', 'LineWidth', 2.5);  % Increased linewidth
ylabel('x_1', 'FontSize', 14, 'FontWeight', 'bold');
title('Time Series of Adimensional Goodwin Oscillator', 'FontSize', 16, 'FontWeight', 'bold');
legend('x_1', 'FontSize', 12, 'Location', 'best');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);  % Enhanced axes

% Second subplot: x2
subplot(3, 1, 2);
plot(t, x(:,2), 'g-', 'LineWidth', 2.5);  % Increased linewidth
ylabel('x_2', 'FontSize', 14, 'FontWeight', 'bold');
legend('x_2', 'FontSize', 12, 'Location', 'best');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);  % Enhanced axes

% Third subplot: x3
subplot(3, 1, 3);
plot(t, x(:,3), 'b-', 'LineWidth', 2.5);  % Increased linewidth
xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('x_3', 'FontSize', 14, 'FontWeight', 'bold');
legend('x_3', 'FontSize', 12, 'Location', 'best');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);  % Enhanced axes
%exportgraphics(gcf, 'fig0.pdf', 'ContentType', 'vector');
%% Enhanced Poincaré Section Plot
figure('Position', [100, 100, 700, 600]);  % Larger figure size
scatter(x_poincare, y_poincare, 15, 'filled', 'MarkerFaceColor', 'g');  % Increased marker size significantly
xlabel('x_1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('x_2', 'FontSize', 16, 'FontWeight', 'bold');
title('Poincaré Section at x_3 = C', 'FontSize', 18, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5);  % Enhanced axes
%%exportgraphics(gcf, 'fig1.pdf', 'ContentType', 'vector');

%% Enhanced 3D Intersection Plot
figure('Position', [100, 100, 800, 700]);  % Larger figure size
plot3(x(:,1), x(:,2), x(:,3), 'b', 'LineWidth', 2.5);  % Increased linewidth
hold on;

 % Plot the plane with transparency - enhanced colors
surf(X1, X2, X3, 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'r');

% Add poincaré points to the 3D plot for clarity
scatter3(x_poincare, y_poincare, C*ones(size(x_poincare)), 50, 'filled', 'MarkerFaceColor', 'g');

% Labels and formatting
xlabel('x_1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('x_2', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('x_3', 'FontSize', 16, 'FontWeight', 'bold');
title('3D Trajectory and Plane x_3 = C', 'FontSize', 18, 'FontWeight', 'bold');
grid on;
axis tight;
view(3);
set(gca, 'FontSize', 14, 'LineWidth', 1.5);  % Enhanced axes

% Add a legend
legend('Trajectory', 'Plane x_3 = C', 'Poincaré Points', 'FontSize', 14, 'Location', 'best');
campos([-10.3704089714954 -48.2274326201413 15.8218496302632]);
camtarget([4.53883034134992 3.83650654612592 4.25124582485103]);
camup([0 0 1]);
%exportgraphics(gcf, 'fig2.pdf', 'ContentType', 'vector');
%% Control Data 
Ampdes=2.2 ;
Freqdes= 2*pi/6; 
Phasedes= 3*pi/2;
Biasdes=3.3;
 %Comparison
reference = Ampdes * sin(Freqdes*t + Phasedes) + Biasdes;
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





