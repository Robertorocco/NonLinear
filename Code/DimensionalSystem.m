%% Brief analysis of system equation in 1.1 
% The study shows that system exibith same behaviur of the previus
close all
clear

%% parameters. All >0

%Clearing Rates (Negative FB)
b1=0.6;
b2=0.4;
b3=0.3;
%Coupling Effects
k1=1;
k2=1;
%Non Linearity function
alpha=1;
K=0.1;
n=100;  %Order of equation (8 stable, 9 unstable)
 


%% Equilibria 

%xi_0 = [0; 0; 0];

%Initial condition for Periodical solution
xi_0 = [1; 10; 10];



handle = @(xi) goodwin_oscillator(xi, b1,b2,b3,k1,k2,alpha,K,n);
x_eq = fsolve(handle, xi_0)

%Dimensional model
%% Solver to see trajectory (PPlane in 3d)

% Time span for the simulation: from t0 to tf
%tspan = [0 250];  %doing this way u let ode decide the time step size
tspan = 0:0.001: 100;  
[t, xi] = ode45(@(t, xi) goodwin_oscillator( xi, b1, b2, b3, alpha, K, n, k1, k2), tspan, xi_0);

xi_01=xi_0/1.5;
[t, xi2] = ode45(@(t, xi2) goodwin_oscillator( xi2, b1, b2, b3, alpha, K, n, k1, k2), tspan, xi_01);


xi_02=xi_0/1.3;
[t, xi3] = ode45(@(t, xi3) goodwin_oscillator( xi3, b1, b2, b3, alpha, K, n, k1, k2), tspan, xi_02);

%Plot
figure;
plot3(xi(:,1), xi(:,2), xi(:,3), 'b-', 'LineWidth', 1); % Blue trajectory
hold on;
plot3(xi2(:,1), xi2(:,2), xi2(:,3), 'r-', 'LineWidth', 1); % Red trajectory
plot3(xi3(:,1), xi3(:,2), xi3(:,3), 'g-', 'LineWidth', 1); % Green trajectory
% Adjust axis limits
xlim([min([xi(:,1); xi2(:,1); xi3(:,1)]), max([xi(:,1); xi2(:,1); xi3(:,1)])]);
ylim([min([xi(:,2); xi2(:,2); xi3(:,2)]), max([xi(:,2); xi2(:,2); xi3(:,2)])]);
zlim([min([xi(:,3); xi2(:,3); xi3(:,3)]), max([xi(:,3); xi2(:,3); xi3(:,3)])]);
% Add labels and title
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
title('3D Trajectories of xi, xi2, and xi3');
grid on;
legend('xi', 'xi2', 'xi3'); % Add a legend
view(3); % Ensure 3D view

figure;
plot3(xi(20000:end,1), xi(20000:end,2), xi(20000:end,3), 'k-', 'LineWidth', 0.5 ,'Color', 'b'); hold on
plot3(xi2(20000:end,1), xi2(20000:end,2), xi2(20000:end,3), 'k-', 'LineWidth', 0.5,'Color', 'r');
plot3(xi3(20000:end,1), xi3(20000:end,2), xi3(20000:end,3), 'k-', 'LineWidth', 0.5,'Color', 'g'); 

xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
grid on;
title('Zoomed Phase Portrait');
% 
% figure;
% 
% plot3(xi(100000:end,1), xi(100000:end,2), xi(100000:end,3), 'k-', 'LineWidth', 0.1 ,'Color', 'b'); hold on
% plot3(xi2(100000:end,1), xi2(100000:end,2), xi2(100000:end,3), 'k-', 'LineWidth', 0.1,'Color', 'r');
% plot3(xi3(100000:end,1), xi3(100000:end,2), xi3(100000:end,3), 'k-', 'LineWidth', 0.1,'Color', 'g'); 
% 
% xlabel('x_1');
% ylabel('x_2');
% zlabel('x_3');
% grid on;
% title('Extra-Zoomed Phase Portrait');

figure;
plot(xi(:,1), xi(:,2)) % 'o-' adds markers and lines
xlabel('x(:,1)') 
ylabel('x(:,2)') 
title('Plot of x(:,1) vs x(:,2)')
grid on

%Highlight periodicy (There is a limit cycle)
% Assume [t, xi] is obtained from ode45, where xi(:,1) is x1(t)
figure;
plot(t, xi(:,3), 'r-', 'LineWidth', 2); 
hold on;

% Detect peaks in x1(t)
[pks, locs] = findpeaks(xi(:,3), t);
plot(locs, pks, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);

xlabel('Time');
ylabel('x_3(t)');
title('Oscillatory Behavior of x_3 with Peak Markers');
legend('x_3(t)','Peaks');
grid on;

%% Poincare
C = 11.296; % Define the plane x3 = C
tolerance = 0.05; % Small tolerance for numerical accuracy

% Initialize figure
figure;
hold on; % Ensure all scatters are on the same plot

% Define datasets and colors
datasets = {xi, xi2, xi3};
colors = {'b', 'r', 'g'};
labels = {'xi', 'xi2', 'xi3'};

for k = 1:length(datasets)
    % Extract current dataset (xi, xi2, or xi3)
    current_data = datasets{k};
    
    % Find indices where x3 crosses C
    crossing_indices = find((current_data(1:end-1,3) - C) .* (current_data(2:end,3) - C) < 0);
    
    % Interpolate crossing points
    x_poincare = [];
    y_poincare = [];
    
    for i = 1:length(crossing_indices)
        idx = crossing_indices(i);
        alpha = (C - current_data(idx,3)) / (current_data(idx+1,3) - current_data(idx,3));
        x_cross = current_data(idx,1) + alpha * (current_data(idx+1,1) - current_data(idx,1));
        y_cross = current_data(idx,2) + alpha * (current_data(idx+1,2) - current_data(idx,2));
        
        x_poincare = [x_poincare; x_cross];
        y_poincare = [y_poincare; y_cross];
    end
    
    % Scatter plot (skip first 4 points to match your original code)
    scatter(x_poincare(1:end), y_poincare(1:end), 5, colors{k}, 'filled');
end

% Add labels, title, and legend
xlabel('x_1');
ylabel('x_2');
title('Poincaré Section at x_3 = C (All Datasets)');
legend(labels);
grid on;

%% Jacobian
J = jacobian(b1,b2,b3,k1,k2,alpha,K,n, x_eq);
eigenvalues = eig(J);




%% functions
function xi_dot = goodwin_oscillator(xi, b1,b2,b3,k1,k2,alpha,K,n)

xi_dot = [-b1*xi(1)+alpha/(1+K*xi(3)^n);...
        -b2*xi(2)+k1*xi(1);...
        -b3*xi(3)+ k2*xi(2);];

end

function J = jacobian(b1,b2,b3,k1,k2,alpha,K,n, xi_eq)

J = [ -b1,  0, -alpha*K*n*xi_eq(3)^(n-1)/(1+K*xi_eq(3)^n)^2;
       k1, -b2,  0;
        0,  k2, -b3 ];


end



