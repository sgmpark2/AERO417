clear, clc, close all
warning off
% Define the properties of the truss
E = 200e9; % Young's modulus in Pascals
A = 6e-4; % Cross-sectional area in m^2
L = 3; % Length of each element in meters
rho =7800; % Density in kg/m^3

abaqus = readtable('abaqus.txt');

% Define the nodes and elements
nodes = [0 0; L 0; 0 L;L L]; % Node coordinates
elements = [1 2; 1 3; 1 4]; % Element connectivity

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);

% Initialize global stiffness and mass matrices
K = zeros(2*numNodes);
M = zeros(2*numNodes);
F = [0; -50e3];

% Assemble global stiffness and mass matrices
for i = 1:numElements
    n1 = elements(i, 1);
    n2 = elements(i, 2);

    % Element length and direction cosines
    x1 = nodes(n1, 1); y1 = nodes(n1, 2);
    x2 = nodes(n2, 1); y2 = nodes(n2, 2);
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    c = (x2 - x1) / L; % Cosine
    s = (y2 - y1) / L; % Sine

    C = (x2 - x1) / L; % Cosine
    S = (y2 - y1) / L; % Sine

    % Element stiffness matrix
    k = (E * A / L) * [c^2, c*s, -c^2, -c*s; 
                             c*s, s^2, -c*s, -s^2; 
                             -c^2, -c*s, c^2, c*s; 
                             -c*s, -s^2, c*s, s^2];

    % Mass per unit length
    %mass_per_length = rho * A;
    m = (rho * A * L / 2) * [1, 0,0,0; 
                                   0, 1,0,0; 
                                   0,0,1, 0; 
                                   0, 0,0,1];

    % Assemble into global matrices
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    K(dof, dof) = K(dof, dof) + k;
    M(dof, dof) = M(dof, dof) + m;
end

% Apply boundary conditions (constrain nodes 2, 3, and 4)
fixedDofs = [2*2-1 2*2 2*3-1 2*3 2*4-1 2*4]; % Constrained DOFs
freeDofs = setdiff(1:2*numNodes, fixedDofs);

% Solve the eigenvalue problem for natural frequencies
K = K(freeDofs, freeDofs);

M = M(freeDofs, freeDofs);




dt = 0.2 / (height(abaqus) - 1);
t_total = 0.2;
t = 0:dt:t_total;
n_steps = length(t);

u = zeros(2,n_steps);
v = zeros(2,1);
a = inv(M) * (F - K * u(:,1));
u(:,2) = (dt^2/2) * a;




%alpha = 0.00001;
%beta = 0.00001;0.00001:0.00001:0.001
difCounter = 1;
alpha_r = 0.00001:0.00001:0.001;
beta_r = 3.000e-05:1e-7:4.000e-05;

for j = alpha_r
    for k = beta_r
        alpha = j;
        beta = k;
        C_reduced = alpha * M + beta * K;
        for i = 2:n_steps-1
            u(:,i+1) = inv(M) * (F - K * u(:,i) - C_reduced * (u(:,i) - u(:,i-1)) / dt) * dt^2 + 2*u(:,i) - u(:,i-1);
        end
        dif(difCounter) = sum(abs(u(1,:) - (abaqus.x_1_1(:,1))'));
        difCounter = difCounter + 1;
    end
end
%dif = abs(dif);


smallest = find(dif == min(dif));


beta_s = mod(smallest, 100);   % 53
alpha_s = (smallest - beta_s) / 100; % y = 9

if alpha_s == 0
    alpha_s = 1;
end

alpha = alpha_r(alpha_s);
beta = beta_r(beta_s);

C_reduced = alpha * M + beta * K;
for i = 2:n_steps-1
    u(:,i+1) = inv(M) * (F - K * u(:,i) - C_reduced * (u(:,i) - u(:,i-1)) / dt) * dt^2 + 2*u(:,i) - u(:,i-1);
end
abaqus = readtable('abaqus.txt');

figure;
subplot(2,1,1);
plot(t,u(1,:));
hold on
plot(abaqus.X,abaqus.x_1_1)
xlabel('Time (s)')
ylabel('Displacement in X (m)')
title('Dynamic Response at node 1 in X direction')

subplot(2,1,2);
plot(t,u(2,:));
hold on
plot(abaqus.X,abaqus.x_1_1_1)
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 1 in Y direction')




