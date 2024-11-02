clear, close all

u1node2 = readtable('u1node2.txt');
u1node3 = readtable('u1node3.txt');
u2node2 = readtable('u2node2.txt');
u2node3 = readtable('u2node3.txt');

figure(1)
subplot(4,1,1)
plot(u1node2.Var1,u1node2.Var2)

subplot(4,1,2)
plot(u1node3.Var1,u1node3.Var2)

subplot(4,1,3)
plot(u2node2.Var1,u2node2.Var2)

subplot(4,1,4)
plot(u2node3.Var1,u2node3.Var2)


E = 209e9;   % Young's Modulus (Pa)
A1 = 1e-3;    % Cross-sectional area (m^2)
A2t5 = 1e-4;    % Cross-sectional area (m^2)
rho = 7800;

N1 = [2060,1000];
N2 = [1700,0];
N3 = [1250,225];
N4 = [0,850];
N5 = [0,1750];

% Define the nodes and elements
nodes = [N1; N2 ;N3; N4; N5]; % Node coordinates
elements = [1 2; 2 3; 3 4;4 5;5 3]; % Element connectivity

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);


% Assemble global stiffness and mass matrices
L1 = (PlaneTrussElementLength(N1(1),N1(2),N2(1),N2(2))) / 1e3;  % Length of elements 1 (m)
L2 = (PlaneTrussElementLength(N2(1), N2(2), N3(1), N3(2))) / 1e3; % Length of element 2 (m)
L3 = (PlaneTrussElementLength(N3(1), N3(2), N4(1), N4(2))) / 1e3; % Length of element 3 (m)
L4 = (PlaneTrussElementLength(N4(1), N4(2), N5(1), N5(2))) / 1e3; % Length of element 4 (m)
L5 = (PlaneTrussElementLength(N5(1), N5(2), N3(1), N3(2))) / 1e3; % Length of element 5 (m)

% Angles of truss elements

Theta1 = AngleCalculator(2060,1000,1700,0);  % Angle of elements 1 (degrees)
Theta2 = AngleCalculator(1700, 0, 1250, 225); % Angle of element 2 (degrees)
Theta3 = AngleCalculator(1250, 225, 0, 850); % Angle of element 3 (degrees)
Theta4 = AngleCalculator(0, 850, 0, 1750); % Angle of element 4 (degrees)
Theta5 = AngleCalculator(0, 1750, 1250, 225); % Angle of element 5 (degrees)

% Element stiffness matrices for the 5 elements
k1 = PlaneTrussElementStiffness(E, A1, L1, Theta1);   % Stiffness of element 1
k2 = PlaneTrussElementStiffness(E, A2t5, L2, Theta2); % Stiffness of element 2
k3 = PlaneTrussElementStiffness(E, A2t5, L3, Theta3); % Stiffness of element 3
k4 = PlaneTrussElementStiffness(E, A2t5, L4, Theta4); % Stiffness of element 4
k5 = PlaneTrussElementStiffness(E, A2t5, L5, Theta5); % Stiffness of element 5

% Initialise global stiffness matrix (10x10 since we have 5 nodes with 2 DOF each)
K = zeros(10, 10);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3
K = PlaneTrussAssemble(K, k3, 3, 4); % Assemble element 3 between nodes 3 and 4
K = PlaneTrussAssemble(K, k4, 4, 5); % Assemble element 4 between nodes 4 and 5
K = PlaneTrussAssemble(K, k5, 3, 5); % Assemble element 5 between nodes 3 and 5

% Element Mass matrices for the 5 elements
m1 = PlaneTrussElementMass(rho,A1,L1);
m2 = PlaneTrussElementMass(rho,A2t5,L2);
m3 = PlaneTrussElementMass(rho,A2t5,L3);
m4 = PlaneTrussElementMass(rho,A2t5,L4);
m5 = PlaneTrussElementMass(rho,A2t5,L5);

% Initialise global Mass matrix (10x10 since we have 5 nodes with 2 DOF each)
M = zeros(10, 10);

% Assembly of global stiffness matrix from element stiffness matrices
M = PlaneTrussMassAssemble(M, m1, 1, 2); % Assemble element 1 between nodes 1 and 2
M = PlaneTrussMassAssemble(M, m2, 2, 3); % Assemble element 2 between nodes 2 and 3
M = PlaneTrussMassAssemble(M, m3, 3, 4); % Assemble element 3 between nodes 3 and 4
M = PlaneTrussMassAssemble(M, m4, 4, 5); % Assemble element 4 between nodes 4 and 5
M = PlaneTrussMassAssemble(M, m5, 3, 5); % Assemble element 5 between nodes 3 and 5

% Apply boundary conditions (constrain nodes 2, 3, and 4)
fixedDofs = [1 2 7 9 10]; % Constrained DOFs
freeDofs = setdiff(1:2*numNodes, fixedDofs);

% Solve the eigenvalue problem for natural frequencies
K = K(freeDofs, freeDofs);
M = M(freeDofs, freeDofs);
F = [5000; 100000; 0;0;0];


dt = 0.0001;
t_total = 0.2;
t = 0:dt:t_total;
n_steps = length(t);

u = zeros(5,n_steps);
v = zeros(5,1);
a = inv(M) * (F - K * u(:,1));
u(:,2) = (dt^2/2) * a;




alpha = 0.0007;
beta = 0.00007;
damping = 0.5;

C_reduced = alpha * M + beta * K;

for i = 2:n_steps-1
    u(:,i+1) = inv(M) * (F - K * u(:,i) - C_reduced * (u(:,i) - u(:,i-1)) / dt) * dt^2 + 2*u(:,i) - u(:,i-1);
end

figure(2);
subplot(2,1,1);
plot(t,u(1,:));
xlabel('Time (s)')
ylabel('Displacement in X (m)')
title('Dynamic Response at node 2 in X direction')

subplot(2,1,2);
plot(t,u(2,:));
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 2 in Y direction')

figure(3);
subplot(2,1,1);
plot(t,u(3,:));
xlabel('Time (s)')
ylabel('Displacement in X (m)')
title('Dynamic Response at node 3 in X direction')

subplot(2,1,2);
plot(t,u(4,:));
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 3 in Y direction')

figure(4);
plot(t,u(5,:));
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 4 in Y direction')
