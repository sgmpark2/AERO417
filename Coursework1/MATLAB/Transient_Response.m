clear, close all
a = readtable('noderesults.txt');
time_Abaqus = table2array(a(:,1));
u1node2 = table2array(a(:,2));
u1node3 = table2array(a(:,3));
u2node2 = table2array(a(:,5));
u2node3 = table2array(a(:,6));
u2node4 = table2array(a(:,7));

figure(1)
subplot(2,1,1)
plot(time_Abaqus,u1node2)
xlabel('Time (s)')
ylabel('Displacement in X (m)')
title('Dynamic Response at node 2 in X direction')
subplot(2,1,2)
plot(time_Abaqus,u2node2)
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 2 in Y direction')

figure(2)
subplot(2,1,1)
plot(time_Abaqus,u1node3)
xlabel('Time (s)')
ylabel('Displacement in X (m)')
title('Dynamic Response at node 3 in X direction')


subplot(2,1,2)
plot(time_Abaqus,u2node3)
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 3 in Y direction')

figure(3)
plot(time_Abaqus,u2node4)
xlabel('Time (s)')
ylabel('Displacement in Y (m)')
title('Dynamic Response at node 4 in Y direction Abaqus')
hold on

%Run the setup function
[E,A1,A2t5,L1,L2,L3,L4,L5,Theta1,Theta2,Theta3,Theta4,Theta5,K,k] = Setup();
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


dt = 1.25 / (13946 - 1);
t_total = 1.25;
t = 0:dt:t_total;
n_steps = length(t);

u = zeros(5,n_steps);
v = zeros(5,1);
a = inv(M) * (F - K * u(:,1));
u(:,2) = (dt^2/2) * a;


alpha = 0.005;
beta = 0.00002;

C_reduced = alpha * M + beta * K;

for i = 2:n_steps-1
    u(:,i+1) = inv(M) * (F - K * u(:,i) - C_reduced * (u(:,i) - u(:,i-1)) / dt) * dt^2 + 2*u(:,i) - u(:,i-1);
end


figure(4);
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

figure(5);
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

%figure(6);
plot(t,u(5,:));
%xlabel('Time (s)')
%ylabel('Displacement in Y (m)')
%title('Dynamic Response at node 4 in Y direction Matlab')
