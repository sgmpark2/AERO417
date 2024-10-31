% Analysis of a 2D Truss Structure using FEM
% This script calculates the global stiffness matrix, node displacements,
% reaction forces, and stresses for a 2D truss structure.

% Material properties and cross-sectional area
E = 200e9;   % Young's Modulus (Pa)
A = 6e-4;    % Cross-sectional area (m^2)

% Lengths of the truss elements
L1and3 = 3;  % Length of elements 1 and 3 (m)
L2 = PlaneTrussElementLength(0, 0, 3, 3); % Length of element 2 (m)

% Angles of the truss elements with respect to the global X-axis
theta1 = 90; % Angle for element 1 (degrees)
theta2 = 45; % Angle for element 2 (degrees)
theta3 = 0;  % Angle for element 3 (degrees)

% Element stiffness matrices for the three elements
k1 = PlaneTrussElementStiffness(E, A, L1and3, theta1); % Stiffness of element 1
k2 = PlaneTrussElementStiffness(E, A, L2, theta2);     % Stiffness of element 2
k3 = PlaneTrussElementStiffness(E, A, L1and3, theta3); % Stiffness of element 3

% Initialise global stiffness matrix (8x8 since we have 4 nodes with 2 DOF each)
K = zeros(8, 8);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, k2, 1, 3); % Assemble element 2 between nodes 1 and 3
K = PlaneTrussAssemble(K, k3, 1, 4); % Assemble element 3 between nodes 1 and 4


% Extract the reduced stiffness matrix for free degrees of freedom (DOF)
k = K(1:2, 1:2); % Node 1 is free (not fixed)

% External force vector (applied forces)
f = [0; -50000]; % Force of 50 kN applied at Node 1 in the negative Y direction

% Solve for displacements at free nodes (Node 1)
d = k \ f;

% Assemble the complete displacement vector (including fixed nodes)
D = [d; 0; 0; 0; 0; 0; 0]; % Nodes 2, 3, and 4 are fixed, hence their displacements are 0

% Compute reaction forces at all nodes
F = K * D;

% Extract displacements for each element to calculate stress
d1 = [D(1); D(2); D(3); D(4)]; % Displacement vector for element 1
d2 = [D(1); D(2); D(5); D(6)]; % Displacement vector for element 2
d3 = [D(1); D(2); D(7); D(8)]; % Displacement vector for element 3

% Calculate stress in each element
sigma1 = PlaneTrussElementStress(E, L1and3, theta1, d1); % Stress in element 1
sigma2 = PlaneTrussElementStress(E, L2, theta2, d2);     % Stress in element 2
sigma3 = PlaneTrussElementStress(E, L1and3, theta3, d3); % Stress in element 3

% Display results
fprintf('---------- Nodal Displacements ----------\n');
fprintf('Node  Displacement X (m)  Displacement Y (m)\n');
for i = 1:4
    fprintf('%4d  %16.3e  %16.3e\n', i, D(2*i-1), D(2*i));
end

fprintf('\n---------- Reaction Forces ----------\n');
fprintf('Node  Reaction Force X (N)  Reaction Force Y (N)\n');
for i = 1:4
    fprintf('%4d  %16.3e  %16.3e\n', i, F(2*i-1), F(2*i));
end

fprintf('\n---------- Element Stresses ----------\n');
fprintf('Element  Stress (Pa)\n');
fprintf('%7d  %16.3e\n', 1, sigma1);
fprintf('%7d  %16.3e\n', 2, sigma2);
fprintf('%7d  %16.3e\n', 3, sigma3);

% End of script