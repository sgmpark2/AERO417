% Analysis of a 1D Bar Element Assembly using FEM
% This script computes the global stiffness matrix, node displacements, 
% and reaction forces for a one-dimensional bar element assembly.

% Material properties and cross-sectional areas
E1and2 = 2e11; % Young's Modulus for elements 1 and 2 (Pa)
E3 = 1e11;     % Young's Modulus for element 3 (Pa)
A1and2 = 6e-4; % Cross-sectional area for elements 1 and 2 (m^2)
A3 = 12e-4;    % Cross-sectional area for element 3 (m^2)
L = 0.6;       % Length of each bar element (m)

% Element stiffness matrices for the three elements
k1 = LinearBarElementStiffness(E1and2, A1and2, L); % Stiffness of element 1
k2 = LinearBarElementStiffness(E1and2, A1and2, L); % Stiffness of element 2
k3 = LinearBarElementStiffness(E3, A3, L);         % Stiffness of element 3


% Initialize global stiffness matrix (4x4 since we have 4 nodes)
K = zeros(4, 4);

% Assembly of global stiffness matrix from element stiffness matrices
K = LinearBarAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = LinearBarAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3
K = LinearBarAssemble(K, k3, 3, 4); % Assemble element 3 between nodes 3 and 4


% Extract the reduced stiffness matrix for free degrees of freedom (DOF)
k = K(2:3, 2:3); % Nodes 2 and 3 are free (i.e., not fixed)

% External force vector (applied forces)
f = [15000; 0]; % Force of 15 kN applied at Node 2 in the positive direction

% Solve for displacements at free nodes (Nodes 2 and 3)
d = k \ f;

% Assemble the complete displacement vector (including fixed nodes)
D = [0; d; 0]; % Nodes 1 and 4 are fixed, hence their displacements are 0

% Compute reaction forces at all nodes
F = K * D;

% Display results
fprintf('---------- Nodal Displacements ----------\n');
fprintf('Node  Displacement (m)\n');
for i = 1:length(D)
    fprintf('%4d  %16.3e\n', i, D(i));
end

fprintf('\n---------- Reaction Forces ----------\n');
fprintf('Node  Reaction Force (N)\n');
for i = 1:length(F)
    fprintf('%4d  %16.3e\n', i, F(i));
end

% End of script

