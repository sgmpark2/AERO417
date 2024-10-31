clc
clear

%Run the setup function
[E,A1,A2t5,L1,L2,L3,L4,L5,Theta1,Theta2,Theta3,Theta4,Theta5,k1,k2,k3,k4,k5] = Setup();

% Initialise global stiffness matrix (8x8 since we have 4 nodes with 2 DOF each)
K = zeros(10, 10);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3
K = PlaneTrussAssemble(K, k3, 3, 4); % Assemble element 3 between nodes 3 and 4
K = PlaneTrussAssemble(K, k4, 4, 5); % Assemble element 4 between nodes 4 and 5
K = PlaneTrussAssemble(K, k5, 3, 5); % Assemble element 5 between nodes 3 and 5

% Extract the reduced stiffness matrix for free degrees of freedom (DOF)
k = K([3,4,5,6,8] , [3,4,5,6,8]); % Node 2 is free (not fixed)

% External force vector (applied forces)
f1 = [5000; 100000; 0;0;0]; % Force of 50 kN applied at Node 2 in the negative Y direction

% Solve for displacements at free nodes (Node 2)
d = k \ f1;

% Assemble the complete displacement vector (including fixed nodes)
D = [ 0;0; d(1);d(2);d(3);d(4);0;d(5);0;0]; % Nodes 1, 4 and 5 are fixed, hence their displacements are 0

% Compute reaction forces at all nodes
F = K * D;

% Extract displacements for each element to calculate stress
d1 = [D(3); D(4); D(1); D(2)]; % Displacement vector for element 1
d2 = [D(3); D(4); D(5); D(6)]; % Displacement vector for element 2
d3 = [D(5); D(6); D(7); D(8)]; % Displacement vector for element 3
d4 = [D(7); D(8); D(9); D(10)]; % Displacement vector for element 4
d5 = [D(5); D(6); D(9); D(10)]; % Displacement vector for element 5

% Calculate stress in each element
sigma1 = PlaneTrussElementStress(E, L1, Theta1, d1); % Stress in element 1
sigma2 = PlaneTrussElementStress(E, L2, Theta2, d2); % Stress in element 2
sigma3 = PlaneTrussElementStress(E, L3, Theta3, d3); % Stress in element 3
sigma4 = PlaneTrussElementStress(E, L4, Theta4, d4); % Stress in element 4
sigma5 = PlaneTrussElementStress(E, L5, Theta5, d5); % Stress in element 5

% Display results
fprintf('---------- Nodal Displacements ----------\n');
fprintf('Node  Displacement X (m)  Displacement Y (m)\n');
for i = 1:5
    fprintf('%4d  %16.3e  %16.3e\n', i, D(2*i-1), D(2*i));
end

fprintf('\n---------- Reaction Forces ----------\n');
fprintf('Node  Reaction Force X (N)  Reaction Force Y (N)\n');
for i = 1:5
    fprintf('%4d  %16.3e  %16.3e\n', i, F(2*i-1), F(2*i));
end

fprintf('\n---------- Element Stresses ----------\n');
fprintf('Element  Stress (Pa)\n');
fprintf('%7d  %16.3e\n', 1, sigma1);
fprintf('%7d  %16.3e\n', 2, sigma2);
fprintf('%7d  %16.3e\n', 3, sigma3);
fprintf('%7d  %16.3e\n', 4, sigma4);
fprintf('%7d  %16.3e\n', 5, sigma5);

% End of script