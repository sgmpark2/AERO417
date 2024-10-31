% Analysis of a Beam using FEM and Plotting Shear and Moment Diagrams
% This script calculates the global stiffness matrix, node displacements,
% reaction forces, and plots shear force and bending moment diagrams
% for a beam using the finite element method.

% Material and section properties
E = 210e9;   % Young's Modulus (Pa)
I = 4e-4;    % Moment of Inertia (m^4)
L = 3;       % Length of each beam element (m)

% Element stiffness matrices for the two beam elements
k1 = BeamElementStiffness(E, I, L); % Stiffness of element 1
k2 = BeamElementStiffness(E, I, L); % Stiffness of element 2

% Initialise global stiffness matrix (6x6 since we have 3 nodes with 2 DOF each)
K = zeros(6, 6);

% Assembly of global stiffness matrix from element stiffness matrices
K = BeamAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = BeamAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3

% Extract the reduced stiffness matrix for free degrees of freedom (DOF)
k = K(3:4, 3:4); % Node 2 is free (not fixed)

% External force vector (applied loads)
f = [-10000; 20000]; % Applied forces: -10 kN at Node 2 in the Y direction and 20 kN.m counterclockwise moment

% Solve for displacements at free nodes (Node 2)
d = k \ f;

% Assemble the complete displacement vector (including fixed nodes)
D = [0; 0; d; 0; 0]; % Nodes 1 and 3 are fixed, hence their displacements are 0

% Compute reaction forces at all nodes
F = K * D;

% Extract displacements for each element to calculate element forces
d1 = [D(1); D(2); D(3); D(4)]; % Displacement vector for element 1
d2 = [D(3); D(4); D(5); D(6)]; % Displacement vector for element 2

% Calculate element forces for each beam element
f1 = BeamElementForces(k1, d1); % Forces in element 1
f2 = BeamElementForces(k2, d2); % Forces in element 2

% Plot shear force diagrams for both elements
figure;
subplot(2,1,1); % Create a subplot for shear force diagrams
BeamElementShearDiagram(f1, L); % Shear force diagram for element 1
title('Shear Force Diagram - Element 1');
subplot(2,1,2); % Create another subplot for element 2
BeamElementShearDiagram(f2, L); % Shear force diagram for element 2
title('Shear Force Diagram - Element 2');

% Plot bending moment diagrams for both elements
figure;
subplot(2,1,1); % Create a subplot for bending moment diagrams
BeamElementMomentDiagram(f1, L); % Bending moment diagram for element 1
title('Bending Moment Diagram - Element 1');
subplot(2,1,2); % Create another subplot for element 2
BeamElementMomentDiagram(f2, L); % Bending moment diagram for element 2
title('Bending Moment Diagram - Element 2');

% Display results
fprintf('---------- Nodal Displacements ----------\n');
fprintf('Node  Displacement Y (m)  Rotation (rad)\n');
for i = 1:length(D)/2
    fprintf('%4d  %16.3e  %16.3e\n', i, D(2*i-1), D(2*i));
end

fprintf('\n---------- Reaction Forces ----------\n');
fprintf('Node  Force Y (N)  Moment (Nm)\n');
for i = 1:length(F)/2
    fprintf('%4d  %16.3e  %16.3e\n', i, F(2*i-1), F(2*i));
end

% End of script