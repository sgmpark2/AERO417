% Symbolic Analysis of a Beam with Distributed Load
% This script solves for the displacements and reactions of a beam 
% subjected to a distributed load using symbolic computation.

% Define symbolic variables
syms E I L w P

% Beam properties and distributed load:
% E  - Young's Modulus (Pa)
% I  - Moment of Inertia (m^4)
% L  - Length of the beam element (m)
% w  - Uniform distributed load (N/m)
% P  - Point load (N) applied at the end of the beam

% Calculate the element stiffness matrix for the beam element
k1 = BeamElementStiffness(E, I, L); % 4x4 stiffness matrix for one beam element

% Initialise the global stiffness matrix for a 2-node beam (4 DOF total)
K = sym(zeros(4, 4)); % 4x4 symbolic zero matrix

% Assembly of the global stiffness matrix from the element stiffness matrix
K = BeamAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2

% Extract the reduced stiffness matrix for free degrees of freedom (DOF)
k = K(3:4, 3:4); % Consider DOFs 3 and 4 (free nodes)

% External force vector (including the distributed load and point load)
f = [(-w*L/2) - P; (w*L^2)/12]; % Symbolic forces applied at the free nodes

% Solve for displacements at free nodes symbolically
d = expand(simplify(k \ f));

% Assemble the complete displacement vector (including fixed nodes)
D = [0; 0; d]; % Displacements at nodes 1 (fixed) and 2 (free)

% Calculate effective element forces using the global stiffness matrix
Fe = K * D; % Symbolic expression for element forces

% Calculate equivalent nodal forces due to the distributed load
F0 = [(-w*L/2); -(w*L^2)/12; (-w*L/2); (w*L^2)/12]; % Symbolic fixed-end forces

% Calculate the final correct global nodal forces and moments
F = expand(simplify(Fe - F0)); % Symbolic reaction forces

% Display the results
disp('---------- Symbolic Nodal Displacements ----------');
pretty(D);

disp('---------- Symbolic Reaction Forces ----------');
pretty(F);

% End of script