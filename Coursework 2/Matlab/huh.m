% Define symbolic variables
syms E A I L P positive
syms u1 u2 theta1 theta2 % Degrees of freedom at nodes

% Stiffness matrix of a plane frame element
K_local = E / L * [A, 0, 0, -A, 0, 0;
                   0, 12*I/L^2, 6*I/L, 0, -12*I/L^2, 6*I/L;
                   0, 6*I/L, 4*I, 0, -6*I/L, 2*I;
                   -A, 0, 0, A, 0, 0;
                   0, -12*I/L^2, -6*I/L, 0, 12*I/L^2, -6*I/L;
                   0, 6*I/L, 2*I, 0, -6*I/L, 4*I]

% Geometric stiffness matrix due to axial force
K_geo = P / L * [0, 0, 0, 0, 0, 0;
                 0, 6/L, 3, 0, -6/L, 3;
                 0, 3, 2*L, 0, -3, L;
                 0, 0, 0, 0, 0, 0;
                 0, -6/L, -3, 0, 6/L, -3;
                 0, 3, L, 0, -3, 2*L];

% Assemble global matrices (assuming single element for simplicity)
K = K_local; % Stiffness matrix
KG = K_geo; % Geometric stiffness matrix

% Eigenvalue problem
[eigenvectors, eigenvalues] = eig(K, KG);

% Critical loads
critical_loads = diag(eigenvalues); % Eigenvalues represent critical loads

% Display critical loads
disp('Critical loads in terms of E, A, I, L:');
disp(critical_loads);
