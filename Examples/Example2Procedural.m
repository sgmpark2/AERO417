% Truss Analysis and Visualisation in MATLAB
% This script analyses a 2D truss structure, calculating node displacements,
% reaction forces, and elemental stresses. It also visualises the original
% and deformed truss.

% Node coordinates [x, y] (in meters)
node = [0 0;  % Node 1: (0, 0)
        0 3;  % Node 2: (0, 3)
        3 3;  % Node 3: (3, 3)
        3 0]; % Node 4: (3, 0)

% Connectivity matrix: Each row represents an element connecting two nodes
conn = [1 2;  % Element 1: Connects Node 1 and Node 2
        1 3;  % Element 2: Connects Node 1 and Node 3
        1 4]; % Element 3: Connects Node 1 and Node 4

% Prescribed degrees of freedom (isolated nodes)
isol = [1 2]; % Node 1's x and y displacements are free to move.

% Material and cross-sectional properties
A = 6e-4;    % Cross-sectional area in m^2
E = 200e9;   % Young's modulus in Pa

% Number of nodes and elements
nn = size(node, 1);  % Number of nodes
ndof = 2 * nn;       % Total degrees of freedom (2 per node)
ne = size(conn, 1);  % Number of elements

% Plot original truss structure
figure; hold on; % Create a new figure and hold the plot
for e = 1:ne
    n1 = conn(e, 1);
    n2 = conn(e, 2);
    
    x1 = node(n1, 1); y1 = node(n1, 2);
    x2 = node(n2, 1); y2 = node(n2, 2);
    
    X = [x1, x2];
    Y = [y1, y2];
    line(X, Y, 'linewidth', 3); % Plot the element as a line between nodes
end

% External force vector (F)
F = zeros(2 * nn, 1); % Initialize force vector
F(1) = 0;             % No force on Node 1 in x-direction
F(2) = -50000;        % Apply a downward force of 50 kN at Node 1 in y-direction

% Global stiffness matrix (K) initialisation
K = zeros(ndof, ndof); % Initialise global stiffness matrix
d = zeros(ndof, 1);    % Initialise displacement vector

% Assembly of global stiffness matrix (K)
for e = 1:ne
    n1 = conn(e, 1);
    n2 = conn(e, 2);
    
    x1 = node(n1, 1); y1 = node(n1, 2);
    x2 = node(n2, 1); y2 = node(n2, 2);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2); % Length of element
    c = (x2 - x1) / L; % Cosine of angle
    s = (y2 - y1) / L; % Sine of angle
    
    % Element stiffness matrix in local coordinates
    ke = (A * E / L) * [c*c, c*s, -c*c, -c*s;
                        c*s, s*s, -c*s, -s*s;
                        -c*c, -c*s, c*c, c*s;
                        -c*s, -s*s, c*s, s*s];
    
    % Global DOF indices for the current element
    sctr = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    
    % Assemble global stiffness matrix by adding element stiffness matrices
    K(sctr, sctr) = K(sctr, sctr) + ke;
end

% Calculate unknown displacements (d) for isolated nodes
d(isol) = K(isol, isol) \ F(isol);

% Display node displacements
fprintf('\n---------- Nodal Displacements ----------');
fprintf('\nNode  X-Displacement (m)  Y-Displacement (m)\n');
for i = 1:nn
    fprintf('%4d  %16.3e  %16.3e\n', i, d(2*i-1), d(2*i));
end

% Calculate reaction forces (f)
f = K * d;

% Display reaction forces at all nodes
fprintf('\n---------- Reactions ----------');
fprintf('\nNode  X-Reaction (N)  Y-Reaction (N)\n');
for i = 1:nn
    fprintf('%4d  %16.3e  %16.3e\n', i, f(2*i-1), f(2*i));
end

% Plot deformed truss structure
for e = 1:ne
    n1 = conn(e, 1);
    n2 = conn(e, 2);
    
    % Original node positions
    x1 = node(n1, 1); y1 = node(n1, 2);
    x2 = node(n2, 1); y2 = node(n2, 2);
    
    % Deformed node positions
    x1new = x1 + d(2*n1-1); y1new = y1 + d(2*n1);
    x2new = x2 + d(2*n2-1); y2new = y2 + d(2*n2);
    
    % Plot the deformed truss element
    X = [x1new, x2new];
    Y = [y1new, y2new];
    line(X, Y, 'linewidth', 3, 'color', [1 0 1]); % Magenta line for deformed shape
end

% Calculate and display elemental strain and stress
fprintf('\n---------- Elemental Strain & Stress ----------');
fprintf('\nElem  Strain          Stress (Pa)\n');
for e = 1:ne
    n1 = conn(e, 1);
    n2 = conn(e, 2);
    
    x1 = node(n1, 1); y1 = node(n1, 2);
    x2 = node(n2, 1); y2 = node(n2, 2);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2); % Length of element
    c = (x2 - x1) / L; % Cosine of angle
    s = (y2 - y1) / L; % Sine of angle
    
    % Strain-displacement matrix (B)
    B = (1 / L) * [-c -s c s];
    
    % Global DOF indices for the current element
    sctr = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    
    % Calculate strain and stress for the element
    Strain = B * d(sctr);
    Stress = E * Strain;
    
    % Display elemental strain and stress
    fprintf('%4d  %16.3e  %16.3e\n', e, Strain, Stress);
end

% End of script