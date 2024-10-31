% Truss FEM Analysis for Free Vibration with Lumped Mass Matrix (3 Elements)

% Define truss properties
E = 200e9; % Young's modulus in Pascals
A = 6e-4; % Cross-sectional area in m^2
L = 3; % Length of each truss member in meters
rho = 7800; % Density in kg/m^3

% Define nodal coordinates and connectivity for 3 elements
nodes = [0 0; 0 L; L L; L 0]; % Coordinates of the nodes
elements = [1 2; 1 3; 1 4]; % Connectivity of the elements

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);

% Global mass and stiffness matrices
M = zeros(numNodes * 2); % Mass matrix
K = zeros(numNodes * 2); % Stiffness matrix

% Assemble global mass and stiffness matrices
for i = 1:numElements
    n1 = elements(i, 1);
    n2 = elements(i, 2);

    % Element length and direction cosines
    x1 = nodes(n1, 1); y1 = nodes(n1, 2);
    x2 = nodes(n2, 1); y2 = nodes(n2, 2);
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    c = (x2 - x1) / L; % Cosine
    s = (y2 - y1) / L; % Sine

    % Element stiffness matrix
    k = (E * A / L) * [c^2, c*s, -c^2, -c*s; 
                             c*s, s^2, -c*s, -s^2; 
                             -c^2, -c*s, c^2, c*s; 
                             -c*s, -s^2, c*s, s^2];

    % Lumped mass matrix for each element
    m = (rho * A * L / 6) * [2, 0,1,0; 
                                   0, 2,0,1; 
                                   1,0,2, 0; 
                                   0, 1,0,2];

    % Assemble into global matrices
    dof = [2*n1-1, 2*n1, 2*n2-1, 2*n2]; % Degrees of freedom
    K(dof, dof) = K(dof, dof) + k;
    M(dof, dof) = M(dof, dof) + m;
end

% Apply external force at node 1 (50 kN)
F = zeros(numNodes * 2, 1); % Global force vector
F(2*1) = -50e3; % Apply force in the Y direction at node 1

% Define boundary conditions (fix node 2,3,4)
fixedNodes = [2;3;4]; % Node 2,3,4 is fixed
fixedDofs = [2*fixedNodes-1; 2*fixedNodes]; % Corresponding DOFs

% Modify K and F to account for boundary conditions
K(fixedDofs, :) = 0; % Set rows to zero
K(fixedDofs, fixedDofs) = zeros(length(fixedDofs)); % Set diagonal to 1
F(fixedDofs) = 0; % Set forces to zero

% Solve eigenvalue problem for free vibration analysis
[eigenVectors, eigenValues] = eig(K, M);
frequencies = sqrt(diag(eigenValues)) ; % Convert to Hz

% Display results
disp('Natural Frequencies (rad/s):');
disp(frequencies);

% Assuming 'eigenVectors' contains the mode shapes from the previous code
% and 'nodes' and 'elements' are defined as before.

% Number of modes to visualize
numModes = 3; % Change this to visualize more modes

% Create a figure for visualization
figure;
hold on;
axis equal;
title('Mode Shapes of the Truss');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');

% Plot the original truss structure
for i = 1:numElements
    n1 = elements(i, 1);
    n2 = elements(i, 2);
    plot([nodes(n1, 1), nodes(n2, 1)], [nodes(n1, 2), nodes(n2, 2)], 'k-', 'LineWidth', 2);
end

% Plot each mode shape
for mode = 1:numModes
    % Extract the mode shape
    modeShape = eigenVectors(:, mode);

    % Scale the mode shape for visualization
    scaleFactor = 0.1; % Adjust this factor for better visualization
    scaledShape = modeShape * scaleFactor;

    % Plot the deformed shape
    for i = 1:numElements
        n1 = elements(i, 1);
        n2 = elements(i, 2);
        % Calculate the deformed positions
        x1 = nodes(n1, 1) + scaledShape(2*n1-1);
        y1 = nodes(n1, 2) + scaledShape(2*n1);
        x2 = nodes(n2, 1) + scaledShape(2*n2-1);
        y2 = nodes(n2, 2) + scaledShape(2*n2);
        plot([x1, x2], [y1, y2], 'r--', 'LineWidth', 1.5); % Deformed shape
    end
end

hold off;