% Define the properties of the truss
E = 200e9; % Young's modulus in Pascals
A = 6e-4; % Cross-sectional area in m^2
L = 3; % Length of each element in meters
rho =7800; % Density in kg/m^3

% Define the nodes and elements
nodes = [0 0; L 0; 0 L;L L]; % Node coordinates
elements = [1 2; 1 3; 1 4]; % Element connectivity

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);

% Initialize global stiffness and mass matrices
K = zeros(2*numNodes);
M = zeros(2*numNodes);

% Assemble global stiffness and mass matrices
for i = 1:numElements
    n1 = elements(i, 1);
    n2 = elements(i, 2);

    % Element length and direction cosines
    x1 = nodes(n1, 1); y1 = nodes(n1, 2);
    x2 = nodes(n2, 1); y2 = nodes(n2, 2);
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    c = (x2 - x1) / L; % Cosine
    s = (y2 - y1) / L; % Sine

    C = (x2 - x1) / L; % Cosine
    S = (y2 - y1) / L; % Sine

    % Element stiffness matrix
    k = (E * A / L) * [c^2, c*s, -c^2, -c*s; 
                             c*s, s^2, -c*s, -s^2; 
                             -c^2, -c*s, c^2, c*s; 
                             -c*s, -s^2, c*s, s^2];

    % Mass per unit length
    %mass_per_length = rho * A;
    m = (rho * A * L / 2) * [1, 0,0,0; 
                                   0, 1,0,0; 
                                   0,0,1, 0; 
                                   0, 0,0,1];

    % Assemble into global matrices
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    K(dof, dof) = K(dof, dof) + k;
    M(dof, dof) = M(dof, dof) + m;
end

% Apply boundary conditions (constrain nodes 2, 3, and 4)
fixedDofs = [2*2-1 2*2 2*3-1 2*3 2*4-1 2*4]; % Constrained DOFs
freeDofs = setdiff(1:2*numNodes, fixedDofs);

% Solve the eigenvalue problem for natural frequencies
[phi, omega2] = eig(K(freeDofs, freeDofs), M(freeDofs, freeDofs));
%[phi, omega2] = eig(K, M);
natural_frequencies = sqrt(diag(omega2)) / (2*pi) ;

% Plot mode shapes
for i = 1:height(natural_frequencies)
    figure;
    mode_shape = phi(:, i);
    plot(nodes(:, 1), nodes(:, 2), 'ko', 'MarkerFaceColor', 'k'); hold on;
    for j = 1:numElements
        n1 = elements(j, 1);
        n2 = elements(j, 2);
        scaleFactor = 0.5; % Adjust this factor for better visualization
        scaledShape = mode_shape * scaleFactor;
        plot([nodes(n1, 1) nodes(n2, 1)], [nodes(n1, 2) nodes(n2, 2)], 'b-');
        % Calculate the deformed positions
        x1 = nodes(n1, 1) + scaledShape(2*n1-1);
        y1 = nodes(n1, 2) + scaledShape(2*n1);
        x2 = nodes(n2, 1) ;
        y2 = nodes(n2, 2);
        plot([x1, x2], [y1, y2], 'r--', 'LineWidth', 1.5); % Deformed shape
    end
    % Scale mode shape for visualization
    %{
    scale = 0.1; % Scale factor for mode shape
    quiver(nodes(:, 1), nodes(:, 2), scale*mode_shape(1:2:end), scale*mode_shape(2:2:end), 0, 'r');
    %}
    title(['Mode Shape for Frequency: ' num2str(natural_frequencies(i)) ' rad/s']);
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    axis equal;
    grid on;
end
fprintf('\n---------- Natural Frequency ----------\n');
fprintf('Natural Frequency (Hz)\n');
fprintf('%7d  %16.3e\n', 1, natural_frequencies(1));
fprintf('%7d  %16.3e\n', 2, natural_frequencies(2));