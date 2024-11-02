% Define the properties of the truss
% Material properties and cross-sectional area
E = 209e9;   % Young's Modulus (Pa)
A1 = 1e-3;    % Cross-sectional area (m^2)
A2t5 = 1e-4;    % Cross-sectional area (m^2)
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


% Assemble global stiffness and mass matrices
L1 = (PlaneTrussElementLength(N1(1),N1(2),N2(1),N2(2))) / 1e3;  % Length of elements 1 (m)
L2 = (PlaneTrussElementLength(N2(1), N2(2), N3(1), N3(2))) / 1e3; % Length of element 2 (m)
L3 = (PlaneTrussElementLength(N3(1), N3(2), N4(1), N4(2))) / 1e3; % Length of element 3 (m)
L4 = (PlaneTrussElementLength(N4(1), N4(2), N5(1), N5(2))) / 1e3; % Length of element 4 (m)
L5 = (PlaneTrussElementLength(N5(1), N5(2), N3(1), N3(2))) / 1e3; % Length of element 5 (m)

% Angles of truss elements

Theta1 = AngleCalculator(2060,1000,1700,0);  % Angle of elements 1 (degrees)
Theta2 = AngleCalculator(1700, 0, 1250, 225); % Angle of element 2 (degrees)
Theta3 = AngleCalculator(1250, 225, 0, 850); % Angle of element 3 (degrees)
Theta4 = AngleCalculator(0, 850, 0, 1750); % Angle of element 4 (degrees)
Theta5 = AngleCalculator(0, 1750, 1250, 225); % Angle of element 5 (degrees)

% Element stiffness matrices for the 5 elements
k1 = PlaneTrussElementStiffness(E, A1, L1, Theta1);   % Stiffness of element 1
k2 = PlaneTrussElementStiffness(E, A2t5, L2, Theta2); % Stiffness of element 2
k3 = PlaneTrussElementStiffness(E, A2t5, L3, Theta3); % Stiffness of element 3
k4 = PlaneTrussElementStiffness(E, A2t5, L4, Theta4); % Stiffness of element 4
k5 = PlaneTrussElementStiffness(E, A2t5, L5, Theta5); % Stiffness of element 5

% Initialise global stiffness matrix (10x10 since we have 5 nodes with 2 DOF each)
K = zeros(10, 10);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3
K = PlaneTrussAssemble(K, k3, 3, 4); % Assemble element 3 between nodes 3 and 4
K = PlaneTrussAssemble(K, k4, 4, 5); % Assemble element 4 between nodes 4 and 5
K = PlaneTrussAssemble(K, k5, 3, 5); % Assemble element 5 between nodes 3 and 5

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
[phi, omega2] = eig(K(freeDofs, freeDofs), M(freeDofs, freeDofs));
natural_frequencies = sqrt(diag(omega2)) / (2*pi) ;


% Plot mode shapes
for i = 1:height(natural_frequencies)
    figure;
    mode_shape = phi(:, i);
    mode_shape([1,5]) = 0;
    plot((nodes(:, 1))/1e3, (nodes(:, 2))/1e3, 'ko', 'MarkerFaceColor', 'k'); hold on;
    for j = 1:numElements
        n1 = elements(j, 1);
        n2 = elements(j, 2);
        scaleFactor = 0.5; % Adjust this factor for better visualization
        scaledShape = mode_shape * scaleFactor;
        plot([(nodes(n1, 1))/1e3 (nodes(n2, 1))/1e3], [(nodes(n1, 2))/1e3 (nodes(n2, 2))/1e3], 'b-');
        % Calculate the deformed positions
        x1 = (nodes(n1, 1))/1e3 + scaledShape(n1);
        y1 = (nodes(n1, 2))/1e3 + scaledShape(n1);
        x2 = (nodes(n2, 1))/1e3 + scaledShape(n2);
        y2 = (nodes(n2, 2))/1e3 + scaledShape(n2);

        if n1 == 4
             x1 = (nodes(n1, 1))/1e3;
        elseif n2 == 4
            x2 = (nodes(n2, 1))/1e3;
        end


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



fprintf('\n---------- Natural Frequency ----------\n');
fprintf('Natural Frequency (Hz)\n');
fprintf('%7d  %16.3e\n', 1, natural_frequencies(1));
fprintf('%7d  %16.3e\n', 2, natural_frequencies(2));
fprintf('%7d  %16.3e\n', 3, natural_frequencies(3));
fprintf('%7d  %16.3e\n', 4, natural_frequencies(4));
fprintf('%7d  %16.3e\n', 5, natural_frequencies(5));
