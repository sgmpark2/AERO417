% Analysis of a Fixed-Fixed Beam with Central Load and Moment
% This script calculates the deflection, reaction forces, and sketches
% the deflection curve of a fixed-fixed beam subjected to a load and moment.

% Beam and Load Parameters
L = 3;           % Length of each beam element (m)
w = 0;           % Distributed load (N/m) - not used in this case
E = 210e9;       % Young's Modulus (Pa)
I = 4e-4;        % Moment of Inertia (m^4)
F2y = -10000;    % Point load at the middle of the beam (N)
M2 = 20000;      % Moment applied at the middle of the beam (N*m)

% Element Stiffness Matrix for a Beam Element
kel = E*I/L^3 * [12 6*L -12 6*L;
                 6*L 4*L^2 -6*L 2*L^2;
                 -12 -6*L 12 -6*L;
                 6*L 2*L^2 -6*L 4*L^2];

% Element Force Vector (No distributed load applied)
fel = w*L/12 * [6; L; 6; -L];

% Connectivity Matrix
elem_con = [1 2; 2 3]; % Element connections between nodes
elem = size(elem_con, 1); % Number of elements
nodes = elem + 1; % Total number of nodes

% Initialise Global Stiffness Matrix and Force Vector
K = zeros(2*nodes); % Global stiffness matrix (6x6 for 3 nodes)
f = zeros(2*nodes, 1); % Global force vector

% Assembly of Global Stiffness Matrix
for i = 1:elem
    node11 = elem_con(i, 1) + (i-1)*1;
    node12 = elem_con(i, 1) + (i-1)*1 + 1;
    node21 = elem_con(i, 2) + (i-1)*1 + 1;
    node22 = elem_con(i, 2) + (i-1)*1 + 2;
    
    % Assemble the stiffness matrix into the global matrix
    K(node11, node11) = K(node11, node11) + kel(1, 1);
    K(node11, node12) = K(node11, node12) + kel(1, 2);
    K(node11, node21) = K(node11, node21) + kel(1, 3);
    K(node11, node22) = K(node11, node22) + kel(1, 4);
    
    K(node12, node11) = K(node12, node11) + kel(2, 1);
    K(node12, node12) = K(node12, node12) + kel(2, 2);
    K(node12, node21) = K(node12, node21) + kel(2, 3);
    K(node12, node22) = K(node12, node22) + kel(2, 4);
    
    K(node21, node11) = K(node21, node11) + kel(3, 1);
    K(node21, node12) = K(node21, node12) + kel(3, 2);
    K(node21, node21) = K(node21, node21) + kel(3, 3);
    K(node21, node22) = K(node21, node22) + kel(3, 4);
    
    K(node22, node11) = K(node22, node11) + kel(4, 1);
    K(node22, node12) = K(node22, node12) + kel(4, 2);
    K(node22, node21) = K(node22, node21) + kel(4, 3);
    K(node22, node22) = K(node22, node22) + kel(4, 4);
end

% Assembly of Global Force Vector
for i = 1:elem
    node11 = elem_con(i, 1) + (i-1)*1;
    node12 = elem_con(i, 1) + (i-1)*1 + 1;
    node21 = elem_con(i, 2) + (i-1)*1 + 1;
    node22 = elem_con(i, 2) + (i-1)*1 + 2;
    
    % Assemble the force vector into the global vector
    f(node11) = f(node11) + fel(1);
    f(node12) = f(node12) + fel(2);
    f(node21) = f(node21) + fel(3);
    f(node22) = f(node22) + fel(4);
end


% Apply the external forces and moments at the middle node
F2 = [F2y; M2]; % Applied force and moment at node 2

% Reduced force vector for free DOFs
fsmall = [f(3); f(4)]; 

% Solve for displacements at free DOFs
soln = K(3:4, 3:4) \ (F2 + fsmall);

% Complete Displacement Vector (including fixed DOFs)
disp = [0; 0; soln(1); soln(2); 0; 0];

% Calculate Reaction Forces
F2 = K * disp - f;

% Plot the Deflection Curve of the Beam
figure
for i = 1:elem
    node11 = elem_con(i, 1) + (i-1)*1;
    node12 = elem_con(i, 1) + (i-1)*1 + 1;
    node21 = elem_con(i, 2) + (i-1)*1 + 1;
    node22 = elem_con(i, 2) + (i-1)*1 + 2;
    
    y = (i-1)*L:0.1:i*L; % Position along the beam
    x = y - (i-1)*L; % Local coordinate within the element
    
    % Calculate the deflection curve for the element
    discurve = disp(node11)*(1 - 3*x.^2/L^2 + 2*x.^3/L^3);
    discurve = discurve + disp(node12)*(x - 2*x.^2/L + x.^3/L^2);
    discurve = discurve + disp(node21)*(3*x.^2/L^2 - 2*x.^3/L^3);
    discurve = discurve + disp(node22)*(-x.^2/L + x.^3/L^2);
    
    % Plot the deflection curve
    plot(y, discurve, 'linewidth', 3);
    hold on;
end

% Label the plot
xlabel('Beam Length (m)');
ylabel('Beam Deflection (m)');
grid on;
title('Deflection of the Beam');

% End of script