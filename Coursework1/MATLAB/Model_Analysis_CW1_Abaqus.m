Abaqus_Mode = readtable('Abaqus Mode Shape.xlsx');

%run the setup function
[E,A1,A2t5,L1,L2,L3,L4,L5,Theta1,Theta2,Theta3,Theta4,Theta5,K,k] = Setup();
rho = 7800;

%Node Geometry setup
N1 = [2060,1000];
N2 = [1700,0];
N3 = [1250,225];
N4 = [0,850];
N5 = [0,1750];

nodes = [N1; N2 ;N3; N4; N5]; % Node coordinates
elements = [1 2; 2 3; 3 4;4 5;5 3]; % Element connectivity
numNodes = size(nodes, 1);
numElements = size(elements, 1);



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
    mode_shape_Abaqus = table2array(Abaqus_Mode(:,i));
    mode_shape([1,5]) = 0;
    %collect the mode shapes from abaqu and calculated
    plot((nodes(:, 1))/1e3, (nodes(:, 2))/1e3, 'ko', 'MarkerFaceColor', 'k'); hold on;
    for j = 1:numElements
        n1 = elements(j, 1);
        n2 = elements(j, 2);
        scaleFactor = 0.5; % Adjust this factor for better visualization
        scaleFactor2 = 0.5;
        scaledShape = mode_shape * scaleFactor;
        scaledShapeAbaqus = mode_shape_Abaqus * scaleFactor2;
        %create the scaled shape of the modes
        plot([(nodes(n1, 1))/1e3 (nodes(n2, 1))/1e3], [(nodes(n1, 2))/1e3 (nodes(n2, 2))/1e3], 'b-');
        % Calculate the deformed positions MATLAB
        x1 = (nodes(n1, 1))/1e3 + scaledShape(n1);
        y1 = (nodes(n1, 2))/1e3 + scaledShape(n1);
        x2 = (nodes(n2, 1))/1e3 + scaledShape(n2);
        y2 = (nodes(n2, 2))/1e3 + scaledShape(n2);
          % Calculate the deformed positions Abaqus
        xA1 = (nodes(n1, 1))/1e3 + (scaledShapeAbaqus(n1) );
        yA1 = (nodes(n1, 2))/1e3 + (scaledShapeAbaqus(n1) );
        xA2 = (nodes(n2, 1))/1e3 + (scaledShapeAbaqus(n2) );
        yA2 = (nodes(n2, 2))/1e3 + (scaledShapeAbaqus(n2) );


        %for the pinned nodes do not allow movement
        if n1 == 4
             x1 = (nodes(n1, 1))/1e3;
             xA1 = (nodes(n1, 1))/1e3;
        elseif n2 == 4
            x2 = (nodes(n2, 1))/1e3;
            xA2 = (nodes(n2, 1))/1e3;
        end



        plot([x1, x2], [y1, y2], 'r--', 'LineWidth', 1.5); % Deformed shape
        plot([xA1, xA2], [yA1, yA2], 'g--', 'LineWidth', 1.5); % Deformed shape
    
    end
    %titles for graphs and labels
    title(['Mode Shape for Frequency: ' num2str(natural_frequencies(i)) ' rad/s']);
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    axis equal;
    grid on;
end

%display results
fprintf('\n---------- Natural Frequency ----------\n');
fprintf('Natural Frequency (Hz)\n');
fprintf('%7d  %16.3e\n', 1, natural_frequencies(1));
fprintf('%7d  %16.3e\n', 2, natural_frequencies(2));
fprintf('%7d  %16.3e\n', 3, natural_frequencies(3));
fprintf('%7d  %16.3e\n', 4, natural_frequencies(4));
fprintf('%7d  %16.3e\n', 5, natural_frequencies(5));
