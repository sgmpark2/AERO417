function [E,A1,L1,L2,Theta1,Theta2,K] = Setup_2()

    % Material properties and cross-sectional area
    E = 209e9;    % Young's Modulus (Pa)
    A1 = 1e-3;    % Cross-sectional area (m^2)

    %Node Coordinates

    N1 = [0,0];
    N2 = [0,1000];
    N3 = [1000,1000];


    
    % Lengths of the truss elements
    L1 = (PlaneTrussElementLength(N1(1),N1(2),N2(1),N2(2))) / 1e3;  % Length of elements 1 (m)
    L2 = (PlaneTrussElementLength(N2(1), N2(2), N3(1), N3(2))) / 1e3; % Length of element 2 (m)
  
    
    % Angles of truss elements
    
    Theta1 = AngleCalculator(2060,1000,1700,0);  % Angle of elements 1 (degrees)
    Theta2 = AngleCalculator(1700, 0, 1250, 225); % Angle of element 2 (degrees)
    
    
    % Element stiffness matrices for the three elements
    k1 = PlaneTrussElementStiffness(E, A1, L1, Theta1);   % Stiffness of element 1
    k2 = PlaneTrussElementStiffness(E, A1, L2, Theta2); % Stiffness of element 2
    
    % Initialise global stiffness matrix (8x8 since we have 4 nodes with 2 DOF each)
    K = zeros(6, 6);
    
    % Assembly of global stiffness matrix from element stiffness matrices
    K = PlaneTrussAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
    K = PlaneTrussAssemble(K, k2, 1, 3); % Assemble element 2 between nodes 2 and 3
    
    
    % Extract the reduced stiffness matrix for free degrees of freedom (DOF)
    fixedDofs = [ 3 4 5 6]; % Constrained DOFs
    freeDofs = setdiff(1:2*3, fixedDofs);
    
    K = K(freeDofs, freeDofs);
end