function [E,A1,A2t5,L1,L2,L3,L4,L5,Theta1,Theta2,Theta3,Theta4,Theta5,k1,k2,k3,k4,k5] = Setup()

    % Material properties and cross-sectional area
    E = 209e9;   % Young's Modulus (Pa)
    A1 = 1e-3;    % Cross-sectional area (m^2)
    A2t5 = 1e-4;    % Cross-sectional area (m^2)

    %Node Coordinates

    N1 = [2060,1000];
    N2 = [1700,0];
    N3 = [1250,225];
    N4 = [0,850];
    N5 = [0,1750];
    
    % Lengths of the truss elements
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
    
    % Element stiffness matrices for the three elements
    k1 = PlaneTrussElementStiffness(E, A1, L1, Theta1);   % Stiffness of element 1
    k2 = PlaneTrussElementStiffness(E, A2t5, L2, Theta2); % Stiffness of element 2
    k3 = PlaneTrussElementStiffness(E, A2t5, L3, Theta3); % Stiffness of element 3
    k4 = PlaneTrussElementStiffness(E, A2t5, L4, Theta4); % Stiffness of element 4
    k5 = PlaneTrussElementStiffness(E, A2t5, L5, Theta5); % Stiffness of element 5
end