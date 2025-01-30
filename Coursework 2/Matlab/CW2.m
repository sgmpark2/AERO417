clc
clear
warning off

%requires the symbols adon for matlab

syms E L P I

% Lengths of the truss elements
L1 = PlaneTrussElementLength(0,0,0,1000)/1e3;  % Length of elements 1 (m)
L2 = PlaneTrussElementLength(0,1000,1000,1000)/1e3; % Length of element 2 (m)

% Angles of the truss elements with respect to the global X-axis
theta1 = 90; % Angle for element 1 (degrees)
theta2 = 0; % Angle for element 2 (degrees)

% Element stiffness matrices for the two elements
k1 = BeamElementStiffness(E,I,L) % Stiffness of element 1
k2 = BeamElementStiffness(E,I,L) % Stiffness of element 2


% Initialise global stiffness matrix (6x6 since we have 3 nodes with 2 DOF each)
K = zeros(6, 6);
K = sym(K);

% Assembly of global stiffness matrix from element stiffness matrices
K = BeamAssemble(K, k1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = BeamAssemble(K, k2, 2, 3); % Assemble element 2 between nodes 2 and 3

%reduce the matrix based on the fixed nodes
k = K([3,4], [3,4])

%sqrt2 cos 45 and sin 45 = 1 and -1
P = P * [1;-1];

%Nodal Displacement
d = k\P

% Element geometric matrices for the two elements
kg1 = GeometricBeamElement(P(2),L) % Stiffness of element 1
kg2 = GeometricBeamElement(P(2),L) % Stiffness of element 2


% Initialise global geometric matrix (6x6 since we have 3 nodes with 2 DOF each)
Kg = zeros(6,6);
Kg = sym(Kg);

% Assembly of global geometric matrix from element stiffness matrices
Kg = BeamAssemble(Kg, kg1, 1, 2); % Assemble element 1 between nodes 1 and 2
Kg = BeamAssemble(Kg, kg2, 2, 3); % Assemble element 2 between nodes 2 and 3


%Reduce Geometric Matrix
Kg = Kg([3,4], [3,4])

%formation of the main matrix
Kmain = k + Kg


%Find Critical Load Expressions
eq = det(Kmain) == 0

%solved values from the Critical Equation
crit = solve(eq, P(1))

%Values from Assignment Values
E_Assigned = 209e9;
I_Assigned = 5e-9;
L_Assigned = [0.7*1,0.5*1]; %Mode shapes in this problem???
P_Assigned = 1;

%Calculate Numerical Critical Load
result = zeros(2, 2);
for i = 1:2
    result([1 2],i) = subs(crit(1), [E,I,L],[E_Assigned,I_Assigned,L_Assigned(i)]);
end
result = double(result);

%include the sqrt2 in the critical load calcluation
result = result / sqrt(2);

%Calculate Numerical Displacement with respect to P
nod_Result = subs(d, [E,I,L,P(1)],[E_Assigned,I_Assigned,1,P_Assigned]);
nod_Result = double(nod_Result);


% Display results
fprintf('---------- Nodal Displacement ----------\n');
fprintf('Node  Displacement X (m)  Displacement Y (m)\n');
fprintf('%4d  %16.3e*P  %16.3e*P\n', 1, nod_Result(1), nod_Result(2));

% Display results
fprintf('\n---------- Critical Loads ----------\n');
fprintf('Mode       Critical Load (N)\n');
fprintf('%4d  %16.3e\n', 1, result(1,1));
fprintf('%4d  %16.3e\n', 2, result(1,2));
fprintf('%4d  %16.3e\n', 3, result(2,1));
fprintf('%4d  %16.3e\n', 4, result(2,2));

function y = GeometricBeamElement(F,L)
%BeamElementStiffness   This function returns the element 
%                       stiffness matrix for a beam   
%                       element with modulus of elasticity E,  
%                       moment of inertia I, and length L.
%                       The size of the element stiffness 
%                       matrix is 4 x 4.
y = F/(30*L)*[36 3*L -36 3*L ;
    3*L 4*L*L -3*L -L*L ;
   -36 -3*L 36 -3*L ;
   3*L -L*L -3*L 4*L*L];
end


