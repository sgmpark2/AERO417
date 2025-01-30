%[E,A1,L1,L2,Theta1,Theta2,K] = Setup_2();

syms E A I L u v f1 f2 f11 f3 P



x = 90*pi/180;
C = cos(x);
S = sin(x);
y1 = round([C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S ;
   -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S]);
k1 = E*A/L* y1;

x = 45*pi/180;
C = cos(x);
S = sin(x);
y2 = 2 * [C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S ;
   -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S];
k2 = E*A/(sqrt(2)*L)* y2;

% Initialise global stiffness matrix (6x6 since we have 3 nodes with 2 DOF each)
K = zeros(6, 6);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, y1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, y2, 1, 3); % Assemble element 2 between nodes 1 and 3

fixedDofs = [ 3 4 5 6]; % Constrained DOFs
freeDofs = setdiff(1:2*3, fixedDofs);

K = K(freeDofs, freeDofs);
disp('Global Equation')
disp((A*E/(2*sqrt(2)* L)) * K )
Ge = (A*E/(2*sqrt(2)* L) * K) * [u;v];
eqn = (A*E/(2*sqrt(2)* L) * K) * [u;v] == [0;-1];
S = solve(eqn);
u = L/(A*E);
v = - L/(A*E);

x = 90*pi/180;
C = cos(x);
S = sin(x);
eqn = [f1;f2] == (A*E/ L) * [1,-1;-1,1] * [C,S,0,0;0,0,C,S] * [u;v;0;0];
S = solve(eqn);
f1 = round(S.f1);
f2 = round(S.f2);
f01 = f1;

x = 45*pi/180;
C = cos(x);
S = sin(x);
eqn = [f11;f3] == (A*E/ L) * [1,-1;-1,1] * [C,S,0,0;0,0,C,S] * [u;v;0;0];
S = solve(eqn);
f1 = round(S.f11);
f3 = round(S.f3);
f02 = f3;

disp('axial force in each member')
disp(f01)
disp(f02)

F01 = f01 * P;
F02 = f02 * P;

x = 90*pi/180;
C = cos(x);
S = sin(x);
y1 = [C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S ;
   -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S];
y11 = [S*S -S*C -S*C S*C ; -S*C C*C S*C -C*C ;
   -S*S S*C S*S -S*C ; S*C -C*C -S*C C*C];



x = 45*pi/180;
C = cos(x);
S = sin(x);
y2 = 2 * [C*C C*S -C*C -C*S ; C*S S*S -C*S -S*S ;
   -C*C -C*S C*C C*S ; -C*S -S*S C*S S*S];
y22 = 2* [S*S -S*C -S*C S*C ; -S*C C*C S*C -C*C ;
   -S*S S*C S*S -S*C ; S*C -C*C -S*C C*C];



% Initialise global stiffness matrix (6x6 since we have 3 nodes with 2 DOF each)
K = zeros(6, 6);
K1 = zeros(6,6);

% Assembly of global stiffness matrix from element stiffness matrices
K = PlaneTrussAssemble(K, y1, 1, 2); % Assemble element 1 between nodes 1 and 2
K = PlaneTrussAssemble(K, y2, 1, 3); % Assemble element 2 between nodes 1 and 3

K1 = PlaneTrussAssemble(K1, y11, 1, 2); % Assemble element 1 between nodes 1 and 2
K1 = PlaneTrussAssemble(K1, y22, 1, 3); % Assemble element 2 between nodes 1 and 3


fixedDofs = [ 3 4 5 6]; % Constrained DOFs
freeDofs = setdiff(1:2*3, fixedDofs);

K = K(freeDofs, freeDofs);
K1 = K1(freeDofs, freeDofs);

