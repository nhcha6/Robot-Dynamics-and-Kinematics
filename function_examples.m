% example code for use of functions

%% xyz_rotation(rotation, order)
% rotate a vector thx about x, thz about z, thy about y (fixed angle
% representation)

syms thx thz thy
rotation = [thx, thy, thz];
order = ['x', 'z', 'y'];
xyz_rotation(rotation, order)

% a frame is rotated thx about x, thz about z and thy about y. Find the
% rotation matrix which takes a point in the rotated frame and expresses it
% in the original frame (euler angles)

syms thx thz thy
rotation = [thx, thy, thz];
% reverse order as euler angles.
order = ['y', 'z', 'x'];
xyz_rotation(rotation, order)

%% dh_matrix(dh_variables)
% calculate the transformation matrix corresponding to the following entry
% in the dh table: ai-1 = 0, alpha = 90deg, di = l1+l2, theta = th1

syms l1 l2 th1
% declare single entry in dh table.
dh_variables = [0, pi/2, l1+l2, th1];
dh_matrix(dh_variables)

%% dk(dh_table)
% given the dh table, calculate the direct kinematics of the system.

syms th1 th2 th3 l1 l2 d3
% declare dh table
dh_table = [0, 0, 0, th1;
            l1, pi/2, 0, th2;
            l2, 0, 0, th3];
% function returns incremental transformation matrices and total
% transformation matrix.
[Tij, T0n] = dk(dh_table)

%% ik_eq(Tij, premult, ef_pos)
% find equations for the actuated variables in terms of the position and 
% orientation of the end effector (a generic form of T0n)
% given the incremental transformation matrix, the function returns the

% system of equations from which the ik solution can be calculated. Tij
% must include final end effector matrix.
eqs = ik_eq(Tij,1,1)

% given the position of the end effector in the final frame is [1 0 0 1]', 
% find equations for the actuated variables in terms of the position
% of the end effector in the base frame

% given the incremental transformation matrix, the function returns the
% system of equations from which the ik solution can be calculated. Tij and
% the end effector position in final frame input to match.
eqs = ik_eq(Tij,1,[1 0 0 1]')

%% Confirm IK Solution
% having found the IK solution from the following transformation matrix:
syms th1 d2 th3 l1 l3
T0n = [sin(th1 + th3), cos(th1 + th3), 0, l3*sin(th1 + th3) + l1*cos(th1) + d2*sin(th1);
       -cos(th1 + th3), sin(th1 + th3), 0, l1*sin(th1) - d2*cos(th1) - l3*cos(th1 + th3);
       0, 0, 1, 0;
       0, 0, 0, 1];
   
% convert T0n to a function
T0nfunc = matlabFunction(T0n);

% select valid test values for actuated variabeles and define T0n test
% test  = [th1_test, th3_test, d2_test]
test = [1;0.5;1]
l1 = 1;
l2 = 1;
% input test variables to function in correct order.
T0nt = vpa(T0nfunc(test(3),l1,l3,test(1),test(2)),3);

% declare variables needed for ik calc from t0nt
r11 = T0nt(1,1);
r21 = T0nt(2,1);
r31 = T0nt(3,1);
r12 = T0nt(1,2);
r22 = T0nt(2,2);
r32 = T0nt(3,2);
r13 = T0nt(1,3);
r23 = T0nt(2,3);
r33 = T0nt(3,3);
x = T0nt(1,4);
y = T0nt(2,4);
z = T0nt(3,4);

% calculate ik solutions and check they are the same as the test variables
% CHANGE IK SOLUTION
ik_sol = zeros(3,2);
a = x-l3*r11;
b = y+l3*r12;
% solution of first variable.
ik_sol(1,1) = atan2(b,a)+atan2(sqrt(a^2+b^2-l1^2),l1);
ik_sol(1,2) = atan2(b,a)-atan2(sqrt(a^2+b^2-l1^2),l1);
% solution of proceeding variables in for loop for multiple solutions.
for i = 1:2
    ik_sol(2,i) = atan2(r11,r12)-ik_sol(1,i);
    ik_sol(3,i) = (x-l3*r11-l1*cos(ik_sol(1,i)))/sin(ik_sol(1,i));
end
% print solution and compare to test.
ik_sol
%% VelProp( T_i_j, dQ, j_type, w00, v00 )
% find the angular and translational velocity of the end effector of a 
% system with supplied DH table and actuated variables.
syms th1 th2 th3 d_th1 d_th2 d_th3 l1 l2 l3

% dh table describing transformation matrices from 0 to ef.
dh = [0 0 0 th1;
      l1 pi/2 0 th2;
      l2 0 0 th3;
      l3 0 0 0];
% direct kinematics converts to incremental and total transformation
% matrices.
[Tij, T0n] = dk(dh);

% declare dQ, with 0 as final variable for ef frame transformation
dQ = [d_th1, d_th2, d_th3 0];

% declare type of joints: 0 revolute, 1 prismatic.
j_type = [0 0 0 0];
[ vii, wii, v0i, w0i ] = VelProp(Tij, dQ, j_type);

vef = simplify(v0i(:,4))
wef = simplify(w0i(:,4))

%% JACOBIAN
% Given vef, find the jacobian in frame 0 and frame 4
J0 = simplify(jacobian(vef,dQ(1:3)))
R40 = transpose(T0n(1:3,1:3));
J4 = simplify(R40*J0)
% solve for singularities
detJ0 = det(J0)
simplify(detJ0 == 0)

%% VELOCITY USING TIME DIFFERENTIATION
syms th1(t) th2(t) th4(t) d2 d3(t) l1 l2 l4
T05 = [cos(th1(t)+th2(t)-th4(t)), sin(th1(t)+th2(t)-th4(t)), 0, l1*cos(th1(t))+l2*cos(th1(t)+th2(t));
       sin(th1(t)+th2(t)-th4(t)), -cos(th1(t)+th2(t)-th4(t)), 0, l1*sin(th1(t))+l2*sin(th1(t)+th2(t));
       0, 0, -1, d2-d3(t)-l4;
       0, 0, 0, 1];
p05 = T05*[0 0 0 1]';
vef = diff(p05,t)

R05 = T05(1:3,1:3);
% Angular velocity is related to R as defined below.
S = simplify(diff(R05,t)*transpose(R05));
w05 = [S(3,2);S(1,3);S(2,1)]

%% lagrange(dh, pcom_i, Rip, m, I, g_vec, Q, dQ, ddQ, j_type)
% Find the dynamic equations of a system given:
%   - the dh table
%   - the location of the com of link i in frame i
%   - the rotation matrix transforming from the principal frame of link i
%     to frame i.
%   - masses of links
%   - intertial tensors of link in principal frame
%   - the gravity vector in the base frame
%   - the actuated variables and their single and double time derivative
%     variables.

% declare required symbolic variables.
syms th1(t) dth1(t) ddth1(t) d2(t) dd2(t) ddd2(t) l1 l2 m1 m2 Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 g

% dh table is known. End effector frame only included if a point mass is
% present.
dh = [0 0 0 th1(t);
      0 pi/2 d2(t) 0];
  
% from dh table also know joint type (0 revolute, 1 prismatic). If a mass
% at the end effector was added we would have to add an additional entry,
% which could be arbitrary if no joint exists.
j_type = [0 1];


% define location of com i in frame i
pcom_i = [0, 0;
          -l1 0;
          0 0];
      
% declare rotation matrices which map from principal frame i to frame i
Rip = cell([1 1 2]);
Rip{1} = [1 0 0;
          0 0 -1;
          0 1 0];
Rip{2} = eye(3);

% masses are known
m = [m1 m2];

% declare principal frame inertial tensors in cell
I = cell([1 1 2]);
I{1} = [Ixx1 0 0;
        0 Iyy1 0;
        0 0 Izz1];
I{2} = [Ixx2 0 0;
    0 Iyy2 0;
    0 0 Izz2];

% declare gravity vector in base frame.
g_vec = [0; -g; 0];

% declare actuated variables and associated time derivatives
Q = [th1(t), d2(t)];
dQ = [dth1(t), dd2(t)];
ddQ = [ddth1(t), ddd2(t)];

% call function, returns lagrange equation, kinetic energy, potential
% energy and dynamic equations.
[L, K, V, tau] = lagrange(dh, pcom_i, Rip, m, I, g_vec, Q, dQ, ddQ, j_type)

% convert to state space
syms th1 d2 dth1 dd2 ddth1 ddd2
Q_ = [th1, d2];
dQ_ = [dth1, dd2];
ddQ_ = [ddth1, ddd2];
for i = 1:numel(Q)
    tau = subs(tau, Q(i), Q_(i));
    tau = subs(tau, dQ(i), dQ_(i));
    tau = subs(tau, ddQ(i), ddQ_(i));
end
M = jacobian(tau,transpose(ddQ_))
G = jacobian(tau,g)*g
V = simplify(expand(tau - M*transpose(ddQ_) - G))

%% cubic_spline(pos1, t_int)
% find the cubic spline coefficients for the following joints space
% trajectory: th1 [0s:0 2s:135 4s:90], th2 [0s:0 2s:-135 4s:-45]

% declare position vectors for th1 and th2
pos1 = [0 135 90];
pos2 = [0 -135 -45];
% declare time vector
t_int = [0 2 4];

% call cubic spline function to return coefficients and other relevant
% variables.
[ p1, t1, cfs1, A1, b1 ] = cubic_spline(pos1, t_int);
[ p2, t2, cfs2, A2, b2 ] = cubic_spline(pos2, t_int);

%% DynamicsNE(Tij, dQ, ddQ, Ci, m, dv00, j_type,fE0, nE0, Ici, w00, dw00 )
% Find the dynamic equations of a system given:
%   - the dh table
%   - the location of the com of link i in frame i
%   - the rotation matrix transforming from the principal frame of link i
%     to frame i.
%   - masses of links
%   - intertial tensors of link in principal frame
%   - the angular/translational accelerational and velocity of the
%   base frame
%   - the external force and moment being applied to the end effector as
%   measured in the end effector frame.
%   - the actuated variables and their single and double time derivative
%     variables.

syms th1 d2 dth1 dd2 ddth1 ddd2 m1 m2 l1 g Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2

% dh table is converted to incremental transformation matrices
dh = [0 0 0 th1;
      0 pi/2 d2 0;
      0 0 0 0];
[Tij, T0n] = dk(dh);

% define velocity and acceleration of actuated variables
dQ = [dth1 dd2];
ddQ = [ddth1 ddd2];

% location of com is known
Ci = [0, 0;
      -l1 0;
       0 0];
% masses are known
m = [m1 m2];

% define acceleration in base frame (up relative to falling
% object).
dv00 = [0; g; 0];

% define types of joints (0 revolute, 1 prismatic)
j_type = [0 1];

% principal inertial tensor is known
Ip = cell([1 1 2]);
Ip{1} = [Ixx1 0 0;
        0 Iyy1 0;
        0 0 Izz1];
Ip{2} = [Ixx2 0 0;
    0 Iyy2 0;
    0 0 Izz2];

% need rotation matrix to convert from principal frame to frame i
Rip = cell([1 1 2]);
Rip{1} = [1 0 0;
          0 0 -1;
          0 1 0];
Rip{2} = eye(3);

% call similarity transform function to convert to frame i
Ici = sim_transform(Ip, Rip);

% external forces and angular acceleration and velocity of base frame are
% zero.
fE0 = [0; 0; 0];
nE0 = [0; 0; 0];
w00 = [0; 0; 0];
dw00 = [0; 0; 0];

% call function
[ E, fii, nii ] = DynamicsNE(Tij, dQ, ddQ, Ci, m, dv00, j_type,fE0, nE0, Ici, w00, dw00 );
% print effort (torque) vector.
E

% convert to state space
M = jacobian(E,transpose(ddQ))
G = jacobian(E,g)*g
V = simplify(expand(E - M*transpose(ddQ) - G))

%% converting a symbolic expression to a function
syms a b c d
T0n = [a, b;c, d];
T0n_func = matlabFunction(T0n)
T0n_func(1,2,3,4)