% Declare the symbolic variables
syms th1 th2 th3 th4 th5 th6 d_th1 d_th2 d_th3 d_th4 d_th5 d_th6 d4 a2

% Alternatively, can define the values to check a solution
% d2 = 0.089159;
% a3 = -0.425;
% a4 = -0.39225;
% d5 = 0.10915;
% d6 = 0.09465;
% d7 = 0.0823;

% th1 = 0;
% th2 = 0;
% th3 = 0;
% th4 = 0;
% th5 = 0;
% th6 = 0;

% dh table describing transformation matrices from 0 to ef.
dh = [0 0 0 th1;
      0 pi/2 0 th2;
      a2 0 0 th3;
      0 pi/2 d4 th4;
      0 -pi/2 0 th5;
      0 pi/2 0 th6];


% direct kinematics converts to incremental and total transformation
% matrices.
[Tij, T0n] = dk(dh);

% declare dQ, with 0 as final variable for ef frame transformation
dQ = [d_th1, d_th2, d_th3, d_th4, d_th5, d_th6, 0];

% declare type of joints: 0 revolute, 1 prismatic.
j_type = [0 0 0 0 0 0];

% calculate the end effector velocity using propagation method.
[ vii, wii, v0i, w0i ] = VelProp(Tij, dQ, j_type);

vef = simplify(v0i(:,6));
wef = simplify(w0i(:,6));

% Given vef, find the jacobian in frame 0
J_pos = simplify(jacobian(vef,dQ(1:6)))
J0_w = simplify(jacobian(wef,dQ(1:6)))