function [ E, fii, nii ] = DynamicsNE( T_i_j, dQ, ddQ, Ci, m, dv00, j_type, ...
                                            fE0, nE0, Ici, w00, dw00 )
%DYNAMICS_NE Dynamics using Newton-Euler algorithm
% Inputs
%   T_i_j   {n+1x1} cell of incremental transformation matrices T_i_j,
%           starting at the limb base frame to the final n-th actuated
%           frame PLUS the final frame (E) where external load is applied. 
%           n also represents the number of mass centres. Matrices can be
%           symbolic.
%   dQ      [1xn] vector of the joint velocity for the given T_i_j
%           frame. Accepts symbolic variables.
%   dQQ     [1xn] vector of the joint acceleration for the given T_i_j
%           frame. Accepts symbolic variables.
%   Ci      [3xn] matrix of location of mass centres, relative to its
%           parent link frame. Each column i is the i-th mass centre
%           location relative to frame i in x y z.
%   m       [1xn] vector representing mass of each mass centre Ci (in kg)
%   dv00    [3x1] vector of acceleration at the base frame, measured at the
%           base frame (i.e. gravity)
%   j_type	(Optional) [1xn] vector indicating joint type for a given
%           frame T_i_j
%           0 = revolute joint
%           1 = prismatic joint
%           If omitted, all revolute joints are assumed.
%   fE0     (Optional) [3x1] vector of external forces at the n+1 (E) frame,
%           measured at the base frame
%   nE0     (Optional) [3x1] vector of external torques at the n+1 (E) 
%           frame, measured at the base frame
%   Ici     (Optional) {n+1} cell of [3x3] inertia tensors for each mass Ci
%   w00     (Optional) [3x1] vector of angular velocity at the base frame,
%           measured at the base frame
%   dw00    (Optional) [3x1] vector of angular acceleration at the base 
%           frame, measured at the base frame
% Output
%   E       [nx1] vector of effort (torque, force) experienced at each
%           actuator
%   fii     [3xn] vector of force exerted on each link, where each column i
%           represents force on link i, measured at frame i
%   nii     [3xn] vector of torque exerted on each link, where each column i
%           represents torque on link i, measured at frame i

nl = length(T_i_j)-1;   % Number of links
% Unit vectors
X = [1 0 0]';
Y = [0 1 0]';
Z = [0 0 1]';

if nargin < 7       % Assume all revolute joints if not given
    j_type = zeros(nl,1);
end
if nargin < 8       % Assume all zero external force if not given
    fE0 = [0 0 0]';
end
if nargin < 9       % Assume all zero external torque if not given
    nE0 = [0 0 0]';
end
if nargin < 10       % Assume point masses if no inertia tensors given
    Ici = arrayfun(@(x) zeros(3), 1:nl, 'Uni', 0)';
end
if nargin < 11      % Assume all zero angular vel/acc at base if not given
    w00 = [0 0 0]';
    dw00 = [0 0 0]';
end

% Initialise all matrices
w = NaN(3,nl,class(dQ));    % Angular velocity
dw = NaN(3,nl,class(dQ));   % Angular acceleration
dv = NaN(3,nl,class(dQ));   % Linear acceleration
dvc = NaN(3,nl,class(dQ));  % Linear acceleration of mass centre
F = NaN(3,nl,class(dQ));    % Force on link exerted by link CoM
N = NaN(3,nl,class(dQ));    % Torque on link exerted by link CoM
E = NaN(nl,1,class(dQ));    % Actuator effort

% ============= YOUR CODE (START) ================
%% Outward iterations
% Use the below commented code (3 lines) if you want
w0 = w00;                   % Angular velocity of base frame
dw0 = dw00;                 % Angular acceleration of base frame
dv0 = dv00;                 % Linear acceleration of base frame

for i = 1:nl
    Tij = T_i_j{i};         % Transformation matrix
    Rij = Tij(1:3,1:3);
    Rji = transpose(Rij);
    Pij = Tij(1:3,4);

    if j_type(i)    % Prismatic joint
        w(:,i) = Rji*w0;
        dw(:,i) = Rji*dw0;
        dv(:,i) = Rji*(cross(dw0,Pij)+cross(w0,cross(w0,Pij))+dv0)+2*cross(w(:,i),dQ(i)*Z)+ddQ(i)*Z;
    else            % Revolute joint
        w(:,i) = Rji*w0+dQ(i)*Z;
        dw(:,i) = Rji*dw0+Rji*cross(w0,dQ(i)*Z)+ddQ(i)*Z;
        dv(:,i) = Rji*(cross(dw0,Pij)+cross(w0,cross(w0,Pij))+dv0);
    end
    
    % calculate acceleration at COM, and inertial force and moment at COM.
    dvc(:,i) = cross(dw(:,i),Ci(:,i))+cross(w(:,i),cross(w(:,i),Ci(:,i)))+dv(:,i);
    F(:,i) = m(i)*dvc(:,i);
    N(:,i) = Ici{i}*dw(:,i)+cross(w(:,i),Ici{i}*w(:,i));
     
    % Use the below commented code (3 lines) if you want
    w0 = w(:,i);            % Angular velocity of this frame to use in next iteration
    dw0 = dw(:,i);          % Angular acceleration of this frame to use in next iteration
    dv0 = dv(:,i);          % Linear acceleration of this frame to use in next iteration
end
w
dw
dv
dvc
F
N

%% Inward iterations
%
fii = NaN(3,nl,class(dQ));    % Force at joint
nii = NaN(3,nl,class(dQ));    % Torque at joint
% Use the below commented code (2 lines) if you want
fj = fE0;
nj = nE0;

for i = nl:-1:1
    % transformation matrices:
    Tij = T_i_j{i+1};         % Transformation matrix
    Rij = Tij(1:3,1:3);
    Pij = Tij(1:3,4);
    
    
    fii(:,i) = Rij*fj+F(:,i);
    nii(:,i) = N(:,i)+Rij*nj+cross(Ci(:,i),F(:,i))+cross(Pij,Rij*fj);
    
    % Use the below commented code (2 lines) if you want
    fj = fii(:,i);
    nj = nii(:,i);
    
    if j_type(i)    % Prismatic joint
        E(i) = fii(:,i).'*Z;
    else            % Revolute joint
        E(i) = nii(:,i).'*Z;
    end
end

end

