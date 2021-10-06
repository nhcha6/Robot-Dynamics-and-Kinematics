function [ vii, wii, v0i, w0i ] = VelProp( T_i_j, dQ, j_type, w00, v00 )
%VELPROP Velocity propagation method for Jacobians and statics.
% Inputs
%   T_i_j   {nx1} cell of incremental transformation matrices T_i_j,
%           starting at the limb base frame. Can be generated using T_i_j
%           or fT_i_j properties of the Bioloid robot. n is the number of
%           frames including the EE frame. Accepts symbolic matrices.
%   dQ      [1xn] vector of the joint velocity for the given T_i_j
%           frame. If a frame is translation only, set = 0. Accepts
%           symbolic variables.
%   j_type  (Optional) [1xn] vector indicating joint type for a given
%           frame T_i_j
%           0 = revolute joint
%           1 = prismatic joint
%           If omitted, all revolute joints are assumed.
%   w00     (Optional) [3x1] vector of base angular speed. Default 0s.
%   v00     (Optional) [3x1] vector of base linear velocity. Default 0s.
%
% Outputs
%   vii     [3xn] matrix of linear velocity vectors, where each column i
%           represents velocity of frame i w.r.t frame i
%   wii     [3xn] matrix same as above but for angular velocity
%   v0i     [3xn] matrix of linear velocity vectors, where each column i
%           represents velocity of frame i w.r.t arm base frame {0}
%   w0i     [3xn] matrix same as above but for angular velocity

n = length(T_i_j)

if nargin < 5
    v00 = [0 0 0]';
end
if nargin < 4 || isempty(w00)
    w00 = [0 0 0]';
end
if nargin < 3 || isempty(j_type)
    j_type = zeros(n,1);
end

wii = zeros(3,n,class(dQ));
vii = zeros(3,n,class(dQ));
% wrt base frame
w0i = zeros(3,n,class(dQ));
v0i = zeros(3,n,class(dQ));

% Initialise R matrix measured at frame {0} (keep for symbolic compatiblity)
R0 = eye(3,class(T_i_j{1})); 

% Variable naming convention
% j == i + 1 

% Use the below commented code (2 lines) if you want
wi = w00;               % Angular velocity of base frame
vi = v00;               % Linear velocity of base frame

for i = 1:n
    % extract required matrices and vectors from Tij
    Tij = T_i_j{i};     % Current transformation matrix
    Rij = Tij(1:3,1:3);
    Pij = Tij(1:3,4);
    R0 = R0*Rij;        % Rotation matrix measured at frame {0} (keep for symbolic compatiblity)
    Rij = transpose(Rij);
    
    if j_type(i)        % Prismatic
        wii(:,i) = Rij*wi;
        vii(:,i) = Rij*(vi+cross(wi,Pij))+dQ(i)*[0;0;1];
    else                % Revolute    
        wii(:,i) = Rij*wi+dQ(i)*[0; 0; 1];
        vii(:,i) = Rij*(vi+cross(wi,Pij));
    end
    
    % Hint: use R0 here to change frame of reference of velocities
    v0i(:,i) = R0*vii(:,i);
    w0i(:,i) = R0*wii(:,i);
    
    
    % Use the below commented code (2 lines) if you want    
    wi = wii(:,i);      % Angular velocity of this frame to use in next iteration
    vi = vii(:,i);      % Linear velocity of this frame to use in the next iteration
end

end
