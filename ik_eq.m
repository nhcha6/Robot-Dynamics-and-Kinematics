function eqs = ik_eq(Tij, premult, ef_pos)
% returns the set of equations from which the ik solution can be
% calculated from the incremental transformation matrices of a system. Can
% set the number of matrices to premultiply the general
% orientation/position matrix by, and if ik is to be completed to target
% position and orientation or position only.

% INPUTS:
%   Tij: incremental transformation matrices (frame i to i+1). End effector
%   transformation matrix only included if using orientation vector.
%   premult: the number of transformation matrices to be premultiplied
%   with the generic orientation/position vector.
%   ef_pos: position of the end effector in the final frame. Set to 1 if
%   we are completing IK for orientation instead of ef position.

% OUTPUT:
%   eqs: a matrix of equations from which the IK solutions are derived.

    % declare generic variables
    syms r11 r12 r13 r21 r22 r23 r31 r32 r33 x y z
            
    % number of transformation matrices given.
    no_joints = size(Tij);
    no_joints = no_joints(3);
    
    % declare initial left hand side matrix
    LHS = [r11, r12, r13, x;
          r21, r22, r23, y;
          r31, r32, r33, z;
          0, 0, 0, 1]; 
    
    % premultiply LHS with prescribed number of transformation matrices.
    for i = 1:premult
        LHS = simplify(Tij{i}^-1*LHS);
    end
    
    % produce remaining RHS matrix.
    RHS = sym(eye(4,4));
    for j = premult+1:no_joints
        RHS = simplify(RHS*Tij{j});
    end
    
    % if working with the full orientation matrix, ef_pos = 1 and LHS = RHS
    if ef_pos == 1
        eqs = LHS==RHS;
    % else ef_pos contains the position of the end effector in the final
    % frame.
    else
        eqs = LHS(1:4,4) == RHS*ef_pos;
end