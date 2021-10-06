function [Tij, T0n] = dk(dh_table)
% takes a dh_table and returns:
%   Tij: a series of transformation matrices from frame i to i+1
%   T0n: the transformation martrix from frame 0 to n (the final frame in
%   the dh table).
    no_joints = size(dh_table);
    no_joints = no_joints(1);

    Tij = sym(zeros(4,4,no_joints));
    T0n = sym(eye(4,4));
    for i = 1:no_joints
        Tij(:,:,i) = dh_matrix(dh_table(i,:));
        T0n = T0n*Tij(:,:,i);
    end
    Tij = num2cell(Tij,[1 2]);
end