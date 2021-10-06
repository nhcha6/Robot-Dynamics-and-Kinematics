function T = dh_matrix(dh_variables)
% takes variables from the DH table and returns the corresponding
% transformation matrix.
% INPUT
%   dh_variables: a single entry in the dh table
% OUTPUT
%   T: the transformation matrix representing the entry in the dh table.
a = dh_variables(1);
alpha = dh_variables(2);
d = dh_variables(3); 
theta = dh_variables(4);

    T = vpa([cos(theta), -sin(theta), 0, a; 
        sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -sin(alpha)*d;
        sin(theta)*sin(alpha), cos(theta)*sin(alpha), cos(alpha), cos(alpha)*d;
        0, 0, 0, 1],5);
end
