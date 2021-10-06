function R = xyz_rotation(rotation, order)
% returns the rotation matrix produced from euler angle/fixed angle.

% rotation = [gamma, beta, alpha]

% fixed angle representation: rotates a point by gamma about x, beta about
% y and alpha about z, and returns the rotation matrix describing its new
% location within the original frame.

% euler angle representation: rotates a frame alpha about z, beta about y
% and gamma about x, and return the rotation matrix for mapping a point in
% the final frame back to the ground frame.

% order is the order of rotation in fixed angle representation (the reverse
% of the euler angle representation).

    gamma = rotation(1);
    beta = rotation(2);
    alpha = rotation(3);

    Rx = [1 0 0; 0 cos(gamma), -sin(gamma); 0, sin(gamma), cos(gamma)];
    Ry = [cos(beta) 0 sin(beta); 0 1, 0; -sin(beta), 0, cos(beta)];
    Rz = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
    R = eye(3,3);
    for i = 1:3
        if order(i)=='x'
            R = Rx*R;
        else if order(i) == 'y'
            R = Ry*R;
            else
                R = Rz*R;
            end
        end
    end
    R = vpa(R,5);
end
