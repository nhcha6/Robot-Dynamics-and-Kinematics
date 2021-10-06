function I_new = sim_transform(I, Rip)
% performs the similarity transform to change the frame an inertial tensor
% is measured in.
% inputs:
%   I: cell containing a series of inertial tensors which need to be
%   converted into a new reference frame.
%   Rip: the rotation matrix which transforms a coordinate from the current
%   reference frame of I to the desired reference frame.
    no_limbs = size(Rip);
    no_limbs = no_limbs(3);
    I_new = cell(size(Rip));
    for i = 1:no_limbs
        I_new{i} = Rip{i}*I{i}*transpose(Rip{i});
    end
end