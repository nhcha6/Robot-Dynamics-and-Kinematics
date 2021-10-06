function [L, K, V, tau] = lagrange(dh, pcom_i, Rip, m, I, g_vec, Q, dQ, ddQ, j_type)
% returns the dynamic equations of a system

% INPUTS:
%   dh: the dh table describing the system. End effector frame will be
%   ingnored unless there is a point mass located there.
%   pcom_i: position of centre of mass i in frame i.
%   Rip: cell of rotation matrices mapping from principal frame i to frame i.
%   m: mass vector.
%   I: cell of principal frame inertial tensors.
%   g_vec: gravity vector in base frame.
%   Q: actuated varibles.
%   dQ: single time derivative of actuated variables.
%   ddQ: double time derivative of actuated variables.
%   j_type: type of joints, 0 revolute, 1 prismatic.

% OUTPUTS
%   L: the Lagrangian Equation (K-V)
%   K: kinetic energy equation
%   V: potential energy equation
%   tau: dynamic equations of system.
    
    % extract number of limbs (better termed number of masses due to
    % potential for point mass on end effector)
    no_limbs = size(pcom_i);
    no_limbs = no_limbs(2); 
    
    % produce incremental transformation matrices and total transformation
    % matrix from dh table.
    [Tij, T0n] = dk(dh);
    
    % differentiate actuation variables wrt time.
    diffQ = sym(zeros(size(Q)));
    for i = 1:no_limbs
        diffQ(i) = diff(Q(i));
    end
       
    % run velocity propagation to calculate vii and wii
    [ vii, wii, v0i, w0i ] = VelProp(Tij, diffQ, j_type);

    % calculate variables used to find V and K
    % velocity of com i in frame 0
    vcom_0 = sym(zeros(size(pcom_i)));
    % position of com i in frame 0
    pcom_0 = sym(zeros(size(pcom_i)));
    % angular velocity of principal frame.
    w_p = sym(zeros(size(pcom_i)));
    % initialise transformation and rotation matrix mapping from frame i to
    % base frame
    R0i = eye(3);
    T0i = eye(4);
    for i = 1:no_limbs
        % update T0i and R0i
        T0i = T0i*Tij{i};
        R0i = R0i*Tij{i}(1:3,1:3);
        % calculate velocity of com i in frame 0 based on the velocity of
        % frame i and angular velocity of frame i
        vcom_0(:,i) = v0i(:,i) + cross(w0i(:,i),R0i*pcom_i(:,i));
        % convert angular velocity of frame i to angular velocit of
        % principal frame
        w_p(:,i) = transpose(Rip{i})*wii(:,i);
        % calculate position of com i in base frame
        temp = T0i*[pcom_i(:,i);1];
        pcom_0(:,i) = temp(1:3);
    end
    % simplify results and print
    vcom_0 = simplify(vcom_0);
    w_p = simplify(w_p);
    pcom_0 = simplify(pcom_0);

    % calculate kinetic and potential energy and lagrange
    K = 0;
    V = 0;
    % iterate through each mass and calculate
    for i = 1:no_limbs
        K = K + 0.5*m(i)*transpose(vcom_0(:,i))*vcom_0(:,i);
        K = K + 0.5*transpose(w_p(:,i))*I{i}*w_p(:,i);
        V = V - m(i)*transpose(g_vec)*pcom_0(:,i);
    end
    L = simplify(K-V);
    
    % take time derivatives of lagrange to find tau
    tau = sym(zeros(no_limbs,1));
    for i = 1:no_limbs
        % replace diff(Q(t),t) with dQ(t)
        L = subs(L,diffQ(i),dQ(i));
        % derive wrt Q
        dLdqi = functionalDerivative(L,Q(i));
        % derive wrt dQ
        dLdqdoti = functionalDerivative(L,dQ(i));
        % replace dQ(t) with diff(Q(t),t)
        dLdqdoti = subs(dLdqdoti, dQ(i), diffQ(i));
        % apply lagrangian equation to get dynamic equation of system.
        tau(i) = diff(dLdqdoti)-dLdqi;
    end
    
    % replace all diff(x(t),t) variables with dx(t) equivalents for output
    % clarity.
    for i = 1:no_limbs
        tau = subs(tau, diffQ(i), dQ(i));
        tau = subs(tau, diff(diffQ(i)), ddQ(i));
        tau = subs(tau, diff(dQ(i)), ddQ(i));
        K = subs(K, diffQ(i), dQ(i));
        V = subs(V, diffQ(i), dQ(i));
    end
end