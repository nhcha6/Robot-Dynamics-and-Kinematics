function [ p, t, cfs, A, b ] = cubic_spline(pos, t_int, plt, t_res )
% CUBIC_SPLINE Trajectory generation with cubic splines. Returns coefficients of cubics.
% Written by Wesley Au
% Inputs
%   pos     [1xn] vector of n-angles representing the desired trajectory
%   t_int   [1xn] vector of cumulative time for n-points starting from 0
%   plt     (Optional) [1x1] scalar; set to 1 to plot trajectory profile
%   t_res   (Optional) [1x1] scalar; time resolution for output trajectory
% Outputs
%   p       [mx1] vector representing cubic spline trajectory
%   t       [mx1] vector of time for each point in p
%   cfs     Rows of cubic coefficients, ordered from 0th to 3rd degree.
%           Each row is a different spline.
%   A       A matrix.
%   b       b vector.

% Initialise variables
n_cubics = length(pos)-1;
A = zeros(4*n_cubics);
b = zeros(4*n_cubics, 1);

% If an incorrect time vector is given, assume 1 second intervals
if t_int(1) ~= 0
    if length(t_int) < length(pos)
        t_int = [0 t_int];
    else
        t_int = t_int - t_int(1);
    end
end

if nargin < 3 || isempty(plt)
    plt = false;
end

if nargin < 4
    t_res = 0.25;
end

tf = diff(t_int);       % Normalise time

%% 2) ========= A AND b MATRIX POPULATION ======================================
% Assume each cubic is time-normalised, i.e. starts at t = 0.
% Therefore, generate coefficients for each cubic such that t = 0
% represents the initial position of each via point.

% declare first two rows as these aren't repeated in the iteration.
A(1,1) = 1;
b(1) = pos(1);
A(2,2) = 1;
% declare last two rows as these aren't repeated in the iteration
Tf = t_int(end)-t_int(end-1);
A(end-1,end-3:end) = [1, Tf, Tf^2, Tf^3];
b(end-1) = pos(end);
A(end,end-2:end) = [1, 2*Tf, 3*Tf^2];
% A matrix and b vector population code here
% loop through and add the 4 constraints added by an additional via point.
for i = 2:n_cubics
    % j is the row of A and b we are up to
    j = 4*(i-1)-1;
    % calculate time-step T
    T = t_int(i)-t_int(i-1);
    % fill in A and b
    A(j,j-2:j+1) = [1, T, T^2, T^3];
    b(j) = pos(i);
    A(j+1,j-1:j+3) = [1, 2*T, 3*T^2, 0, -1];
    b(j+1) = 0;
    A(j+2,j:j+4) = [2, 6*T, 0, 0, -2];
    b(j+2) = 0;
    A(j+3,j+2) = 1;
    b(j+3)=pos(i);
end
% ==========================================

%% 2A) ======== SOLVE COEFFICIENTS =============================================
m = rref([A b]);                        % same as solving via inverse(A) method
cfs = reshape(m(:,end),4,[])';

%% 3) ========= CREATE TRAJECTORY ==============================================
% Create position and velocity trajectories
% Initialise variables
p = cell(n_cubics,1);
t = p;
for i = 1:n_cubics
    c = cfs(i,end:-1:1);                % Obtain cubic vector
    n = fix(tf(i)/t_res)+1;             % Number of points after initial time to tf
    t1 = linspace(0,tf(i),n);           % Generate time vector between 0 to tf
    p1 = polyval(c, t1);                % Evaluate polynomial for time vector t1
    % Append to vectors p, v and t
    if i == 1
        % If first spline, use all points
        p{i} = p1;
        t{i} = t1 + t_int(i);
    else
        % Otherwise, truncate first value
        p{i} = p1(2:end);
        t{i} = t1(2:end) + t_int(i);
    end
    
end
p = [p{:}]';
t = [t{:}]';


%% 4) ========= PLOT CODE ======================================================
% Runs only if plt == 1
if plt
    f = figure;
    ax = axes(f);
    n = 21;                                    % Number of points to plot
    set(ax, 'ColorOrder', lines(n_cubics));    % Plot subsequent lines at different colours
    
    ylbl = {'Position units', 'Velocity units', 'Acceleration units'};
    
    for i = 1:3
        % Plot position profile
        subplot(3,1,i)
        hold all
        for j = 1:n_cubics
            t1 = linspace(0,tf(j),n);
            c = cfs(j,end:-1:1);
            for k = 2:i
                c = polyder(c);
            end
            v = polyval(c, t1);                 % Evaluate polynomial with t1
            x = t1 + t_int(j);                      % Time adjust x-axis
            plot(x,v,'LineWidth',2)

        end
        if i == 1
            title('Trajectory profile')
            plot(t_int,pos,'o','MarkerEdgeColor','k')       % Mark via points
        end
        ylabel(ylbl{i})
        grid on

    end

    xlabel('Time units')
    
end

end
