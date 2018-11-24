function Exam4
clc; clear;
P5
end

%% Problem 1
function P1
%%%%% Given %%%%%%%
k_att = 0.48;
b_att = 3.96;
q_goal = [-0.65; 1.99; 1.07];
k_rep = 0.49;
b_rep = 1.59;
q = [-1.74; -1.23; 0.83];
d = 1.74;
dgrad = [4.35; 0.31; 4.83];

%%%%% Computing %%%%%%%
if ( norm(q-q_goal) <= b_att )
    gradf = k_att*(q-q_goal);
else
    gradf = k_att*b_att*(q-q_goal)/norm(q-q_goal);
end

if d <= b_rep
    fgrad = gradf + (-k_rep*((1/d)-(1/b_rep))*(1/d)^2)*dgrad;
else
    fgrad = gradf + zeros(size(q_goal));
end

%%%%% Solution %%%%%%%
mat2str(fgrad)
end

%% Problem 2
function P2
%%%%% Given %%%%%%%
q = [-1.85; 1.04; -0.20];
k_des = 0.56;
b_des = 1.85;
fgrad = [4.30; 0.43; -0.47];

%%%%% Computing %%%%%%%
if norm(k_des*fgrad) <= b_des
    q = q - k_des*fgrad;
else
    q = q - b_des*(fgrad/norm(fgrad));
end

%%%%% Solution %%%%%%%
mat2str(q)
end

%% Problem 3
function P3
%%%%% Given %%%%%%%
k_att = 0.63;
b_att = 4.04;
q_goal = [-1.93; -0.22; 1.79];
q = [-1.09; -1.26; -0.80];

%%%%% Computing %%%%%%%
if ( norm(q-q_goal) <= b_att )
    gradf = k_att*(q-q_goal);
else
    gradf = k_att*b_att*(q-q_goal)/norm(q-q_goal);
end

%%%%% Solution %%%%%%%
mat2str(gradf)
end

%% Problem 4
function P4
% Parameters (method of collision avoidance)
params.r = 1.00;
params.k_att = 0.76;
params.b_att = 0.96;
params.k_rep = 0.98;
params.b_rep = 0.92;
params.k_des = 0.03;
params.b_des = 0.23;

% Parameters (simulation)
params.t0 = 0.00;
params.t1 = 10.00;
params.dt = 0.05;

% Start and goal position
o_start = [-4.41; 0.79; -3.48];
o_goal = [4.57; 0.61; -3.62];

% Obstacles
obst = {};
obst = AddObstacle_Sphere(obst, [2.30; -3.67; -0.18], 1.30);
obst = AddObstacle_Sphere(obst, [-3.65; -3.71; -1.66], 1.14);
obst = AddObstacle_Sphere(obst, [-2.66; -2.64; -3.05], 1.13);
obst = AddObstacle_Sphere(obst, [-2.95; -0.24; -0.91], 1.41);
obst = AddObstacle_Sphere(obst, [2.69; -2.71; -1.27], 1.18);
obst = AddObstacle_Sphere(obst, [-2.90; -1.49; -2.34], 0.91);
obst = AddObstacle_Sphere(obst, [-0.06; -3.43; -2.36], 1.25);
obst = AddObstacle_Sphere(obst, [-1.94; 3.74; -3.64], 0.27);
obst = AddObstacle_Sphere(obst, [0.04; 2.93; -0.03], 0.86);
obst = AddObstacle_Sphere(obst, [3.92; 1.85; -0.95], 0.66);
    
    % Planner (initialize o_desired)
    o_desired = o_start;
    
    % Data to store
    data.t = [params.t0];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(params.t1/params.dt)
        
        % Planner (update o_desired)
        o_desired = planner(o_desired, o_goal, obst, params);
        
        % Store data
        data.t(:, end+1) = data.t(:, end) + params.dt;
        data.o_desired(:, end+1) = o_desired;
        data.o_goal(:, end+1) = o_goal;
        data.obst{:, end+1} = obst;

    end
end

function o_desired = planner(o_desired, o_goal, obst, params)
% Rename o_desired, o_goal, and the radius of the bounding volume for convenience
 q = o_desired;
 q_goal = o_goal;
 r = params.r;
 
  % Get attractive part of gradient
 if ( norm(q-q_goal) <= params.b_att )
    gradf = params.k_att*(q-q_goal);
 else
    gradf = params.k_att*params.b_att*(q-q_goal)/norm(q-q_goal);
 end
 
 % Get repulsive part of gradient
 for i=1:length(obst)
     % Get center and radius of i'th spherical obstacle
     p = obst{i}.p;
     s = obst{i}.s;
     % Compute distance and gradient of distance to this obstacle
     d = norm(q-p)-(s+r);
     dgrad = (q-p)/norm(q-p);
     %Compute repulsive gradient for this obstacle and add it to total gradient
    if ( d <= params.b_rep )
        gradf = gradf + params.k_rep*(1/d - 1/params.b_rep)*(-1)/d^2*dgrad;
    end
 end
 
 % Take a step
 if ( params.k_des*gradf <= params.b_des )
    q = q - params.k_des*gradf;
 else
    q = q - params.b_des*(gradf/norm(gradf));
 end
 
 % Recover o_desired

    o_desired = q;
    mat2str(q)
end

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end

%% Problem 5
function P5
% Parameters (method of collision avoidance)
params.r = 0.33;
params.k_att = 2.00;
params.b_att = 1.00;
params.k_rep = 2.00;
params.b_rep = 1.00;
params.k_des = 0.10;
params.b_des = 0.10;

% Parameters (simulation)
params.t0 = 0.00;
params.t1 = 20.00;
params.dt = 0.05;

% Start and goal position
o_start = [-5.48; -0.10; -2.10];
o_goal = [3.83; -0.10; -1.90];

% Obstacles
obst = {};
obst = AddObstacle_Sphere(obst, [0.00; 0.95; -2.00], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; 0.29; -1.10], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; -0.77; -1.44], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; -0.77; -2.56], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; 0.29; -2.90], 0.51);

    % Planner (initialize o_desired)
    o_desired = o_start;
    
    % Data to store
    data.t = [params.t0];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    % Define b_rep_bounds
    b_rep_bounds = [2,10];
    res = 1;
    while (res == 1)
        
    if abs(b_rep_bounds(2) - b_rep_bounds(1)) > 1e-10
        params.b_rep = 0.5*(b_rep_bounds(1)+b_rep_bounds(2));
        res = 1;
    else
        res = 0;
    end
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(params.t1/params.dt)
        % Planner (update o_desired)
        o_desired = planner2(o_desired, o_goal, obst, params);

        % Store data
        data.t(:, end+1) = data.t(:, end) + params.dt;
        data.o_desired(:, end+1) = o_desired;
        data.o_goal(:, end+1) = o_goal;
        data.obst{:, end+1} = obst;
    end
    
    error = norm(o_desired - o_goal);
    if error < 1e-10
       b_rep_bounds(1) = params.b_rep;
    else
       b_rep_bounds(2) = params.b_rep;
    end
    
    o_desired = o_start;
    end
    b_rep_bounds
end

function o_desired = planner2(o_desired, o_goal, obst, params)
% Rename o_desired, o_goal, and the radius of the bounding volume for convenience
 q = o_desired;
 q_goal = o_goal;
 r = params.r;
        

 
  % Get attractive part of gradient
 if ( norm(q-q_goal) <= params.b_att )
    gradf = params.k_att*(q-q_goal);
 else
    gradf = params.k_att*params.b_att*(q-q_goal)/norm(q-q_goal);
 end
 
 % Get repulsive part of gradient
 for i=1:length(obst)
     % Get center and radius of i'th spherical obstacle
     p = obst{i}.p;
     s = obst{i}.s;
     % Compute distance and gradient of distance to this obstacle
     d = norm(q-p)-(s+r);
     dgrad = (q-p)/norm(q-p);
     %Compute repulsive gradient for this obstacle and add it to total gradient
    if ( d <= params.b_rep )
        gradf = gradf + params.k_rep*(1/d - 1/params.b_rep)*(-1)/d^2*dgrad;
    end
 end
 
 % Take a step
 if ( params.k_des*gradf <= params.b_des )
    q = q - params.k_des*gradf;
 else
    q = q - params.b_des*(gradf/norm(gradf));
 end
 
 % Recover o_desired
 o_desired = q;   
end

function obst = AddObstacle_Sphere2(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end