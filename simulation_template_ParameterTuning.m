function simulation_template
    clc; clear;
    %% Changing Parameters
    
    % Parameters (method of collision avoidance)
    params.r = 0.26;
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
    o_start = [-5.00; -0.10; -2.10];
    o_goal = [5.00; -0.10; -1.90];

    % Obstacles
    obst = {};
    obst = AddObstacle_Sphere(obst, [0.00; 1.00; -2.00], 0.51);
    obst = AddObstacle_Sphere(obst, [0.00; 0.31; -1.05], 0.51);
    obst = AddObstacle_Sphere(obst, [0.00; -0.81; -1.41], 0.51);
    obst = AddObstacle_Sphere(obst, [0.00; -0.81; -2.59], 0.51);
    obst = AddObstacle_Sphere(obst, [0.00; 0.31; -2.95], 0.51);

    %% Initialization Code
    
    % Planner (initialize o_desired)
    o_desired = o_start;
    
    % Data to store
    data.t = [params.t0];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    %% What is the largest value of b_rep > 1 quadrotor no stuck (4.4.3)
    %{
    % Define b_rep_bounds
    %     If the final output is equal to the min b_rep_bound, 
    %     lower your b_rep_max
    %
    %     If the final output is equal to the max b_rep_bound, 
    %     increase your b_rep_min
   
    b_rep_min = 1;
    b_rep_max = 10;
    
    b_rep_bounds = [b_rep_min b_rep_max];

    while 1
        % Use Bisection Method to compute b_rep 
        if abs(b_rep_bounds(2) - b_rep_bounds(1)) > 1e-10
            params.b_rep = 0.5*(b_rep_bounds(1)+b_rep_bounds(2));
        else
            break
        end

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

        % Decide on the what b_rep to use next by changing b_rep_bounds
        error = norm(o_desired - o_goal);
        if error < 1e-10
           b_rep_bounds(1) = params.b_rep;
        else
           b_rep_bounds(2) = params.b_rep;
        end

        % Set Quadrotor to starting position
        o_desired = o_start;
    end
    
    disp(b_rep_bounds);
    %}
    
    %% What is the smallest value of b_rep > 1 quadrotor no stuck (4.4.2)
    %{
    % Define b_rep_bounds
    %     If the final output is equal to the max b_rep_bound, 
    %     lower your b_rep_max
    
    b_rep_min = 1;
    b_rep_max = 5;
    
    b_rep_bounds = [b_rep_min b_rep_max];
    
    while 1
        
        % Use Bisection Method to compute b_rep
        if abs(b_rep_bounds(2) - b_rep_bounds(1)) > 1e-10
            params.b_rep = 0.5*(b_rep_bounds(1)+b_rep_bounds(2));
        else
            break
        end

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

        % Decide on the what b_rep to use next by changing b_rep_bounds
        error = norm(o_desired - o_goal);
        if error < 1e-10
           b_rep_bounds(2) = params.b_rep;
        else
           b_rep_bounds(1) = params.b_rep;
        end

        % Set Quadrotor to starting position
        o_desired = o_start;
    end
    disp(b_rep_bounds);
    %}
    
    %% What is the largest value of b_rep < 1 quadrotor no stuck (4.4.1)
    %{
    % Define b_rep_bounds
    %     If the b_rep_min stays constant for multiple iterations
    %       lower b_rep_max to be equal or slightly above that value that
    %       b_rep_min got stuck at and rerun. 
    
    b_rep_min = 0;
    b_rep_max = 1;
    
    b_rep_bounds = [b_rep_min b_rep_max];
    
    while 1
    
        % Use Bisection Method to compute b_rep
        if abs(b_rep_bounds(2) - b_rep_bounds(1)) > 1e-10
            params.b_rep = 0.5*(b_rep_bounds(1)+b_rep_bounds(2));
        else
            break
        end

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
        
        % Decide on the what b_rep to use next by changing b_rep_bounds
        error = norm(o_desired - o_goal);
        if error < 1e-10
           b_rep_bounds(1) = params.b_rep;
        else
           b_rep_bounds(2) = params.b_rep;
        end
        
        disp(b_rep_bounds);
        
        % Set Quadrotor to starting position
        o_desired = o_start;
    end
    %}
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
        % Compute repulsive gradient for this obstacle and add it to total gradient
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

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end


