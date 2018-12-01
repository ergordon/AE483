function simulation_template
% Parameters (method of collision avoidance)
params.r = 0.31;
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
obst = AddObstacle_Sphere(obst, [0.00; 0.95; -2.00], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; 0.00; -1.05], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; -0.95; -2.00], 0.51);
obst = AddObstacle_Sphere(obst, [0.00; -0.00; -2.95], 0.51);
    
    %%%%%%%%%%%%%%%%%%%%
    %
    % TO ADD A SPHERICAL OBSTACLE:
    %
    %   p       center of obstacle (3x1)
    %   s       radius of obstacle (1x1)
    %
    %   obst = AddObstacle_Sphere(obst, p, s);
    %
    % TO ACCESS THE CENTER AND RADIUS OF OBSTACLE "i":
    %
    %   obst{i}.p
    %   obst{i}.s
    %
    % TO GET THE NUMBER OF OBSTACLES:
    %
    %   length(obst)
    %
    %%%%%%%%%%%%%%%%%%%%
    
    % Planner (initialize o_desired)
    o_desired = o_start;
    
    % Data to store
    data.t = [params.t0];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    params.b_rep=.9%.2294551
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
        (norm(o_desired - o_goal))
    % Comment the next line to suppress the visualization (it's not
    % necessary to solve this problem)
    % visualize(data, params);

end

function o_desired = planner(o_desired, o_goal, obst, params)
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

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CODE TO VISUALIZE RESULTS (CAN BE IGNORED)

function visualize(data, params)
    
    fig = [];
    for i=1:length(data.t)
        fig = UpdateFigure(fig, data.t(i), data.o_desired(:, i), ...
                           data.o_goal, params.r, data.obst{:, i});
    end
    
end

function fig = UpdateFigure(fig, t, o_desired, o_goal, r, obst)

    if isempty(fig)

        % Create figure
        clf;
        set(gcf,'renderer','opengl');
        set(gcf,'color','w');
        title(sprintf('t = %6.2f', t));

        % Create axis
        axis equal;
        axis(8*[-1 1 -1 1 -1 0]);
        axis manual;
        hold on;

        % Reverse the y and z axes to get the "z down" view
        set(gca, 'ydir', 'reverse');
        set(gca, 'zdir', 'reverse');

        % Draw lights and labels
        lighting flat
        light('Position',[0 -2 -1])
        light('Position',[0 -2 1])
        xlabel('x');
        ylabel('y');
        zlabel('z');

        % Create unit sphere
        [fig.xSphere,fig.ySphere,fig.zSphere]=sphere(16);
        [m,n]=size(fig.xSphere);

        % Create template array for sphere color
        c = ones(m,n,3);

        % Draw obstacles
        for i=1:length(obst)
            if (obst{i}.type == 1)
                c(:,:,1) = 0.75;
                c(:,:,2) = 0.25;
                c(:,:,3) = 0.25;
                fig.obst(i) = ...
                    DrawBubble([], obst{i}. p,obst{i}.s, ...
                                   fig.xSphere, fig.ySphere, fig.zSphere, ...
                                   c, 0.8);
            end
        end

        % Draw desired position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.25;
        c(:,:,3) = 0.75;
        fig.qbub = ...
            DrawGlowingBubble([], o_desired, r, ...
                              fig.xSphere, fig.ySphere, fig.zSphere, ...
                              c, 0.4);

        % Draw goal position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.75;
        c(:,:,3) = 0.25;
        fig.gbub = ...
            DrawBubble([], o_goal, r, ...
                       fig.xSphere, fig.ySphere, fig.zSphere, ...
                       c, 0.2);

        % Make everything look pretty
        h = findobj('Type','surface');
        set(h,'FaceLighting','gouraud',...
              'FaceColor','interp',...
              'EdgeColor',[.4 .4 .4],...
              'LineStyle','none',...
              'BackFaceLighting','lit',...
              'AmbientStrength',0.4,...
              'DiffuseStrength',0.6,...
              'SpecularStrength',0.5);
        material default

    else

        title(sprintf('t = %6.2f', t));
        for i=1:length(obst)
            if (obst{i}.type == 1)
                fig.obst(i) = ...
                    DrawBubble(fig.obst(i), obst{i}.p, obst{i}.s, ...
                               fig.xSphere, fig.ySphere, fig.zSphere);
            end
        end
        fig.qbub = DrawGlowingBubble(fig.qbub, o_desired, r, ...
                                     fig.xSphere, fig.ySphere, fig.zSphere);
        fig.gbub = DrawBubble(fig.gbub, o_goal, r, ...
                              fig.xSphere, fig.ySphere, fig.zSphere);

    end

    drawnow;

end

function bubble = DrawGlowingBubble(bubble,o,r,x,y,z,c,a)
    if (isempty(bubble))
        bubble.surface = surf(o(1)+r*x,o(2)+r*y,o(3)+r*z,c,'FaceAlpha',a);
        bubble.light = light('position',o,'style','local');
    else
        set(bubble.surface,'xdata',o(1)+r*x,'ydata',o(2)+r*y,'zdata',o(3)+r*z);
        set(bubble.light,'position',o);
    end
end

function bubble = DrawBubble(bubble,o,r,x,y,z,c,a)
    if (isempty(bubble))
        bubble = surf(o(1)+r*x,o(2)+r*y,o(3)+r*z,c,'FaceAlpha',a);
    else
        set(bubble,'xdata',o(1)+r*x,'ydata',o(2)+r*y,'zdata',o(3)+r*z);
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

