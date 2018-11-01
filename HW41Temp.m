function simulation_template

    k_att = 1.35;
k_rep = 0.90;
k_des = 0.03;
q_goal = [4.13; -0.17; -3.41];
q = [-4.61; 0.77; -2.37];
r = 1.50;
p = [0.02; 0.52; -3.49];
s = 0.79;
t0 = 0.00;
t1 = 10.00;
dt = 0.05;

    % Create variables to keep track of time and desired position
    data.t = [t0];
    data.q = [q];
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(t1/dt)
        
        % Update desired position
        % 
        % q = ...
        dq = norm(p-q) - (s+r)
        fgrad = k_att * (q-q_goal) + k_rep*(1/dq^3)*((p-q)/norm(p-q))
        q = q-k_des*fgrad

        
        % Store time and desired position
        data.t(:, end+1) = data.t(:, end) + dt;
        data.q(:, end+1) = q;

    end
    
    mat2str(q)
    
    % Put other stuff to visualize into data
    data.q_goal = q_goal;
    data.r = r;
    data.p = p;
    data.s = s;
    
    % Comment the next line to suppress the visualization (it's not
    % necessary to solve this problem)
    visualize(data);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CODE TO VISUALIZE RESULTS (CAN BE IGNORED)

function visualize(data)
    
    fig = [];
    for i=1:length(data.t)
        fig = UpdateFigure(fig, data.t(i), data.q(:, i), data.q_goal, data.r, data.p, data.s);
    end
    
end

function fig = UpdateFigure(fig, t, q, q_goal, r, p, s)

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

        % Draw obstacle
        c(:,:,1) = 0.75;
        c(:,:,2) = 0.25;
        c(:,:,3) = 0.25;
        fig.pbub = DrawBubble([], p, s, fig.xSphere, fig.ySphere, fig.zSphere, c, 0.8);

        % Draw desired position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.25;
        c(:,:,3) = 0.75;
        fig.qbub = DrawGlowingBubble([], q, r, fig.xSphere, fig.ySphere, fig.zSphere, c, 0.4);

        % Draw goal position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.75;
        c(:,:,3) = 0.25;
        fig.gbub = DrawBubble([], q_goal, r, fig.xSphere, fig.ySphere, fig.zSphere, c, 0.2);

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
        fig.pbub = DrawBubble(fig.pbub, p, s, fig.xSphere, fig.ySphere, fig.zSphere);
        fig.qbub = DrawGlowingBubble(fig.qbub, q, r, fig.xSphere, fig.ySphere, fig.zSphere);
        fig.gbub = DrawBubble(fig.gbub, q_goal, r, fig.xSphere, fig.ySphere, fig.zSphere);

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

