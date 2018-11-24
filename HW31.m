function HW31
clc;clear;
p4
end

function p1
    % Find a symmetric matrix N that satisfies
    %     x'M x = x'Nx
    % for all x

    M = [1 7 1 1; -1 -3 4 10; -2 0 0 -2; 5 10 -10 -7];
    mat2str((M+M')/2);
end

function p2
    %Mark all of the statements that are true
    M1 = [6.0 7.0; 7.0 7.0];
    M2 = [4.0 6.0; 6.0 8.0];
    M3 = [1.0 6.5; 6.5 2.0];
    M4 = [1.0 5.0; 5.0 2.0];
    M5 = [8.0 4.0; 4.0 6.0];

    M=M1;

    eig_M = eig(M);
    flag = 0;
    for i = 1:rank(M)
        if eig_M(i) <= 0 
        flag = 1;
        end
    end
    if flag == 1
        disp('the matrix is not positive definite')
        else
        disp('the matrix is positive definite')
    end

end

function p3 
    F1 = [7.0 2.0 -0.5 3.5; 2.0 6.0 5.5 4.5; -0.5 5.5 8.0 2.5; 3.5 4.5 2.5 8.0];
    g1 = [3.0; -0.5; -4.5; 5.0];
    h1 = [-1.5];
    F2 = [8.0 0.5 3.0 5.0; 0.5 5.0 2.0 4.0; 3.0 2.0 1.0 4.5; 5.0 4.0 4.5 6.0];
    g2 = [-3.5; -3.5; 4.0; -4.5];
    h2 = [0.5];

    x = sym('x', [4 1],'real');

    %Minimizer Function
    minimizer1 = vpa(0.5*x'*F1*x + g1'*x + h1);
    minimizer2 = vpa(0.5*x'*F2*x + g2'*x + h2);

    %Minimizer
    sol1 = solve(jacobian(minimizer1,x)==0,x);
    sol2 = solve(jacobian(minimizer2,x)==0,x);
    %x1 = [sol1.x1; sol1.x2];
    %x2 = [sol2.x1; sol2.x2];
    %x1 = [sol1.x1; sol1.x2 ; sol1.x3];
    %x2 = [sol2.x1; sol2.x2; sol2.x3];
    x1 = [sol1.x1; sol1.x2 ; sol1.x3; sol1.x4];
    x2 = [sol2.x1; sol2.x2; sol2.x3; sol2.x4];


    %Minimum Values
    MV1 = vpa(subs(minimizer1,x,x1));
    MV2 = vpa(subs(minimizer2,x,x2));

    disp('Minimizer 1')
        disp(mat2str(double(x2)))
        disp(MV2)
    disp('Minimizer 2')
        disp(mat2str(double(x1)))
        disp(MV1)
end

function p4
    Ad = [-4.4 2.6 2.6; 4.8 4.0 1.5; 0.7 0.3 -3.1];
    Bd = [-4.3 -3.4; 2.1 4.8; 3.5 -0.6];
    Cd = [4.5 1.0 3.2];

    syms P K real
    xd0 = sym('xd0', [3 1], 'real');
    xd1 = sym('xd1', [3 1], 'real');
    yd1 = sym('yd1', [3 1], 'real');
    ud0 = sym('ud0', [2 1], 'real');

    %Minimize ud(0), xd(1)
    minFun = @(xd1,yd1,ud0) 0.5*(Cd*xd1-yd1)'*(Cd*xd1-yd1)+0.5*ud0'*ud0;
    %Such That
    x_d1 = @(xd0,ud0) Ad*xd0+Bd*ud0;

    minFun1 = @(xd1,yd1) 0.5*(Cd*xd1-yd1)'*(Cd*xd1-yd1);
    minFun2 = @(ud0) 0.5*ud0'*ud0;

    jacobian(minFun1(xd1,yd1),xd1) 
    jacobian(minFun2(ud0),ud0)

    %with a minimizer
    ud0 = @(xd0,yd1) -K*(Cd*Ad*xd0 - yd1);
    %minimum value has the form
    minval = @(xd0,yd1) 0.5*(Cd*Ad*xd0-yd1)'*P*(Cd*Ad*xd0-yd1);

    sol = solve(vpa(jacobian(0.5*((Cd*(Ad*xd0+Bd*ud0)-yd1))'*(Cd*(Ad*xd0+Bd*ud0)-yd1)+0.5*ud0'*ud0,xd0))==0,ud0)
    XX = [subs(sol.ud01,[xd0;yd1],[1;0;0;0;0;0]);...
          subs(sol.ud01,[xd0;yd1],[0;1;0;0;0;0]);...
          subs(sol.ud01,[xd0;yd1],[0;0;1;0;0;0])]
    YY = subs(sol.ud01,[xd0;yd1],[0;0;0;1;0;0])


    K = lqr(Ad,Bd,eye(3),eye(2))
end