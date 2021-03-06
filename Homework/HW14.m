function hw14
clc; clear; 
P1
end

function [th1,th2,th3] = P1()
    %ZYX
    R_1in0 = [0.821156019723 0.454239466298 0.345498623051; 0.041697129882 -0.651524878455 0.757480482992; 0.569178478715 -0.607603357472 -0.553944057966];
    syms c1 s1 c2 s2 c3 s3
    R =  z(c1,s1)*y(c2,s2)*x(c3,s3)
    th2 = -asin(R_1in0(3,1))
end

function w011 = P4()
    theta1 = 0.90;
    theta2 = 0.69;
    theta3 = -1.07;
    theta1dot = 1.00;
    theta2dot = 0.79;
    theta3dot = -0.23;

    syms c1 s1 c2 s2 c3 s3

    R =  y(c2,s2)*x(c3,s3)

    w011 = R*[theta1dot; 0; 0];
end

%Make Skew Matrix
function S = skew(w)
    S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

%% Rotation Matrix Maker
function X = x(c,s)
    X = [1 0 0; 0 c -s; 0 s c];
end

function Y = y(c,s)
    Y= [c 0 s; 0 1 0; -s 0 c];
end

function Z = z(c,s)
    Z = [c -s 0; s c 0; 0 0 1];
end