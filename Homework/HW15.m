function hw15
clc; clear;
P5
end

function x = P2()
    o_1in0 = [0.85185365; 0.44955071; -0.86297332];
    R_1in0 = [0.12233314 -0.42942110 0.89478049; -0.12359739 0.88794231 0.44303740; -0.98476306 -0.16479069 0.05554945];
    v_01in0 = [-0.44550600; -0.94736152; -0.33312745];
    w_01in1 = [-0.56710857; 0.73383419; -0.74890073];
    v_in0 = [-0.57025088; -0.63420489; -0.77029285];
    v_in0_dot = [0.33338284; -1.21014159; 0.93525262];


    mat2str(R_1in0*(skew(w_01in1)*v_in0+v_in0_dot))
end

function X = P4()
    R_1in0 = [-0.05057978 0.64653627 -0.76120466; -0.92407314 0.25882904 0.28124076; 0.37885422 0.71763387 0.58435529];
    w_01in1 = [-0.27382255; -0.16958297; -0.30432453];
    o_1in0 = [0.16794702; -0.73973547; -0.69218719];
    v_01in0 = [-0.00605263; -0.87903643; -0.31473725];

    mat2str(o_1in0-R_1in0*skew(w_01in1)*v_01in0)
end

function X = P5()
    R_1in0 = [0.95797351 0.04678506 -0.28301575; 0.26547766 -0.51834199 0.81292262; -0.10866631 -0.85389270 -0.50897829];
    w_10in0 = [-0.41227882; 0.40207768; -0.78628338];

    X = mat2str(-R_1in0'*w_10in0)
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