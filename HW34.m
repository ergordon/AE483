function HW34
clc; clear;
Ad = [0.5 0.6 -0.4; -0.4 0.3 -0.6; -0.7 0.8 0.5];
Bd = [0.0; -0.4; -0.7];
Qd = [0.8 0.0 0.0; 0.0 0.7 0.0; 0.0 0.0 0.4];
Rd = [0.5];

[K2, P2] = dlqr(Ad,Bd,Qd,Rd);
[K, P] = myDLQR(Ad,Bd,Qd,Rd);

mat2str(P2)
mat2str(K2)
mat2str(P)
mat2str(K)
end

function [K,P] = myDLQR(Ad,Bd,Qd,Rd)
tol = 0.1;
P = eye(size(Ad));
n = 1;
while(1)
    K = inv(Rd+Bd'*P*Bd)*Bd'*P*Ad;
    Pnew = Qd+Ad'*P*Ad-Ad'*P*Bd*K;
    if norm(P-Pnew) < tol
        P = Pnew;
        return
    else
        P=Pnew;
        n=n+1;
    end
end
end