%{
A1 = [-0.24480328 0.82909520 -0.08621209; 0.01146697 0.08751085 1.10328562; 1.07060318 -0.04335231 -0.03380335];
A2 = [0.65538327 -0.73484583 -0.17456911; 0.68717188 0.48419935 0.54161407; -0.31347659 -0.47492378 0.82230154];
A3 = [-0.42078512 -0.42372897 -0.80211822; -0.38626443 -0.71636016 0.58105758; -0.82081647 0.55433012 0.13776227];
A4 = [0.09823347 0.07508717 0.95467921; 0.65507671 -0.68037222 -0.06725167; 0.62546561 0.66956714 0.01345332];
A5 = [-0.08072661 0.25060643 -0.96471738; -0.99647131 -0.04260740 0.07231555; -0.02298136 0.96715098 0.25316167];


disp('A1')
det(A1)
A1'*A1

disp('A2')
det(A2)
A2'*A2

disp('A3')
det(A3)
A3'*A3

disp('A4')
det(A4)
A4'*A4

disp('A5')
det(A5)
A5'*A5

%}

%{

a = [6; 10; -3];

ahat = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]

%}


v_in5 = [-5; 6; 6];
syms a b c
v_in5'*[a;b;c]
