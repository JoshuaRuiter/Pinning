%{
c1 = sym("c1");
c2 = sym("c2");
B = GetB(6,2,[c1 c2]);
r = RootSystem("B",2,6);
u = sym("u");
v = sym(["v1" "v2"]);
alpha = [1 -1 0 0 0 0];
beta = [1 0 0 0 0 0];
disp("In SO_62, the Pinning X12 and X1 are")
X_alpha_12 = X_SO(6,r,B,alpha,u)
X_alpha_1 = X_SO(6,r,B,beta,v)
%}
c3 = sym('c',[1 3]);
B = GetB(7,2,c3);
r = RootSystem("B",2,7);
u = sym("u");
v = sym(["v1" "v2" "v3"]);
alpha = [1 -1 0 0 0 0 0];
beta = [1 0 0 0 0 0 0];
%disp("In SO_62, the Pinning X12 and X1 are")
X_alpha_12 = X_SO(7,r,B,alpha,u)
X_alpha_1 = X_SO(7,r,B,beta,v)