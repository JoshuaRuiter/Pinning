%SL3
r= RootSystem("A",2,3);
alpha = [1 -1 0];
beta = [0 1 -1];
u= sym("u");
v= sym("v");
disp("The Pinning X12 and X32 are")
X_SL(3,0,0,alpha, u)
X_SL(3,0,0,beta,v)
disp("the commutator coefficient of alpha and beta is")
CommutatorCoefficientSL(3,r,alpha,beta,1,1,u,v)
Commutator(X_SL(3,0,0,alpha, u),X_SL(3,0,0,beta,v))

%SO62
c1 = sym("c1");
c2 = sym("c2");
B = [0 0 0 0 1 0;
     0 0 c1 0 0 1;
     0 0 0 c2 0 0;
     1 0 0 0 0 0;
     0 1 0 0 0 0];
r = RootSystem("B",2,2);
u = sym("u");
v = sym(["v1" "v2"]);
alpha = [1 -1];
beta = [1 0];
disp("The Pinning X12 and X2 are")
X_SO(6,r,B,alpha,u)
X_SO(6,r,B,beta,v)