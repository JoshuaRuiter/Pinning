function mat = LieX_SO_1(n,q)
%LIEX Summary of this function goes here
%   Detailed explanation goes here
assert(n == 2*q+2);
I = eye(n);
S = sym(I);
s = sym('s',[1 q]);
for i=1:q
    S(i,i) = s(i);
    S(n-q+i,n-q+i) = s(i)^(-1);
end
c = sym('c',[1 2]);
X = sym('X%d%d', [n n]);
for i=1:q
    X(i,q+1) = -c(1)*X(q+1, q+2+i);
    X(i,q+2) = -c(2)*X(q+2, q+2+i);
    X(q+2+i,q+1) = -c(1)*X(q+1, i);
    X(q+2+i,q+2) = -c(2)*X(q+2, i);
    for j=1:q
        X(i,q+2+j) = sym(0);
        X(q+2+i,j) = sym(0);
        X(q+2+i,q+2+j) = -X(i,j);
    end
end

X(q+1,q+1) = sym(0);
X(q+2,q+2) = sym(0);
X(q+2,q+1) = -c(1)/c(2)*X(q+1,q+2);

sXs = S*X*S^(-1);
mat = sXs;

end

