function mat = LieGroup_SO(n,q)
%LIEX Summary of this function goes here
%   Detailed explanation goes here
assert(n >= 2*q+2);
diff = n-2*q;
I = eye(n);
S = sym(I);
s = sym('s',[1 q]);
for i=1:q
    S(i,i) = s(i);
    S(q+i,q+i) = s(i)^(-1);
end
c = sym('c',[1 diff]);
C = diag(c);
X = sym('X%d%d', [n n]);

X(1:q,2*q+1:n) = -transpose(X(2*q+1:n,q+1:2*q))*C;
X(q+1:2*q,2*q+1:n) = -transpose(X(2*q+1:n,1:q))*C;
X(q+1:2*q,q+1:2*q) = -transpose(X(1:q,1:q));

for i = 1:q
    X(i,q+i) = 0;
    X(i,q+i:2*q) = - X(i:q,q+i);
    X(q+i,i) = 0;
    X(q+i,i:q) = - X(q+i:2*q,i);
end

% for i=1:q
%     X(i,2*q+1) = -c(1)*X(2*q+1, q+i);
%     X(i,2*q+2) = -c(2)*X(2*q+2, q+i);
%     X(q+2+i,q+1) = -c(1)*X(q+1, i);
%     X(q+2+i,q+2) = -c(2)*X(q+2, i);
%     X(i,q+2+i) = sym(0);
%     X(q+2+i,i) = sym(0);
%     X(q+2+1:n,q+2+1:n) = -transpose(X(1:q,1:q));
    %for j=1:q
        %X(q+2+i,q+2+j) = -X(i,j);
    %end

%     for k=i:q
%         X(k,q+2+i) = -X(i,q+2+k);
%         X(q+2+k, i) = -X(q+2+i, k);
% 
%     end
% 
% 
% end
% 
for i= 1:diff
    X(2*q+i,2*q+i) = sym(0);
end
% X(q+2,q+1) = -c(1)/c(2)*X(q+1,q+2);

X(2*q+1:n,2*q+1:n)= -C^(-1)*transpose(X(2*q+1:n,2*q+1:n))*C;

sXs = S*X*S^(-1);
mat = sXs;

end

