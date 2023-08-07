n = 7;
q = 2;

% Setting up the root system
Phi = RootSystem('B',q,n);

% Building the form matrix B
if n > 2*q
    vec_C = sym('c',[1,n-2*q]);
else
    vec_C = [];
end
Form = NIForm(n,q,1,vec_C,0,'symmetric bilinear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up roots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = zeros(1,n);
alpha(1) = 1;

beta = zeros(1,n);
beta(1) = 1;
beta(2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating root space elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms s
syms t
sym('u',[n-2*q,1]);
sym('v',[n-2*q,1]);

X_alpha_u = X_SO(n,Phi,Form,alpha,u);
X_beta_s = X_SO(n,Phi,Form,beta,s);
X_alpha_minus_beta_v = X_SO(n,Phi,Form,alpha-beta,v);
X_beta_minus_alpha_v = X_SO(n,Phi,Form,beta-alpha,v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating a Weyl group element
% and conjugation by it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r
W_r = W_SO(n,Phi,Form,beta-2*alpha,r)
conjugation1 = W_r*X_alpha_u*W_r^(-1)

iota_1_t = sym(zeros(n-2*q,1));
iota_1_t(1) = t;
X_alpha_iota_1_t = X_SO(n,Phi,Form,alpha,iota_1_t);

conjugation2 = W_r*X_alpha_iota_1_t*W_r^(-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating some commutators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
commutator1 = simplify(Commutator(X_alpha_u,X_beta_minus_alpha_v));

% syms N12;
% X_SO(n,Phi,Form,2*alpha-beta,N12)
% commutator1 = simplify(Commutator(X_beta_t, X_gamma_u))
% 
% sum11 = beta+gamma;
% assert(isequal(sum11,alpha))
% dim_V_sum11 = RootSpaceDimensionSO(n,Phi,sum11);
% N11 = sym('N11_',[dim_V_sum11,1]);
% X_sum11_N11 = X_SO(n,Phi,Form,sum11,N11);
% 
% sum12 = beta+2*gamma;
% dim_V_sum12 = RootSpaceDimensionSO(n,Phi,sum12);
% N12 = sym('N12_',[dim_V_sum12,1]);
% X_sum12_N12 = X_SO(n,Phi,Form,sum12,N12);
% 
% vec_1 = ones(1,dim_V_alpha);
% w_alpha_1 = simplify(W_SO(n,Phi,Form,alpha,vec_1));
% v = sym('v',[dim_V_alpha,1]);
% w_alpha_v = simplify(W_SO(n,Phi,Form,alpha,v));
% h_alpha_v = simplify(w_alpha_v*w_alpha_1^(-1));
% 
% delta = beta-alpha;
% dim_V_delta = RootSpaceDimensionSO(n,Phi,delta);
% assert(dim_V_delta==n-2*q)
% w = sym('w',[dim_V_delta,1]);
% X_delta_w = X_SO(n,Phi,Form,delta,w);
% 
% X_beta_t;
% commutator2 = simplify(Commutator(X_alpha_s, X_delta_w));
