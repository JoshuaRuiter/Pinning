% This script is intended to help answer questions about
% the "multipliable roots" in the BC_q root system,
% that is, the short roots. These are "multipliable," meaning that
% for such a root, a multiple of it is also a root.
% More concretely, for a short root alpha in BC_q, 2alpha is also a root.

q = 2;
n = 6;
BC_q = RootSystem('BC',q,n);

vec_C = sym('c',[1,q]);
s = GenericTorusElementSU(n,q,vec_C);
s_inverse = s^(-1);
H = GetH(n,q,vec_C);

for i=1:length(BC_q.RootList)
    alpha = BC_q.RootList{i};
    if dot(alpha,alpha)==1 
        % alpha is a short root, so it is a multipliable root
        assert(BC_q.IsRoot(2*alpha))
        
        dim_V_alpha = RootSpaceDimensionSU(BC_q,alpha);
        assert(dim_V_alpha==4);
        u = sym('u',[1,dim_V_alpha]);
        X_alpha_u = X_SU(n,BC_q,H,alpha,u);
        conjugation = s*X_alpha_u*s_inverse;

        alpha_of_s = PinnedGroup.CharacterEval(n,alpha,s);
        X_alpha_alpha_of_s_u = X_SU(n,BC_q,H,alpha,alpha_of_s*u);

        dim_V_2alpha = RootSpaceDimensionSU(BC_q,2*alpha);
        assert(dim_V_2alpha==1);
        v = sym('v',[1,dim_V_2alpha]);
        X_2alpha_v = X_SU(n,BC_q,H,2*alpha,v);

        alpha
        X_alpha_u
        conjugation
        X_alpha_alpha_of_s_u
        X_2alpha_v

        assert(SymbolicIsEqual(conjugation,X_alpha_alpha_of_s_u))
        

    end
end

function mat = GenericTorusElementSU(MatrixSize, RootSystemRank, DiagonalValues)
    % given a vector [t_1, t_2, t_3, ... t_q] of length q = RootSystemRank
    % build the symbolic diagonal matrix 
    % diag(t_1, ... t_q, t_1^(-1), ... t_q^(-1), 1, ... 1)
    % Note: t1, ... t_q, etc. should be "purely real"
    
    assert(length(DiagonalValues)==RootSystemRank);

    % Add an assumption that each entry of DiagonalValues is real
    for i=1:length(DiagonalValues)
        t_i = DiagonalValues(i);
        assumeAlso(t_i, 'real');
    end

    mat = SymbolicZeros(MatrixSize);
    for i=1:RootSystemRank
        mat(i,i) = DiagonalValues(i);
        mat(RootSystemRank+i,RootSystemRank+i) = DiagonalValues(i)^(-1);
    end
    for i=2*RootSystemRank+1:MatrixSize
        mat(i,i) = 1;
    end
end