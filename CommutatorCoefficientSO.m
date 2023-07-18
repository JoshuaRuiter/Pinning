function coeff = CommutatorCoefficientSO(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % Given two roots alpha and beta, output the coefficient
    % N_{ij}^{alpha beta} (u,v)
    % that arises in the commutator of X(alpha,u) and X(beta,v)
    
    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(Root_System.IsRoot(beta))
    
    function d = compare(alpha,beta)
        d = find(alpha)<find(beta);
    end
    
    if ~Root_System.IsRoot(alpha+beta)
        % If alpha+beta is not a root, then the coefficient is zero
        % which is the same as saying that the commutator of X(alpha,u)
        % and X(beta,v) is X(alpha+beta,0)=I
        coeff = 0;
    else
        q = Root_System.Rank;
        diff = MatrixSize - 2*q;
        c = sym('c',[1 diff]);
        if (i==1 && j==1 && length(u) == diff && length(v) == diff && sum(alpha)== 1 && sum(beta)==1)
    
        % Placeholder
        coeff = -compare(alpha,beta)*(sum(c*u*v));
        
        elseif (i==1 && j==1 && sum(alpha) == 1 && sum(beta) == 0 && length(u) == diff && length(v) == 1)
            coeff = sum(-u*v);
        elseif (i==1 && j==1 && sum(alpha) == 0 && sum(beta) == 1 && length(u) == 1 && length(v) == diff)
            coeff = sum(u*v);
        elseif (i==2 && j==1 && sum(alpha) == 1 && sum(beta) == 0 && length(u) == diff && length(v) == 1)
            coeff = sum(-v^2*c/2*u);
        elseif (i==1 && j==2 && sum(alpha) == 0 && sum(beta) == 1 && length(u) == 1 && length(v) == diff)
            coeff = sum(v^2*c/2*u);
    
        end
        % [Xai, Xaj] = Xai+aj(N11) -- N11(u,v) =  e( c1u1v1 + ...
        % e is a function
    end
end