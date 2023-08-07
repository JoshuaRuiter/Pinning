function val = CommutatorCoefficientSL(MatrixSize,Root_System,Form,alpha,beta,p,q,u,v)
    % Given two roots alpha and beta, output the sign coefficient (1 or -1)
    % that arises in the commutator of X(alpha,u) and X(beta,v)

    % validating inputs
    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(Root_System.IsRoot(beta))
    assert(p==1); % for the type A root system, p and q can only ever be 1
    assert(q==1); % for the type A root system, p and q can only ever be 1
    
    % If alpha+beta is not a root, the coefficient is zero
    val = 0;

    if ~Root_System.IsRoot(alpha+beta)
        % If alpha+beta is not a root, then the coefficient is zero
        % which is the same as saying that the commutator of X(alpha,u)
        % and X(beta,v) is X(alpha+beta,0)=Identity
        val = 0;

    else
        % In this case, alpha+beta is a root, aka an element of Phi(SL_n,T)
        % In this case, the value is either 1 or -1
        assert(Root_System.IsRoot(alpha+beta))

        % We can always write alpha=alpha_ij
        % and beta=alpha_kl for some i,j,k,l
        % This finds the index of the first entry of alpha which is 1
        i = find(alpha==1); 
        j = find(alpha==-1);
        k = find(beta==1);
        l = find(beta==-1);

        % If alpha+beta is a root, then there must be exactly one
        % overlapping index between two indices corresponding
        % to oppositely signed entries of alpha and beta.
        % In other words, there are exactly two possibilities:
        % (1) j==k and i,j,l are distinct
        % (2) i==l and i,j,k are distinct
        assert(j==k || i==l)

        % The coefficient should be uv in case (1)
        % and -uv in case (2)
        if j==k
            val = u*v;
        elseif i==l
            val = -u*v;
        else
            % This should be impossible, throw an error if you get here
            assert(false, "A sum of two roots in A_(n-1) is a root, but" + ...
                "the roots are not of the right form.");
        end
    end
end