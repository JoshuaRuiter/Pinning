function val = WeylGroupConjugationCoefficientSL(MatrixSize,root_system,FormMatrix,alpha,beta,u)
    % The function computes the sign (+1 or -1) that appears as a result of 
    % conjugating X_beta(v) by a Weyl group element

    % That is, 
    % w_alpha(u)*X_beta(v)*w_alpha(u)^(-1) = X_{sig_alpha(beta)}(cv)
    % and this function computes c, which is always 1 or -1

    % Note that this sign is NOT very important for theoretical purposes.
    % In terms of proving theorems, it is only useful to know that c is
    % always 1 or -1, not usually useful to know exactly when it is 1 or -1
    % But anyway, we can calculate it here.

    % We can always write alpha=alpha_ij=alpha_i-alpha_j
    % and beta=alpha_kl=alpha_k-alpha_l for some i,j,k,l
    % This finds the index of the first entry of alpha which is 1
    i = find(alpha==1);
    %j = find(alpha==-1);
    k = find(beta==1);
    l = find(beta==-1);

    % I don't understand why it is -1 when i=k or i=l
    % and 1 otherwise, but it does pass the tests.
    val = u;
    if i==k || i==l
        val = -u;
    end
end
