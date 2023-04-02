function mat = X_SL(MatrixSize, root_system, FormMatrix, alpha, u)
    % Given a positive integer n, a root alpha in the root system of SL_n,
    % and symbolic variable u, output the associated element X_alpha(u) 
    % of the root subgroup

%     % Validating inputs
%     RootSystem = GetRootsSL(n);
%     assert(IsRoot(alpha,RootSystem));

    % In general, this should be the exponential of LieX(alpha,u)
    % It just so happens that the matrix exponential of LieX(alpha,u)
    % is just I + LieX(alpha,u) for the SLn case because 
    % LieX(alpha,u)^2 = 0
    mat = eye(MatrixSize) + LieX_SL(MatrixSize,root_system,FormMatrix,alpha,u);
end