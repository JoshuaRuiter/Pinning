function mat = W_SL(MatrixSize,root_system,FormMatrix,alpha,u)
    % Return the Weyl group element w_alpha(u) associated with alpha
    % Note that this is undefined when u=0, so maybe there should be
    % something like assert(u~=0)
    mat = X_SL(MatrixSize,root_system,FormMatrix,alpha,u)*...
        X_SL(MatrixSize,root_system,FormMatrix,-alpha,-u^(-1))*...
        X_SL(MatrixSize,root_system,FormMatrix,alpha,u);
end