function mat = W_SU(MatrixSize,Root_System,FormMatrix,alpha,u)
    % Return the Weyl group element w_alpha(u) associated with alpha
    % Note that this is undefined when u is the zero vector

    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(length(u) == RootSpaceDimensionSU(Root_System,alpha))

    if length(u) == 1
        % If u is length 1, alpha should be a long root
        assert(dot(alpha,alpha)==4)
        u_prime = -u^(-1);

    elseif length(u) == 2
        % If u is a vector of length 2, which happens when alpha is a medium
        % length root, switch to viewing u as a complex number (element of the
        % field L)
        assert(dot(alpha,alpha)==2)
        if isa(u,'sym')
            assumeAlso(u(1),'real')
            assumeAlso(u(2),'real')
        end
        u_complex = u(1) + 1i*u(2);
        
        % u_prime is the vector version of -(u_complex)^(-1)
        % where u_complex = u(1) + 1i*u(2);
        u_prime = sym('u_prime',[1,2]);
        u_prime(1) = -u(1)/(u(1)^2 + u(2)^2);
        u_prime(2) = u(2)/(u(1)^2 + u(2)^2);
        u_prime_complex = u_prime(1) + 1i*u_prime(2);
        assert(simplify(-u_complex^(-1) - u_prime_complex)==0)

    elseif length(u) == 4
        % This only happens in the non-quasisplit case
        % Not sure yet what w_alpha(u) should look like in this case
        % INCOMPLETE

    end

    mat = X_SU(MatrixSize, Root_System, FormMatrix, alpha, u)* ...
        X_SU(MatrixSize, Root_System, FormMatrix, -alpha, u_prime)*...
        X_SU(MatrixSize, Root_System, FormMatrix, alpha, u);
end