function mat = W_SU(MatrixSize,Root_System,FormMatrix,alpha,u)
    % Return the Weyl group element w_alpha(u) associated with alpha
    % Note that this is undefined when u is the zero vector

    n = MatrixSize;
    q = Root_System.Rank;
    dim_V_alpha = RootSpaceDimensionSU(n,Root_System,alpha);

    assert(Root_System.VectorLength == n)
    assert(length(FormMatrix) == n)
    assert(length(alpha) == n)
    assert(Root_System.IsRoot(alpha))
    assert(length(u) == dim_V_alpha)

    switch dot(alpha,alpha)
        
        case 4
            % alpha is a long root
            % u should have length 1
            assert(dim_V_alpha==1)
            assert(length(u)==1)
            u_prime = -u^(-1);

        case 2
            % alpha is a medium root
            % u should have length 2
            assert(dim_V_alpha==2)
            assert(length(u)==2)

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

            % Validating that u_prime is the vector version of -(u_complex)^(-1)
            u_prime_complex = u_prime(1) + 1i*u_prime(2);
            assert(simplify(-u_complex^(-1) - u_prime_complex)==0)

        case 1
            % alpha is a short root
            % u should have length n-2*q
            assert(dim_V_alpha==n-2*q)

            % This only happens in the non-quasisplit case
            % I don't know yet what w_alpha(u) should look like in this case
            % INCOMPLETE
            u_prime = zeros([1,n-2*q]); % placeholder

    end

    mat = X_SU(MatrixSize, Root_System, FormMatrix, alpha, u)* ...
        X_SU(MatrixSize, Root_System, FormMatrix, -alpha, u_prime)*...
        X_SU(MatrixSize, Root_System, FormMatrix, alpha, u);
end