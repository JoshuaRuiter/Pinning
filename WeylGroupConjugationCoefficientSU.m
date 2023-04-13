function coeff = WeylGroupConjugationCoefficientSU(MatrixSize,Root_System,FormMatrix,alpha,beta,v)
    % given two roots alpha, beta
    % compute the coefficient phi(v)
    % which is the input to X_{sig_alpha(beta)}( )
    % appearing on the right hand side of the formula
    % w_alpha(1)*X_beta(v)*w_alpha(1) = X_{sig_alpha(beta)}( phi(v) )

    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(Root_System.IsRoot(beta))
    reflected_root = RootSystem.ReflectRoot(alpha,beta);
    assert(RootSpaceDimensionSU(Root_System,beta) == length(v))
    assert(RootSpaceDimensionSU(Root_System,reflected_root) == RootSpaceDimensionSU(Root_System,beta))

    if IsLong(beta)

        assert(length(v) == 1)

        % We can write beta = +/- 2 alpha_j
        j = find(beta~=0);
        eps_j = beta(j)/2;
        alpha_j = zeros(1,Root_System.VectorLength);
        alpha_j(j) = 1;
        assert(isequal(beta,2*eps_j*alpha_j))
       
        if IsLong(alpha)
            % alpha and beta are long
            % coeff should be +/- v in this case
            % We can write alpha = +/- 2 alpha_i
            i = find(alpha~=0);
            eps_i = alpha(i)/2;
            alpha_i = zeros(1,Root_System.VectorLength);
            alpha_i(i) = 1;
            assert(isequal(alpha,2*eps_i*alpha_i))
            
            if i==j
                coeff = -v;
            else
                coeff = v;
            end

        elseif IsMedium(alpha)
            % alpha medium, beta is long

            % We can write alpha = +/- alpha_i +/- alpha_k
            ik = find(alpha~=0);
            i = ik(1);
            k = ik(2);
            eps_i = alpha(i);
            eps_k = alpha(k);
            alpha_i = zeros(1,Root_System.VectorLength);
            alpha_i(i) = 1;
            alpha_k = zeros(1,Root_System.VectorLength);
            alpha_k(k) = 1;
            assert(isequal(alpha,eps_i*alpha_i + eps_k*alpha_k))

            coeff = v*(-1)*eps_i*eps_k;

        else
            % alpha is short
            % only happens in non-quasisplit case
            % INCOMPLETE
            assert(false,'Non-quasisplit case incomplete for Weyl group conjugation coefficients.')
        end

    elseif IsMedium(beta)

        assert(length(v)==2)
        if isa(v,'sym')
            assumeAlso(v(1),'real')
            assumeAlso(v(2),'real')
        end
        v_complex = v(1)+1i*v(2);
        v_complex_conjugate = conj(v_complex);
        v_complex_conjugate_vector = [real(v_complex_conjugate),imag(v_complex_conjugate)];

        % We can make beta = +/- alpha_j +/- alpha_k
        jk = find(beta~=0);
        j = jk(1);
        k = jk(2);
        eps_j = beta(j);
        eps_k = beta(k);
        alpha_j = zeros(1,Root_System.VectorLength);
        alpha_j(j) = 1;
        alpha_k = zeros(1,Root_System.VectorLength);
        alpha_k(k) = 1;
        assert(isequal(beta,eps_j*alpha_j + eps_k*alpha_k))

        if IsLong(alpha)
            % alpha is long and beta is medium
             % alpha and beta are long
            % coeff should be +/- v in this case
            % We can write alpha = +/- 2 alpha_i
            i = find(alpha~=0);
            eps_i = alpha(i)/2;
            alpha_i = zeros(1,Root_System.VectorLength);
            alpha_i(i) = 1;
            assert(isequal(alpha,2*eps_i*alpha_i))

            if i==min(j,k)
                coeff = v_complex_conjugate_vector;
            else
                coeff = v;
            end
            coeff = coeff*eps_i*eps_j*eps_k;
    
        elseif IsMedium(alpha)
            % alpha and beta are medium length

            % We can make alpha = eps_p*alpha_p + eps_q*alpha_q
            pq = find(alpha~=0);
            p = pq(1);
            q = pq(2);
            eps_p = alpha(p);
            eps_q = alpha(q);
            alpha_p = zeros(1,Root_System.VectorLength);
            alpha_p(p) = 1;
            alpha_q = zeros(1,Root_System.VectorLength);
            alpha_q(q) = 1;
            assert(isequal(alpha,eps_p*alpha_p + eps_q*alpha_q))

            if eps_j*eps_k*eps_p*eps_q == 1
                coeff = -v;        
            else
                coeff = -v_complex_conjugate_vector;
            end

        else
            % alpha is short
            % only happens in non-quasisplit case
            % INCOMPLETE
            assert(false,'Non-quasisplit case incomplete for Weyl group conjugation coefficients.')
        end

    else 
        % beta is short
        % only happens in non-quasisplit case
        % INCOMPLETE
        assert(false,'Non-quasisplit case incomplete for Weyl group conjugation coefficients.')
    end

    assert(length(coeff)==RootSpaceDimensionSU(Root_System,reflected_root))

end

function bool = IsShort(alpha)
    bool = (dot(alpha,alpha)==1);
end
function bool = IsMedium(alpha)
    bool = (dot(alpha,alpha)==2);
end
function bool = IsLong(alpha)
    bool = (dot(alpha,alpha)==4);
end