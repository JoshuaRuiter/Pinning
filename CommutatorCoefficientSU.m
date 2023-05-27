function N = CommutatorCoefficientSU(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % Compute the commutator coefficient N_{ij}^{alpha beta}(u,v)
    % for the group SU_{n,q}(L,h)
    % That is, N_{ij}^{alpha beta}(u,v) is the input to
    % X_{i*alpha+j*beta} appearing on the right hand side of the formula
    % for the commutator [X_alpha(u), X_beta(v)]

    n = MatrixSize;

    % Validating inputs
    assert(i>=1);
    assert(j>=1);
    assert(strcmp(Root_System.Type,'C') || strcmp(Root_System.Type,'BC'))
    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(Root_System.IsRoot(beta))
    assert(length(u) == RootSpaceDimensionSU(n,Root_System,alpha))
    assert(length(v) == RootSpaceDimensionSU(n,Root_System,beta))

    % If alpha and beta are proportional roots, then the commutator formula
    % does not apply, and we should throw an error when trying to compute 
    % the commutator coefficient
    assert(~RootSystem.IsProportionalRoot(alpha,beta));

    % If i*alpha + j*beta is not a root, then the coefficient should be
    % zero, and there is nothing else to do
    if ~Root_System.IsRoot(i*alpha+j*beta)
        N = 0;
        return;
    end

    combos = Root_System.LinearCombos(alpha,beta);
    assert(Root_System.IsRoot(alpha+beta))
    assert(length(combos)>=1)

    if Root_System.Type == 'C'
        % Quasisplit case
        % We know that alpha and beta are roots in C_q, (q = rank of root system) 
        % and that alpha + beta is a root
        assert(Root_System.IsRoot(alpha+beta))
        assert(Root_System.IsRoot(i*alpha+j*beta))

        % In this case, there are only two root lengths: medium and long
        % The possibilities are:
        %   (1) alpha, beta are both medium, and i=j=1, and
        %       (1a) and alpha+beta is medium
        %       (1b) and alpha+beta is long
        %   (2) one is medium, one is long, and alpha+beta is medium.
        %       In this case, it is possible to have i=j=1, or
        %       one of i,j is 1 and the other is 2.
        %
        % Notably, it is not possible for both alpha and beta to be long
        % and for alpha+beta to be a root
        assert(i==1 || i==2)
        assert(j==1 || j==2)
        assert(i+j <= 3)
        assert(length(combos) <= 2)

        if IsMedium(alpha) && IsMedium(beta)
            % both alpha and beta are medium length,
            % so the sum alpha+beta could be long or medium
            assert(i==1)
            assert(j==1)
            assert(IsMedium(alpha+beta) || IsLong(alpha+beta))
            assert(length(combos)==1)

            N = Commutator_Coefficient_Medium_Medium_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsMedium(alpha) && IsLong(beta)
            % alpha and beta are different lengths, and their sum is a root
            % so it must be medium length, and there should be two linear
            % combinations: alpha+beta, and 2*alpha+beta
            assert(length(combos)==2)
            assert(Root_System.IsRoot(alpha+beta))
            assert(Root_System.IsRoot(2*alpha+beta))
            assert(IsMedium(alpha+beta))
            assert(IsLong(2*alpha+beta))

            N = Commutator_Coefficient_Medium_Long_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsLong(alpha) && IsMedium(beta)
            % alpha and beta are different lengths, and their sum is a root
            % so it must be medium length, and there should be two linear
            % combinations: alpha+beta, and alpha+2*beta
            % alpha+2*beta is long
            assert(length(combos)==2)
            assert(Root_System.IsRoot(alpha+beta))
            assert(Root_System.IsRoot(alpha+2*beta))
            assert(IsMedium(alpha+beta))
            assert(IsLong(alpha+2*beta))

            % In this case, just reverse the roles of alpha and beta
            % which also reverses the roles of i and j, and u and v
            % and makes the coefficient negative
            N = -Commutator_Coefficient_Medium_Long_Quasisplit(MatrixSize,Root_System,Form,beta,alpha,j,i,v,u);

        else
            % This should be impossible
            assert(false,'The sum of two roots is a root but not of the expected length.')
        end

    elseif Root_System.Type == 'BC'
        % Non quasisplit case
        % We know that alpha and beta are roots in BC_q, (q = rank of root system)
        % and that alpha + beta is a root
        assert(Root_System.IsRoot(i*alpha+j*beta))

        % In this case, there are three root lengths: short, medium, and long
        % There are many possibilities:
        %   (1) alpha, beta are both medium, and
        %       (1a) and alpha+beta is medium
        %       (1b) and alpha+beta is long
        %   (2) one is medium, one is long, and alpha+beta is medium
        %   (3) alpha, beta are both short, and alpha+beta is medium
        %   (4) one is short, one is medium, and alpha+beta short
        % 
        % Note that we do NOT need to consider the situation of a short
        % root plus itself to get a long root, since in that case alpha and
        % beta would be proportional, in which case the commutator formula
        % does not apply.
        % 
        % Cases which are impossible
        %   (a) Both roots are long
        %   (b) One is short and one is long

        if IsMedium(alpha) && IsMedium(beta)
            % both alpha and beta are medium length,
            % so the sum alpha+beta could be long or medium
            assert(IsLong(alpha+beta) || IsMedum(alpha+beta))

            % Hopefully this case ends up being basically identical to the
            % case with the same root lengths in the quasisplit case
            % In that case, use
            % N = Commutator_Coefficient_Quasisplit_Two_Medium(alpha,beta,i,j,u,v);

            % Otherwise, use
            N = Commutator_Coefficient_Medium_Medium(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsMedium(alpha) && IsLong(beta)
            % alpha and beta are medium and long, and their sum is a root
            % so it must be medium length
            assert(IsMedium(alpha+beta))

            % Hopefully this case ends up being basically identical to the
            % case with the same root lengths in the quasisplit case
            % In that case, use
            % N = Commutator_Coefficient_Quasisplit_Medium_Long(alpha,beta,i,j,u,v);

            % Otherwise, use
            N = Commutator_Coefficient_Medium_Long(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsLong(alpha) && IsMedium(beta)
            % alpha and beta are medium and long, and their sum is a root
            % so it must be medium length
            assert(IsMedium(alpha+beta))

            % In this case, just reverse the roles of alpha and beta
            % which also reverses the roles of i and j, and u and v
            % and makes the coefficient negative of whatever it would be in
            % the previous case

            % If it works out to be the same as the quasisplit case, use
            % N = Commutator_Coefficient_Quasisplit_Medium_Long(beta,alpha,j,i,v,u);

            % Otherwise, use
            N = -Commutator_Coefficient_Medium_Long(MatrixSize,Root_System,Form,beta,alpha,j,i,v,u);

        elseif IsShort(alpha) && IsShort(beta)
            % alpha and beta are both short, and their sum is a root
            % and they are not proportional, in particular alpha is not
            % equal to beta, so their sum is medium
            assert(IsMedium(alpha+beta))

            % This is unlike anything in the quasisplit case
            N = Commutator_Coefficient_Short_Short(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsShort(alpha) && IsMedium(beta)
            % alpha is short and beta is medium, and their sum is a root
            % so it must be short
            assert(IsShort(alpha+beta))

            % This is unlike anything in the quasisplit case
            N = Commutator_Coefficient_Short_Medium(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsMedium(alpha) && IsShort(beta)
            % alpha is short and beta is medium, and their sum is a root
            % so it must be short
            assert(IsShort(alpha+beta))

            % Just reuse the previous case, but reverse the roles of alpha
            % and beta, i and j, u and v, and then make the coefficient the
            % negative of what it would have been
            N = -Commutator_Coefficient_Short_Medium(MatrixSize,Root_System,Form,beta,alpha,j,i,v,u);

        else
            % This should be impossible
            assert(false,'The sum of two roots is a root but not of the expected length.')
        end
    end
    
    % Validating output
    % N should always be a vector whose length is the dimension of the root
    % space for i*alpha+j*beta
    assert(length(N) == RootSpaceDimensionSU(n,Root_System,i*alpha+j*beta))

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

function N = Commutator_Coefficient_Medium_Medium_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % alpha and beta are both medium length
    % the sum can be long or medium

    assert(IsMedium(alpha))
    assert(IsMedium(beta))
    assert(IsLong(alpha+beta) || IsMedium(alpha+beta))

    % In this case, there should only be one integral linear combination of
    % alpha and beta which is a root, namely alpha+beta
    combos = Root_System.LinearCombos(alpha,beta);
    assert(length(combos)==1)
    assert(i==1)
    assert(j==1)

    % In this case, it must be that alpha and beta have the form
    % alpha = +/- alpha_m +/- alpha_n
    % beta = +/- alpha_p +/- alpha_q
    mn = find(alpha);
    assert(length(mn)==2)
    m = mn(1);
    n = mn(2);
    assert(m<n)
    pq = find(beta);
    assert(length(pq)==2)
    p = pq(1);
    q = pq(2);
    assert(p<q);

    alpha_m = zeros(1,Root_System.VectorLength);
    alpha_m(m) = 1;
    alpha_n = zeros(1,Root_System.VectorLength);
    alpha_n(n) = 1;
    alpha_p = zeros(1,Root_System.VectorLength);
    alpha_p(p) = 1;
    alpha_q = zeros(1,Root_System.VectorLength);
    alpha_q(q) = 1;
    eps_m = alpha(m);
    eps_n = alpha(n);
    eps_p = beta(p);
    eps_q = beta(q);
    assert(isequal(alpha,eps_m*alpha_m+eps_n*alpha_n));
    assert(isequal(beta,eps_p*alpha_p+eps_q*alpha_q));

    if IsLong(alpha+beta)
        % alpha and beta are both medium, and their sum is long
        % This is only possible if there are exactly two distinct indices among (m,n,p,q)
        % Since m~=n and p~=q, this means that the set {alpha, beta} has
        % the form
        % {+/- (alpha_a + alpha_b), +/- (alpha_c - alpha_d)}
        % So there are two cases to analyze
        %   (1) alpha = +/-(alpha_m + alpha_n), and beta = +/-(alpha_m - alpha_n)
        %   (2) alpha = +/-(alpha_m - alpha_n), and beta = +/-(alpha_m - alpha_n)

        if sum(alpha)==0

            % Out first goal is to write alpha in the form
            % alpha_sign*(alpha_m - alpha_n)
            % with m < n

            % Note that we overwrite the previous values of m and n
            m = find(alpha==1);
            n = find(alpha==-1);
            assert(m~=n)

            if m < n
                % alpha = alpha_sign(alpha_m - alpha_n) and m < n
                alpha_sign = 1;
            else % m > n
                % alpha = alpha_sign(alpha_m - alpha_n) and m < n
                % swap the values of m and n
                temp = m;
                m = n;
                n = temp;
                alpha_sign = -1;
            end

            % Now alpha = alpha_sign(alpha_m - alpha_n)
            alpha_m = zeros(1,Root_System.VectorLength);
            alpha_m(m) = 1;
            alpha_n = zeros(1,Root_System.VectorLength);
            alpha_n(n) = 1;
            assert(isequal(alpha,alpha_sign*(alpha_m-alpha_n)))

            % Our next goal is to write beta in the form 
            % beta_sign*(alpha_m + alpha_n)
            assert(abs(sum(beta))==2)
            if sum(beta)==2
                beta_sign = 1;
            else
                beta_sign = -1;
            end
            assert(isequal(beta,beta_sign*(alpha_m+alpha_n)))

            % the product of the signs on alpha and beta
            % determines whether u should be conjugated
            % Since u and v are represented as a vectors of length 2, first we
            % convert to a "complex" representation
            P = Form.PrimitiveElement;
            u_complex = complexify(u,P);
            v_complex = complexify(v,P);
            u_bar = conjugate(u_complex,P);
            v_bar = conjugate(v_complex,P);

            if Form.Epsilon == 1
                % Hermitian case
                if alpha_sign == 1
                    if beta_sign == 1
                        N = u_bar*v_complex - u_complex*v_bar;
                    else
                        assert(beta_sign == -1)
                        N = u_bar*v_bar - u_complex*v_complex;
                    end
                else
                   assert(alpha_sign == -1)
                    if beta_sign == 1
                        N = u_complex*v_complex - u_bar*v_bar;
                    else
                        assert(beta_sign == -1)
                        N = u_complex*v_bar - u_bar*v_complex;
                    end
                end
                N = N/P;
            else
                if alpha_sign == 1
                    if beta_sign == 1
                        N = u_complex*v_bar + u_bar*v_complex;
                    else
                        N = -u_complex*v_complex - u_bar*v_bar;
                    end
                else
                    if beta_sign == 1
                        N = u_complex*v_complex + u_bar*v_bar;
                    else
                        N = -u_complex*v_bar - u_bar*v_complex;
                    end
                end
            end

        elseif abs(sum(alpha))==2

            assert(sum(beta)==0)

            % In this case, just reverse the roles of alpha and beta
            % and re-run the calculation. Reversing the roles of alpha 
            % and beta in this way adds a negative sign.
            N = -Commutator_Coefficient_Medium_Medium_Quasisplit(MatrixSize,Root_System,Form,beta,alpha,j,i,v,u);

        else
            % This should be impossible
            assert(false,"A medium root in C_n has the wrong form.")
        end
    
    elseif IsMedium(alpha+beta)
        % alpha and beta are both medium, and their sum is medium
        % This is only possible if n>=6 (equivalently q>=3)

        % In this case, the N will always have the form
        % uv, -uv, conj(u) v, -conj(u) v, u conj(v), -u conj(v), etc.
        % u may be conjugated or not
        % v may be conjugated or not
        % the overall sign may be 1 or -1
        % and all three of these happen independently
        P = Form.PrimitiveElement;
        u_complex = complexify(u,P);
        v_complex = complexify(v,P);
        u_bar = conjugate(u_complex,P);
        v_bar = conjugate(v_complex,P);

        % Recall that
        % alpha = +/- alpha_m +/- alpha_n
        % beta = +/- alpha_p +/- alpha_q
        % with m < n and p < q
        assert(isequal(alpha,eps_m*alpha_m+eps_n*alpha_n));
        assert(isequal(beta,eps_p*alpha_p+eps_q*alpha_q));
        assert(m<n)
        assert(p<q)

        % In this case, there must be exactly 3 distinct indices among (m,n,p,q)
        % Furthermore, the matching indices have opposite signs
        % In other words, we can rewrite alpha and beta as
        % alpha = eps_r*alpha_r + eps_s*alpha_s
        % beta = -eps_r*alpha_r + eps_t*alpha_t
        if m == p || m == q
            r = m;
            s = n;
            if m == p
                t = q;
            else
                t = p;
            end
        elseif n == p || n == q
            r = n;
            s = m;
            if n == p
                t = q;
            else
                t = p;
            end
        else
            assert(false,'A sum of two medium roots does not have the expected form.');
        end

        % alpha = eps_r*alpha_r + eps_s*alpha_s
        % beta = -eps_r*alpha_r + eps_t*alpha_t
        alpha_r = zeros(1,Root_System.VectorLength);
        alpha_r(r) = 1;
        alpha_s = zeros(1,Root_System.VectorLength);
        alpha_s(s) = 1;
        alpha_t = zeros(1,Root_System.VectorLength);
        alpha_t(t) = 1;
        eps_r = alpha(r);
        eps_s = alpha(s);
        eps_t = beta(t);
        assert(beta(r)==-eps_r)
        assert(isequal(alpha,eps_r*alpha_r+eps_s*alpha_s))
        assert(isequal(beta,-eps_r*alpha_r+eps_t*alpha_t))

        % Placeholder
        N_complex = 0;

        r
        s
        t

        if alpha(r)==-alpha(s) && beta(r)==-beta(t)
        %if sum(alpha)==0 && sum(beta)==0
        %if eps_r == -eps_s && -eps_r == -eps_t
            % The entries of each root have opposite signs
            % This behaves essentially like a copy of SL_n
            N_complex = - eps_r*u_complex*v_complex;
            
        elseif r < s && r < t
            % r is smaller than both s and t
            if alpha(r)==alpha(s) && beta(r)==beta(t)
            %if eps_r == eps_s && -eps_r == eps_t
                %N = alpha(r)*vectorConjugate(complexProduct(u,v));
                N_complex = eps_r*u_bar*v_bar;
            elseif (t-s)*beta(r)*beta(t) > 0
                % Entries with the same sign are spaced farther apart than
                % entries with the opposite sign.
                %N = beta(r)*complexProduct(u,v);
                N_complex = beta(r)*u_complex*v_complex;
            else
                %N = beta(r)*vectorConjugate(complexProduct(u,v));
                N_complex = eps_r*u_bar*v_bar;
            end

        elseif r > s && r > t
            % r is bigger than both s and t
            if alpha(r)==alpha(s) && beta(r)==beta(t)
                %N = alpha(r)*complexProduct(u,v);
                N_complex = alpha(r)*u_complex*v_complex;
            elseif s < t
                %N = beta(r)*complexProduct(u,vectorConjugate(v));
                N_complex = beta(r)*u_complex*v_bar;
            else 
                %N = beta(r)*complexProduct(vectorConjugate(u),v);
                N_complex = beta(r)*u_bar*v_complex;
            end
            
        else
            % In this case, r is between s and t

            if (t-s)*beta(r)*beta(t)>0 && alpha(s)==beta(t)
                %N = beta(r)*complexProduct(u,v);
                N_complex = beta(r)*u_complex*v_complex;
            elseif s < t
                %N = beta(r)*alpha(s)*beta(t)*complexProduct(u,vectorConjugate(v));
                %N_complex = (-1)*beta(r)*alpha(s)*beta(t)*u_complex*v_bar;
                N_complex = beta(r)*alpha(s)*beta(t)*u_complex*v_bar;
            else 
                %N = beta(r)*alpha(s)*beta(t)*complexProduct(vectorConjugate(u),v);
                N_complex = beta(r)*alpha(s)*beta(t)*u_bar*v_complex;
            end
        end

        N = uncomplexify(N_complex,P);

    else
        % This should be impossible
        assert(false,'The sum of two medium roots is a root but neither medium nor long.')
    end

    assert(length(N) == RootSpaceDimensionSU(n,Root_System,i*alpha+j*beta))
end
function N = Commutator_Coefficient_Medium_Long_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    
    P = Form.PrimitiveElement;
    eps = Form.Epsilon;
    if eps == 1
        P_eps = P;
    else % eps == -1
        P_eps = 1;
    end

    % alpha is medium, beta is long
    assert(IsMedium(alpha))
    assert(IsLong(beta))
    assert(length(u)==2)
    assert(length(v)==1)
    assert(IsMedium(alpha+beta))

    % In this case, alpha+beta and 2*alpha+beta are roots
    assert(Root_System.IsRoot(alpha+beta))
    assert(Root_System.IsRoot(2*alpha+beta))
    combos = Root_System.LinearCombos(alpha,beta);
    assert(length(combos)==2)
    assert(i==1 || i==2)
    assert(j==1)

    % We can always write beta in the form -2*eps_s*alpha_s
    % for some p and some +/- sign eps_s
    s = find(beta~=0);
    assert(length(s)==1)
    eps_s = 1;
    if sum(beta)==2
        eps_s = -1;
    end
    alpha_s = zeros(1,Root_System.VectorLength);
    alpha_s(s) = 1;
    assert(isequal(beta,-2*eps_s*alpha_s))

    % We can always write alpha in the form eps_s*alpha_s + eps_t*alpha_t
    % for some independent signs eps_s and eps_t
    % and some positive integer indices s and t
    st = find(alpha~=0);
    assert(length(st)==2);
    if st(1)==s
        t = st(2);
    else
        t = st(1);
    end
    assert(s~=t)
    eps_t = alpha(t);
    alpha_t = zeros(1,Root_System.VectorLength);
    alpha_t(t) = 1;
    assert(isequal(alpha, eps_s*alpha_s + eps_t*alpha_t))

    % Finally, do the calculations for N
    if i==1 && j==1
        % alpha is medium, beta is long, 
        % and the sum under consideration is alpha+beta
        assert(IsMedium(i*alpha+j*beta))
        assert(RootSpaceDimensionSU(MatrixSize,Root_System,i*alpha+j*beta)==2);

        assert(length(v)==1);
        assert(length(u)==2)
        u_complex = complexify(u,P);
        
        % c_st trackes the order of p and q
        % if s < t, then c_st = 1
        % if s > t, then c_st = -1
        % c_pq determines if u should be conjugated or not
        c_st = 1;
        if s > t
            c_st = -1;
            u_prime = u_complex;
        else
            c_st = 1;
            u_prime = conjugate(u_complex,P);
        end

        if Form.Epsilon == 1
            % Hermitian case
            N_complex = (-1)*eps_s*P_eps*v*u_prime;
            if abs(sum(alpha))==2
                N_complex = c_st*N_complex;
            end
        else 
            % Form.Epsilon == -1
            % Skew hermitian case
            if abs(sum(alpha)) == 2
                N_complex = eps_s*v*u_prime;
            else
                % sum(alpha) == 0
                N_complex = (-1)*eps_s*v*u_prime;
            end
        end
        N = uncomplexify(N_complex,P);

    elseif i==2 && j==1
        % alpha is medium, beta  is long,
        % and the sum under consideration is 2*alpha+beta
        assert(IsLong(i*alpha+j*beta))
        u_complex = complexify(u,P);
        u_bar = conjugate(u_complex,P);

        if Form.Epsilon == 1
            % Hermitian case
            N = v*u_complex*u_bar;
        else
            % Form.Epsilon == -1
            % Skew hermitian case
            N = v*u_complex*u_bar;
            if abs(sum(alpha)) == 2
                N = (-1)*N;
            end
        end
    else
        % This should be impossible
        assert(false,"A medium and long root in C_n have an unexpected linear combination.")
    end

    assert(length(N)==RootSpaceDimensionSU(MatrixSize,Root_System,i*alpha+j*beta))

end
function N = Commutator_Coefficient_Medium_Medium(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % alpha and beta are both medium length
    % the sum can be long or medium
    assert(IsMedium(alpha))
    assert(IsMedium(beta))
    assert(IsLong(alpha+beta) || IsMedium(alpha+beta))

    if IsLong(alpha+beta)
        % alpha and beta are both medium, and their sum is long
        % INCOMPLETE
        N = 0;
    
    elseif IsMedium(alpha+beta)
        % alpha and beta are both medium, and their sum is medium
        % INCOMPLETE
        N = 0;

    else
        % This should be impossible
        assert(false,'The sum of two medium roots is a root but neither medium nor long.')
    end
end
function N = Commutator_Coefficient_Medium_Long(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % alpha is medium, beta is long
    assert(IsMedium(alpha))
    assert(IsLong(beta))
    assert(IsMedium(alpha+beta))

    % INCOMPLETE
    N = 0;
end
function N = Commutator_Coefficient_Short_Short(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % alpha and beta are short, 
    % so the sum must be medium
    assert(IsShort(alpha))
    assert(IsShort(beta))
    assert(IsMedium(alpha+beta))

    % INCOMPLETE
    N = 0;
end
function N = Commutator_Coefficient_Short_Medium(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % alpha is short, beta is medium
    % so their sum must be short
    assert(IsShort(alpha))
    assert(IsMedium(beta))
    assert(IsShort(alpha+beta))

    % INCOMPLETE
    N = 0;
end
