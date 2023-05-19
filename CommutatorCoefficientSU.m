function N = CommutatorCoefficientSU(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % Compute the commutator coefficient N_{ij}^{alpha beta}(u,v)
    % for the group SU_{n,q}(L,h)
    % That is, N_{ij}^{alpha beta}(u,v) is the input to
    % X_{i*alpha+j*beta} appearing on the right hand side of the formula
    % for the commutator [X_alpha(u), X_beta(v)]

    n = MatrixSize;
    q = Root_System.Rank;
    P = Form.PrimitiveElement;

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
            assert(IsMedium(alpha+beta))
            assert(Root_System.IsRoot(2*alpha+beta))
            assert(length(combos)==2)
            N = Commutator_Coefficient_Medium_Long_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v);

        elseif IsLong(alpha) && IsMedium(beta)
            % alpha and beta are different lengths, and their sum is a root
            % so it must be medium length, and there should be two linear
            % combinations: alpha+beta, and alpha+2*beta
            % alpha+2*beta is long
            assert(IsMedium(alpha+beta))
            assert(Root_System.IsRoot(alpha+2*beta))
            assert(IsLong(alpha+2*beta))
            assert(length(combos)==2)

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

    assert(RootSpaceDimensionSU())

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
    assert(m~=n)
    pq = find(beta);
    assert(length(pq)==2)
    p = pq(1);
    q = pq(2);
    assert(p~=q)

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
            u_complex = u(1) + u(2)*P;
            v_complex = v(1) + v(2)*P;
            
            u_prime = u_complex;
            if alpha_sign*beta_sign == 1
                u_prime = conjugate(u_complex,P);
            end

            N = beta_sign*Tr(u_prime*v_complex);

        elseif abs(sum(alpha))==2

            assert(sum(beta)==0)

            % In this case, just reverse the roles of alpha and beta
            % and re-run the calculation. Reversing the roles of alpha 
            % and beta in this way adds a negative sign.
            N = -Commutator_Coefficient_Medium_Medium_Quasisplit(Root_System,beta,alpha,j,i,v,u);

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

        % In this case, there must be exactly 3 distinct indices among (m,n,p,q)
        % Furthermore, the matching indices have opposite signs
        % We want to keep track of which entry is nonzero in both roots
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
            assert(false);
        end

        % In this case, N should be a vector of length 2
        % The commutator coefficient depends on the order of r, s, and t
        if alpha(r)==-alpha(s) && beta(r)==-beta(t)
            % The entries of each root have opposite signs
            % Here, the commutator coefficient is just the sign of beta(r)
            N = beta(r)*complexProduct(u,v);
        elseif r < s && r < t
            if alpha(r)==alpha(s) && beta(r)==beta(t)
                N = alpha(r)*vectorConjugate(complexProduct(u,v));
            elseif (t-s)*beta(r)*beta(t) > 0
                % Entries with the same sign are spaced farther apart than
                % entries with the opposite sign.
                N = beta(r)*complexProduct(u,v);
            else
                N = beta(r)*vectorConjugate(complexProduct(u,v));
            end
        elseif r > s && r > t
            if alpha(r)==alpha(s) && beta(r)==beta(t)
                N = alpha(r)*complexProduct(u,v);
            elseif s < t
                N = beta(r)*complexProduct(u,vectorConjugate(v));
            else 
                N = beta(r)*complexProduct(vectorConjugate(u),v);
            end
        else
            % In this case, r is between s and t
            if (t-s)*beta(r)*beta(t)>0 && alpha(s)==beta(t)
                N = beta(r)*complexProduct(u,v);
            elseif s < t
                N = beta(r)*alpha(s)*beta(t)*complexProduct(u,vectorConjugate(v));
            else 
                N = beta(r)*alpha(s)*beta(t)*complexProduct(vectorConjugate(u),v);
            end
        end

    else
        % This should be impossible
        assert(false,'The sum of two medium roots is a root but neither medium nor long.')
    end

    assert(length(N) == RootSpaceDimensionSU(n,Root_System,i*alpha+j*beta))
end
function N = Commutator_Coefficient_Medium_Long_Quasisplit(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    
    P = Form.PrimitiveElement;

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

    % We can always write beta in the form -2*epsilon*alpha_p
    % for some p and some +/- sign epsilon
    p = find(beta~=0);
    assert(length(p)==1)
    epsilon = 1;
    if sum(beta)==2
        epsilon = -1;
    end
    alpha_p = zeros(1,Root_System.VectorLength);
    alpha_p(p) = 1;
    assert(isequal(beta,-2*epsilon*alpha_p))

    % We can always write alpha in the form epsilon*alpha_p + omega*alpha_q
    % for some independent signs epsilon and omega
    % and some positive integer indices p and q
    pq = find(alpha~=0);
    assert(length(pq)==2);
    if pq(1)==p
        q = pq(2);
    else
        q = pq(1);
    end
    assert(p~=q)
    omega = alpha(q);
    alpha_q = zeros(1,Root_System.VectorLength);
    alpha_q(q) = 1;
    assert(isequal(alpha, epsilon*alpha_p + omega*alpha_q))

    % Finally, do the calculations for N
    if i==1 && j==1
        % alpha is medium, beta is long, 
        % and the sum under consideration is alpha+beta
        assert(IsMedium(i*alpha+j*beta))
        
        % c_pq trackes the order of p and q
        % if p < q, then c_pq = 1
        % if p > q, then c_pq = -1
        c_pq = 1;
        if p > q
            c_pq = -1;
        end

        % c_pq determines if u should be conjugated or not
        assert(length(u)==2)
        u_complex = u(1)+u(2)*P;
        if c_pq == 1
            %u(2) = -u(2);
            u_prime = conjugate(u_complex,P);
        else
            u_prime = u;
        end

        % omega is a scalar +/-1
        % u is a vector of length 2
        % v is a scalar
        % so the multiplication below is well defined
        N_complex = omega*v*u_prime;

        % Now convert N to vector version
        N_real_part = subs(N_complex,P,0);
        N_imag_part = (N_complex - N_real_part)/P;
        N = [N_real_part,N_imag_part];

    elseif i==2 && j==1
        % alpha is medium, beta  is long,
        % and the sum under consideration is 2*alpha+beta
        assert(IsLong(i*alpha+j*beta))
        u_complex = u(1) + u(2)*P;
        u_squared_norm = u_complex*conjugate(u_complex,P);
        N = (-1)*epsilon*omega*v*u_squared_norm;

    else
        % This should be impossible
        assert(false,"A medium and long root in C_n have an unexpected linear combination.")
    end

    alpha
    beta
    i
    j
    N

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


function uv = complexProduct(u,v)
    uv = [u(1)*v(1)-u(2)*v(2),u(1)*v(2)+u(2)*v(1)];
end
function uOut = vectorConjugate(uIn)
    uOut = [uIn(1),-uIn(2)];
end