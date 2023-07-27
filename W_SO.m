function mat = W_SO(MatrixSize,Root_System,Form,alpha,u)
    % Return the Weyl group element w_alpha(u) associated with alpha
    % Note that this is undefined when u=0, so maybe there should be
    % something like assert(u~=0)
    
    assert(dot(alpha,alpha)==1 || dot(alpha,alpha)==2)
    assert(strcmpi(Root_System.Type,'B'))
    assert(length(u)==RootSpaceDimensionSO(MatrixSize,Root_System,alpha))

    n = MatrixSize;
    q = Root_System.Rank;
    vec_C = Form.AnisotropicPartVector;
    
    if dot(alpha,alpha)==2
        % alpha is a long root
        % so the root space is one-dimensional
        assert(length(u)==1)

        if abs(sum(alpha))==2
            v = u^(-1);
        elseif abs(sum(alpha))==0
            v = -u^(-1);
        else
            % This should be impossible
            assert(false)
        end

    elseif dot(alpha,alpha)==1
        % alpha is a short root
        assert(length(u)==n-2*q)
        Qu = sym(0);
        for i=1:n-2*q
            Qu = Qu+vec_C(i)*u(i)^2;
        end
        v = 2*u/Qu;
    end

    mat = X_SO(MatrixSize,Root_System,Form,alpha,u)*...
        X_SO(MatrixSize,Root_System,Form,-alpha,v)*...
        X_SO(MatrixSize,Root_System,Form,alpha,u);

    % Why does W_SO have this form, and how was v determined?
    % Let w = W_SO(alpha,u)
    % If you assume that 
    %   (1) w^2 = identity
    %   (2) w is a product of the form X_alpha(u)*X_{-alpha}(v)*X_alpha(u)
    % then it must be the case that w has the form above.

    % Furthermore, it should also be the case that
    %   (3) w normalizes the torus S, meaning that
    %       w*s*w^(-1) is in S for any s in S
    % and this condition is met by the w as defined above.

    % Unforunately, this seems not to behave as I expect when computing a
    % conjugation of the form
    % w*X_beta(v)*w^(-1)
    % It is definitely X_{sig_alpha(beta)} ( something )
    % but that something is not +/- v, rather that something is
    % a more complicated linear transformation of v,
    % albeit a linear transformation which when written as a matrix
    % does seem to square to the identity

end