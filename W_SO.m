function mat = W_SO(MatrixSize,Root_System,FormMatrix,alpha,u)
    % Return the Weyl group element w_alpha(u) associated with alpha
    % Note that this is undefined when u=0, so maybe there should be
    % something like assert(u~=0)
    
    % Extract anisotropic part of form matrix
    c1 = FormMatrix(2,2);
    c2 = FormMatrix(3,3);

    % Extract values from the vector u
    assert(length(u)==2);
    u1 = u(1);
    u2 = u(2);

    % Compute the values of a new vector v
    Qu = c1*u1^2+c2*u2^2;
    v1 = 2*u1/Qu;
    v2 = 2*u2/Qu;
    v = [v1, v2];

    mat = X_SO(MatrixSize,Root_System,FormMatrix,alpha,u)*...
        X_SO(MatrixSize,Root_System,FormMatrix,-alpha,v)*...
        X_SO(MatrixSize,Root_System,FormMatrix,alpha,u);

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