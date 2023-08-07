function coeff = WeylGroupConjugationCoefficientSO(MatrixSize,root_system,FormMatrix,alpha,beta,v)
    % given two roots alpha, beta
    % compute the coefficient phi(v)
    % which is the input to X_{sig_alpha(beta)}( )
    % appearing on the right hand side of the formula
    % w_alpha(1)*X_beta(v)*w_alpha(1) = X_{sig_alpha(beta)}( phi(v) )

    % Extract anisotropic part of form matrix
    c1 = FormMatrix(2,2);
    c2 = FormMatrix(3,3);

    % Extract values from the vector u
    assert(length(v)==2);
    v1 = v(1);
    v2 = v(2);
    
    % OLD VERSION
%     if isequal(alpha,beta)
%         coeff1 = (2*(c1-c2)*v1 + 4*c2*v2)/(c1+c2)^2;
%         coeff2 = (-2*(c1-c2)*v2 + 4*c1*v1)/(c1+c2)^2;
%     else
%         coeff1 = (c1-c2)*v1/2 + c2*v2;
%         coeff2 = -(c1-c2)*v2/2 + c1*v1;
%     end

    % NEW VERSION
    if isequal(alpha,beta)
        coeff1 = 2*v1/c1;
        coeff2 = -2*v2/c1;
    else
        coeff1 = c1*v1/2;
        coeff2 = -c1*v2/2;
    end
    
    coeff = [coeff1,coeff2];
end
