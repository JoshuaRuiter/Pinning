function my_vector = uncomplexify(scalar,PrimitiveElement)
    assert(length(scalar)==1);

    real_part = 0;
    imag_part = 0;

    % get the list of coefficients of scalar, viewed as a polynomial in P,
    % in descending order
    coefficient_list = coeffs(scalar,PrimitiveElement,'all');
    degree = polynomialDegree(scalar,PrimitiveElement);
    virtual_degree = degree;

    if mod(degree,2)==0
        % pad the coefficient list with an extra zero at the start
        coefficient_list = [0,coefficient_list];
        virtual_degree = virtual_degree + 1;
    end

    % Now there should always be an even length list of coefficients
    assert(mod(length(coefficient_list),2)==0);

    for m=1:length(coefficient_list)
        if mod(m,2)==0
            real_part = real_part + coefficient_list(m)*PrimitiveElement^(virtual_degree-m+1);
        else
            imag_part = imag_part + coefficient_list(m)*PrimitiveElement^(virtual_degree-m);
        end
    end
    my_vector = [real_part,imag_part];

end