function value = HomDefectCoefficientSU(MatrixSize,Root_System,Form,alpha,u,v)
    n = MatrixSize;
    q = Form.Index;
    eps = Form.Epsilon;
    Phi = Root_System;
    P = Form.PrimitiveElement;
    assert(Root_System.IsRoot(2*alpha));
    assert(abs(sum(alpha))==1);
    assert(RootSpaceDimensionSU(n,Phi,alpha)==2*(n-2*q));
    assert(RootSpaceDimensionSU(n,Phi,2*alpha)==1);

    if n==2*q
        % Quasisplit case
        value = 0;
    elseif n > 2*q
        % Non quasipslit case
        vec_C = Form.AnisotropicPartVector;
        assert(length(vec_C)==n-2*q)

        % u and v are vectors of length 2*(n-2q) of "real" quantities
        % Convert u and v to vectors of length (n-2q) of "complex" quantities
        u_complex = sym(zeros(1,n-2*q));
        u_bar = sym(zeros(1,n-2*q));
        v_complex = sym(zeros(1,n-2*q));
        v_bar = sym(zeros(1,n-2*q));
        for j=1:n-2*q
            u_complex(j) = u(2*j-1) + u(2*j)*P;
            u_bar(j) = u(2*j-1) - u(2*j)*P;
            v_complex(j) = v(2*j-1) + v(2*j)*P;
            v_bar(j) = v(2*j-1) - v(2*j)*P;
        end

        % Finally compute the value of the coefficient which is the whole purpose of this script
        value = 0;
        for j=1:n-2*q
            value = value - vec_C(j)*(u_bar(j)*v_complex(j)- u_complex(j)*v_bar(j))/(2*P);
        end

        if eps == -1
           value = (-1)*sum(alpha)*P*value;
        end

    else
        % This should be impossible
        assert(false,'The Witt index of a form is exceeding half the dimension.')
    end

    assert(length(value)==RootSpaceDimensionSU(n,Phi,2*alpha))

end