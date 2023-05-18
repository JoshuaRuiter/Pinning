function value = HomDefectCoefficientSU(MatrixSize,Root_System,Form,alpha,u,v)
   
    assert(Root_System.IsRoot(2*alpha));

    n = MatrixSize;
    q = Form.Index;
    Phi = Root_System;
    P = Form.PrimitiveElement;
    dim_V_2alpha = RootSpaceDimensionSU(n, Phi, alpha);

    if n==2*q
        % Quasisplit case
        value = sym(zeros(1,dim_V_2alpha));
    elseif n > 2*q
        % In terms of the quadratic field extension, output should be
        % value = conjugate(v,P)*u - conjugate(u)*v;
        % However, we have to write this out in terms of vectors
        
        assert(RootSpaceDimensionSU(n,Phi,alpha)==2);
        assert(length(u)==2);

        u_complex = u(1)+u(2)*P;
        v_complex = v(1)+v(2)*P;
        u_bar = conjugate(u,P);
        v_bar = conjugate(v,P);

        value_complex = u_complex*v_bar + u_bar*v_complex;
        real_part = subs(value_complex,P,0);
        imag_part = (value_complex - real_part)/P;
        value = [real_part,imag_part];

    else
        % This should be impossible
        assert(false,'The Witt index of a form is exceeding half the dimension.')
    end

end