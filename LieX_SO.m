function mat = LieX_SO(MatrixSize, Root_System, FormMatrix, alpha, v)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    assert(Root_System.IsRoot(alpha));
    mat = SymbolicZeros(MatrixSize);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specific to SO_4_1
    c1 = FormMatrix(2,2);
    c2 = FormMatrix(3,3);
    if isequal(alpha,[1,0,0,0])
        mat(1,2) = -c1*v(1);
        mat(1,3) = -c2*v(2);
        mat(2,4) = v(1);
        mat(3,4) = v(2);
    elseif isequal(alpha,[-1,0,0,0])
        mat(2,1) = v(1);
        mat(3,1) = v(2);
        mat(4,2) = -c1*v(1);
        mat(4,3) = -c2*v(2);
    else
        assert(false,'LieX_SO not written for cases other than n=4,q=1 yet');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
end