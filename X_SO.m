function mat = X_SO(MatrixSize, Root_System, FormMatrix, alpha, v)
%X_SO Summary of this function goes here
%   Detailed explanation goes here
    assert(Root_System.IsRoot(alpha));
    mat = expm(LieX_SO(MatrixSize, Root_System, FormMatrix, alpha, v));
end

