function mat = X_SU(MatrixSize, Root_System, NIForm, alpha, u)
    % Given a root alpha in Phi(SU4,S)
    % and given a variable u
    % output the associated matrix X_alpha(u) in SU_2n(L,H)
    assert(Root_System.IsRoot(alpha));
    mat = expm(LieX_SU(MatrixSize, Root_System, NIForm, alpha, u));
end