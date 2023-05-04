function dim = RootSpaceDimensionSU(Root_System,alpha)

    % The dimension is 1 for long roots, 2 for medium roots, 
    % and 2*(n-2*q) for short roots

    % dot(alpha,alpha) gives the squared length of alpha,
    % which is 1 for short roots, 2 for medium roots, and 4 for long roots

    assert(strcmpi(Root_System.Type,'BC') ...
        || strcmpi(Root_System.Type,'C'))
    assert(Root_System.IsRoot(alpha))
    n = Root_System.VectorLength;
    q = Root_System.Rank;

    switch dot(alpha,alpha)
        case 1
            assert(strcmpi(Root_System.Type,'BC'))
            dim = 2*(n-2*q);
        case 2
            dim = 2;
        case 4
            dim = 1;
        otherwise
            assert(false,'A root in BC_q or C_q has the wrong length.')
    end
end