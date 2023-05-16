function RunSUTests(n,q,eps)

    % Use the PinnedGroup class to verify a pinning of SU_{n,q}

    MatrixSize = n;
    RootSystemRank = q;
    syms d;
    PrimitiveElement = sqrt(d);

    if eps == 1
        NameString = strcat('special unitary group of size',{' '},num2str(n),...
            {' associated to a hermitian form with Witt index '},num2str(q));
    elseif eps == -1
        NameString = strcat('special unitary group of size',{' '},num2str(n),...
            {' associated to a skew-hermitian form with Witt index '},num2str(q));
    end

    % Setting up the root system
    if n==2*q
        % Quasisplit case, use C_q root system
        Type = 'C';
    else
        % Non-quasisplit case, use BC_q root system
        Type = 'BC';
    end
    root_system = RootSystem(Type,RootSystemRank,MatrixSize);

    % Building the matrix of the (skew-)hermitian form
    if n > 2*q
        vec_C = sym('c',[1,n-2*q]);
        if eps == 1
            % In the hermitian case, entries of C must be "purely real"
            for i=1:n-2*q
                vec_C(i) = QuadraticExtensionElement(vec_C(i),0,PrimitiveElement);
            end
        elseif eps == -1
            % In the skew-hermitian case, entries of C must be "purely imaginary"
            for i=1:n-2*q
                vec_C(i) = QuadraticExtensionElement(0,vec_C(i),PrimitiveElement);
            end
        end
    else
        vec_C = [];
    end
    if eps == 1
        label = 'hermitian';
    elseif eps == -1
        label = 'skew-hermitian';
    end
    Form = NIForm(n,q,eps,vec_C,PrimitiveElement,label);
    Form.Matrix % Remove later, for debugging

    RootSpaceDimension = @RootSpaceDimensionSU;
    RootSpaceMap = @LieX_SU;
    RootSubgroupMap = @X_SU;
    WeylGroupMap = @W_SU;
    GenericTorusElementMap = @GenericTorusElementSU;
    IsGroupElement = @IsInSU;
    IsTorusElement = @IsTorusElementSU;
    IsLieAlgebraElement = @IsIn_little_su;
    CommutatorCoefficientMap = @CommutatorCoefficientSU;
    WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSU;

    SU_n_q = PinnedGroup(NameString,MatrixSize,root_system,Form,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        CommutatorCoefficientMap,WeylGroupCoefficientMap);

    SU_n_q.RunTests();

end

function bool = IsTorusElementSU(MatrixSize, RootSystemRank, MatrixToTest)

    q = RootSystemRank;
    n = MatrixSize;

    % Check for correct matrix dimensions
    bool = (length(MatrixToTest)==n) && SymbolicIsDiag(MatrixToTest);

    % Check that the first q entries and the next q diagonal entries are
    % respective inverses
    for i=1:q
        bool = bool && (MatrixToTest(i,i) == MatrixToTest(q+i,q+i)^(-1));
    end

    % Check that the last block is the identity
    for i=2*q+1:n
        bool = bool && (MatrixToTest(i,i) == 1);
    end
end
function bool = IsIn_little_su(MatrixSize,MatrixToTest,Form)
    zero_quad = QuadraticExtensionElement(0,0,Form.PrimitiveElement);

    % Compute trace of the matrix
    my_trace = zero_quad;
    for i=1:MatrixSize
        my_trace = my_trace + MatrixToTest(i,i);
    end

    bool = (length(MatrixToTest)==MatrixSize && ...
        eq(my_trace,zero_quad) && ...
        SymbolicIsEqual(transpose(conj(MatrixToTest))*Form.Matrix,-Form.Matrix*MatrixToTest));
end
function bool = IsInSU(MatrixSize,MatrixToTest,Form)
    bool = (length(MatrixToTest)==MatrixSize && ...
        det(MatrixToTest) == 1 && ...
        SymbolicIsEqual(ctranspose(MatrixToTest)*Form.Matrix*MatrixToTest,Form.Matrix));
end
function mat = GenericTorusElementSU(MatrixSize, RootSystemRank, DiagonalValues)
    % given a vector [t_1, t_2, t_3, ... t_q] of length q = RootSystemRank
    % build the symbolic diagonal matrix 
    % diag(t_1, ... t_q, t_1^(-1), ... t_q^(-1), 1, ... 1)
    % Note: t1, ... t_q, etc. should be "purely real"
    
    assert(length(DiagonalValues)==RootSystemRank);

    % Add an assumption that each entry of DiagonalValues is real
    for i=1:length(DiagonalValues)
        t_i = DiagonalValues(i);
        assumeAlso(t_i, 'real');
    end

    mat = sym(zeros(MatrixSize));
    for i=1:RootSystemRank
        mat(i,i) = DiagonalValues(i);
        mat(RootSystemRank+i,RootSystemRank+i) = DiagonalValues(i)^(-1);
    end
    for i=2*RootSystemRank+1:MatrixSize
        mat(i,i) = 1;
    end
end

