function RunSUTests(n,q,eps)

    % Use the PinnedGroup class to verify a pinning of SU_{n,q}

    MatrixSize = n;
    RootSystemRank = q;
    % syms d;
    % PrimitiveElement = sqrt(d);
    syms P;
    PrimitiveElement = P;

    if eps == 1
        P_eps = P;
    else % eps == -1
        P_eps = 1;
    end


    if eps == 1
        NameString = strcat('special unitary group \nof size',{' '},num2str(n),...
            {' associated to a hermitian form with Witt index '},num2str(q));
    elseif eps == -1
        NameString = strcat('special unitary group \nof size',{' '},num2str(n),...
            {' associated to a skew-hermitian form with Witt index '},num2str(q));
    end

    if n==2*q
        NameString = strcat('quasisplit',{' '},NameString);
    elseif n > 2*q
        NameString = strcat('non-quasisplit',{' '},NameString);
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

        % In the hermitian case, entries of C must be "purely real"
        % In the skew-hermitian case, entries of C must be "purely imaginary"
        if eps == -1
            vec_C = vec_C*P;
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

    RootSpaceDimension = @RootSpaceDimensionSU;
    RootSpaceMap = @LieX_SU;
    RootSubgroupMap = @X_SU;
    WeylGroupMap = @W_SU;
    GenericTorusElementMap = @GenericTorusElementSU;
    IsGroupElement = @IsInSU;
    IsTorusElement = @IsTorusElementSU;
    IsLieAlgebraElement = @IsIn_little_su;
    HomDefectCoefficientSU = @HomDefectCoefficientSU;
    CommutatorCoefficientMap = @CommutatorCoefficientSU;
    WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSU;

    SU_n_q = PinnedGroup(NameString,MatrixSize,root_system,Form,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        HomDefectCoefficientSU, CommutatorCoefficientMap,...
        WeylGroupCoefficientMap);

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
    conjugated_matrix = conjugate(MatrixToTest,Form.PrimitiveElement);
    conj_transpose = transpose(conjugated_matrix);
    bool = (length(MatrixToTest)==MatrixSize && ...
        trace(MatrixToTest)==0 && ...
        isequal(conj_transpose*Form.Matrix,-Form.Matrix*MatrixToTest));
end
function bool = IsInSU(MatrixSize,MatrixToTest,Form)
    conjugated_matrix = conjugate(MatrixToTest,Form.PrimitiveElement);
    conj_transpose = transpose(conjugated_matrix);
    bool = (length(MatrixToTest)==MatrixSize && ...
        det(MatrixToTest) == 1 && ...
        SymbolicIsEqual(conj_transpose*Form.Matrix*MatrixToTest,Form.Matrix));
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

