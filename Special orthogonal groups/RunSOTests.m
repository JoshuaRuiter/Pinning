function RunSOTests(n,q)

    % Use the PinnedGroup class to verify a pinning of SO_{n,q}

    NameString = strcat('special orthogonal group of size',{' '},num2str(n),{' '},'with Witt index',{' '},num2str(q));
    MatrixSize = n;

    % Setting up the root system
    Root_System = RootSystem('B',q,n);

    % Building the form matrix B
    if n > 2*q
        vec_C = sym('c',[1,n-2*q]);
    else
        vec_C = [];
    end
    Form = NIForm(n,q,1,vec_C,0,'symmetric bilinear');

    RootSpaceDimension = @RootSpaceDimensionSO;
    RootSpaceMap = @LieX_SO;
    RootSubgroupMap = @X_SO;
    WeylGroupMap = @W_SO;
    GenericTorusElementMap = @GenericTorusElementSO;
    IsGroupElement = @IsInSO;
    IsTorusElement = @IsTorusElementSO;
    IsLieAlgebraElement = @IsIn_little_so;
    HomDefectCoefficientSO = @HomDefectCoefficientSO;
    CommutatorCoefficientMap = @CommutatorCoefficientSO;
    WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSO;

    SO_n_q = PinnedGroup(NameString,MatrixSize,Root_System,Form,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        HomDefectCoefficientSO, CommutatorCoefficientMap,WeylGroupCoefficientMap);

    SO_n_q.RunTests();

end

function bool = IsIn_little_so(MatrixSize, MatrixToTest, Form)
    % Return true if X belongs to the special orthogonal Lie algebra defined by B
    bool = length(MatrixToTest)==MatrixSize &&...
        SymbolicIsEqual(transpose(MatrixToTest)*Form.Matrix,-Form.Matrix*MatrixToTest);
end
function bool = IsInSO(MatrixSize, MatrixToTest, Form)
    % Return true if X=MatrixToTest belongs to the special orthogonal group 
    % defined by B=Form.Matrix
    bool = length(MatrixToTest)==MatrixSize && ...
        SymbolicIsEqual(transpose(MatrixToTest)*Form.Matrix*MatrixToTest, Form.Matrix);
end
function myMatrix = GenericTorusElementSO(MatrixSize, RootSystemRank, vec_t)
    assert(length(vec_t)==RootSystemRank);
    myMatrix = sym(eye(MatrixSize));
    for i=1:length(vec_t)
        myMatrix(i,i) = vec_t(i);
        myMatrix(RootSystemRank + i, RootSystemRank + i) = vec_t(i)^(-1);
    end
    % for i=1:MatrixSize - 2*RootSystemRank
    %     myMatrix(RootSystemRank+i,RootSystemRank+i) = 1;
    % end
    assert(length(myMatrix)==MatrixSize)
    assert(det(myMatrix)==1)
end
function bool = IsTorusElementSO(MatrixSize, RootSystemRank, MatrixToTest)
    
    n = MatrixSize;
    q = RootSystemRank;
    T = MatrixToTest;

    % Check for correct matrix dimensions
    bool = (length(T)==n) && SymbolicIsDiag(T);

    % Check that the first q entries and the next q diagonal entries are
    % respective inverses
    for i=1:q
        bool = bool && (T(i,i) == T(q+i,q+i)^(-1));
    end

    % Check that the final block is the identity
    for i=2*q+1:n
        bool = bool && (T(i,i) == 1);
    end
        
end
function num = HomDefectCoefficientSO(MatrixSize,RootSystem,Form,u,v)  %#ok<INUSD> 
    num = 0;
end