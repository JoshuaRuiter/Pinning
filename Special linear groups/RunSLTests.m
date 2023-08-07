function RunSLTests(n)

    NameString = strcat('special linear group of size',{' '},num2str(n));
    MatrixSize = n;
    RootSystemRank = n-1;
    Root_System = RootSystem('A',RootSystemRank,n);
    RootSpaceDimension = @RootSpaceDimensionSL;
    RootSpaceMap = @LieX_SL;
    RootSubgroupMap = @X_SL;
    WeylGroupMap = @W_SL;
    GenericTorusElementMap = @GenericTorusElementSL;
    FormMatrix = eye(n); % This is irrelevant for SL_n
    IsGroupElement = @IsInSL;
    IsTorusElement = @IsTorusElementSL;
    IsLieAlgebraElement = @IsIn_little_sl;
    HomDefectCoefficientMap = @HomDefectCoefficientSL;
    CommutatorCoefficientMap = @CommutatorCoefficientSL;
    WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSL;

    SLn = PinnedGroup(NameString,MatrixSize,Root_System,FormMatrix,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        HomDefectCoefficientMap, CommutatorCoefficientMap,WeylGroupCoefficientMap);

    SLn.RunTests()
end

function bool = IsTorusElementSL(MatrixSize, RootSystemRank, MatrixToTest)
    % Return true if X belongs to the diagonal torus of SL_n
    bool = (length(MatrixToTest)==MatrixSize && SymbolicIsDiag(MatrixToTest) && det(MatrixToTest)==1);
end
function bool = IsInSL(MatrixSize, MatrixToTest, FormMatrix)
    bool = (length(MatrixToTest)==MatrixSize && det(MatrixToTest)==1);
end
function bool = IsIn_little_sl(MatrixSize, MatrixToTest, FormMatrix)
    bool = (length(MatrixToTest)==MatrixSize && trace(MatrixToTest)==0);
end
function myMatrix = GenericTorusElementSL(MatrixSize, RootSystemRank, vec_t)
    % return a symbolic diagonal matrix with entries vec_T along the
    % diagonal

    % RootSystemRank is n-1, where n is the size of the matrix to return
    assert(MatrixSize == RootSystemRank + 1);
    assert(length(vec_t)==RootSystemRank);
    
    % append an extra entry to the end
    vec_t(MatrixSize) = 1;
    for i=1:(MatrixSize-1)
        vec_t(MatrixSize) = vec_t(MatrixSize)*vec_t(i)^(-1);
    end

    assert(length(vec_t)==MatrixSize);
    myMatrix = diag(vec_t);
    assert(length(myMatrix)==MatrixSize);
end
function dim = RootSpaceDimensionSL(MatrixSize,Root_System,alpha) %#ok<INUSD> 
    dim = 1;
end

function num = HomDefectCoefficientSL(MatrixSize,RootSystem,Form,u,v)  %#ok<INUSD> 
    num = 0;
end