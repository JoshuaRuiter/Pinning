% Output all root space maps for SO(7,2)

n=7;
q=2;

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
HomDefectCoefficientSO_label = @HomDefectCoefficientSO;
CommutatorCoefficientMap = @CommutatorCoefficientSO;
WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSO;

SO_n_q = PinnedGroup(NameString,MatrixSize,Root_System,Form,...
    RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
    IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
    HomDefectCoefficientSO_label, CommutatorCoefficientMap,WeylGroupCoefficientMap);

for i=1:length(Root_System.RootList)
    alpha = Root_System.RootList{i};
    dim_V_alpha = RootSpaceDimensionSO(n,Root_System,alpha);
    v = sym('v',[dim_V_alpha,1]);

    alpha
    LieX_SO(n,Root_System,Form,alpha,v);
    X_SO(n,Root_System,Form,alpha,v)
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
    
    % Check for correct matrix dimensions
    bool = (length(MatrixToTest)==MatrixSize) && SymbolicIsDiag(MatrixToTest);

    % Check that the first q entries and the last q diagonal entries are
    % respective inverses

    for i=1:RootSystemRank
        bool = bool && (MatrixToTest(i,i) == MatrixToTest(MatrixSize+1-i,MatrixSize+1-i)^(-1));
    end

    % Check that the middle block is the identity
    for i=RootSystemRank+1:MatrixSize-RootSystemRank
        bool = bool && (MatrixToTest(i,i) == 1);
    end
        
end
function num = HomDefectCoefficientSO(MatrixSize,RootSystem,Form,u,v)  %#ok<INUSD> 
    num = 0;
end