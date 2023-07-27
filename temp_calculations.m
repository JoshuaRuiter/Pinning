% Calculating some example commutator coefficients 
% for special orthogonal groups

n = 7;
q = 2;

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

alpha = [1,0,0,0,0,0,0];
dim_V_alpha = RootSpaceDimensionSO(MatrixSize,Root_System,alpha);
u = sym('u',[dim_V_alpha,1]);
X_alpha_u = X_SO(MatrixSize, Root_System, Form, alpha, u);

beta = [0,1,0,0,0,0,0];
dim_V_beta = RootSpaceDimensionSO(MatrixSize,Root_System,beta);
v = sym('v',[dim_V_beta,1]);
X_beta_v = X_SO(MatrixSize, Root_System, Form, beta, v);

commutator = Commutator(X_alpha_u,X_beta_v)

sum = alpha+beta;
dim_V_sum = RootSpaceDimensionSO(MatrixSize,Root_System,sum);
w = sym('w',[dim_V_sum,1]);
X_sum_w = X_SO(MatrixSize, Root_System, Form, sum, w)

for i=1:length(Root_System.RootList)
    alpha = Root_System.RootList{i};
    dim_V_alpha = RootSpaceDimensionSO(MatrixSize,Root_System,alpha);
    u = sym('u',[dim_V_alpha,1]);
    X_alpha_u = X_SO(MatrixSize, Root_System, Form, alpha, u);
    for j=1:length(Root_System.RootList)
        beta = Root_System.RootList{j};
        dim_V_beta = RootSpaceDimensionSO(MatrixSize,Root_System,beta);
        v = sym('v',[dim_V_beta,1]);
        X_beta_v = X_SO(MatrixSize, Root_System, Form, beta, v);

        if not(RootSystem.IsProportionalRoot(alpha,beta)) ...
                && Root_System.IsRoot(alpha+beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Two Long roots %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if dot(alpha,alpha)==2 && dot(beta,beta)==2
%                 sum = alpha+beta;
%                 syms N
%                 X_sum_N = X_SO(MatrixSize,Root_System,Form,sum,N);
%                 position = find(X_sum_N==N);
%                 commutator = Commutator(X_alpha_u,X_beta_v);
%                 commutator(position)
%             end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Two Short roots %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if dot(alpha,alpha)==1 && dot(beta,beta)==1
%                 sum = alpha+beta;
%                 syms N
%                 X_sum_N = X_SO(MatrixSize,Root_System,Form,sum,N);
%                 position = find(X_sum_N==N);
%                 commutator = Commutator(X_alpha_u,X_beta_v);
%                 commutator(position)
%             end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Two roots of different lengths %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note we just do the case where alpha is short
% and beta is long, since the other case is symmetrical.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if dot(alpha,alpha)==1 && dot(beta,beta)==2
% 
%                 sum11 = alpha+beta;
%                 sum21 = 2*alpha+beta;
%                 assert(Root_System.IsRoot(sum11));
%                 assert(Root_System.IsRoot(sum21));
%                 combos = Root_System.LinearCombos(alpha,beta);
%                 assert(length(combos)==2); % There are exactly two linear combos
%                 assert(dot(sum11,sum11)==1); % alpha+beta is short
%                 assert(dot(sum21,sum21)==2); % 2*alpha+beta is long
%                 assert(RootSpaceDimensionSO(MatrixSize,Root_System,sum11)==n-2*q)
%                 assert(RootSpaceDimensionSO(MatrixSize,Root_System,sum21)==1)
% 
%                 N11 = sym('N11_',[n-2*q,1]); % N11 is a vector of length n-2q
%                 syms N21 % N21 is a scalar
% 
%                 X_sum11_N11 = X_SO(MatrixSize,Root_System,Form,sum11,N11);
%                 X_sum21_N21 = X_SO(MatrixSize,Root_System,Form,sum21,N21);
%                 commutator = Commutator(X_alpha_u,X_beta_v);
% 
%                 pos_N21 = find(X_sum21_N21==N21);
%                 pos_N11_1 = find(X_sum11_N11==N11(1));
%                 pos_N11_2 = find(X_sum11_N11==N11(2));
% 
%                 fprintf("N11 = ")
%                 commutator(pos_N11_1)
%                 commutator(pos_N11_2)
%                 commutator(pos_N21)
%             end

        end

    end
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