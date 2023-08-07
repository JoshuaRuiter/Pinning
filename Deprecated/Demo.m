function Demo()
    % Some code for demonstrating various features
    %SpecialLinearDemo(3)
    SpecialUnitaryDemo(5,2,1)
end

function SpecialLinearDemo(n)
    NameString = strcat('special linear group of size',{' '},num2str(n));
    MatrixSize = n;
    RootSystemRank = n-1;
    Root_System = RootSystem('A',RootSystemRank,n);
    RootSpaceDimension = @RootSpaceDimensionSL;
    RootSpaceMap = @LieX_SL;
    RootSubgroupMap = @X_SL;
    WeylGroupMap = @W_SL;
    GenericTorusElementMap = @GenericTorusElementSL;
    I_n = eye(n);
    FormMatrix = I_n; % This is irrelevant for SL_n
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

    SLn
    Phi = SLn.Root_System
    %SLn.Root_System.RootList

    syms u
    syms v
    alpha = SLn.Root_System.RootList{1}
%     beta = SLn.Root_System.RootList{2}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic root spaces demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     root_space = LieX_SL(n,Phi,I_n,alpha,u)
%     root_subgroup = X_SL(n,Phi,I_n,alpha,u)

%     X_alpha_u = X_SL(n,Phi,I_n,alpha,u)
%     X_alpha_v = X_SL(n,Phi,I_n,alpha,v)
%     product = X_alpha_u*X_alpha_v
%     X_alpha_u_plus_v = X_SL(n,Phi,I_n,alpha,u+v)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weyl group conjugation demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     X_alpha_u = X_SL(n,Phi,I_n,alpha,u)
%     w_beta = W_SL(n,Phi,I_n,alpha,1)
%     
%     w_X_w_inv = w_beta*X_alpha_u*w_beta^(-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commutator formula demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = SLn.Root_System.RootList{4}
    sum = alpha + beta
    X_alpha_u = X_SL(n,Phi,I_n,alpha,u)
    X_beta_v = X_SL(n,Phi,I_n,beta,v)
    commutator = Commutator(X_alpha_u,X_beta_v)
    coefficient = CommutatorCoefficientSL(n,Phi,I_n,alpha,beta,1,1,u,v)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Necessary functions to run special linear group demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function SpecialUnitaryDemo(n,q,eps)

    % Use the PinnedGroup class to verify a pinning of SU_{n,q}

    MatrixSize = n;
    RootSystemRank = q;
    syms P;
    PrimitiveElement = P;

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
    H = Form.Matrix;

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
        WeylGroupCoefficientMap)

    Phi = SU_n_q.Root_System
    H

    syms u1
    syms u2
    u = [u1,u2];
    syms v1
    syms v2
    v = [v1,v2];

    alpha = Phi.RootList{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic root spaces demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    root_space = LieX_SU(n,Phi,Form,alpha,u)
    root_subgroup = X_SU(n,Phi,Form,alpha,u)

    X_alpha_u = X_SU(n,Phi,Form,alpha,u)
    X_alpha_v = X_SU(n,Phi,Form,alpha,v)
    product = X_alpha_u*X_alpha_v
    X_alpha_u_plus_v = X_SU(n,Phi,Form,alpha,u+v)
    are_they_equal = isequal(product,X_alpha_u_plus_v)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Necessary functions to run special unitary group demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%