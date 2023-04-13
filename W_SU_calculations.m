function W_SU_calculations()
    % SU_{4,2}
    n = 4;
    q = 2;
    
    NameString = strcat('special unitary group of size',{' '},num2str(n),{' '},'with Witt index',{' '},num2str(q));
    MatrixSize = n;
    RootSystemRank = q;
    
    % Setting up the root system
    if n==2*q
        % Use the C_q root system
        Type = 'C';
    else
        % Use the BC_q root system
        Type = 'BC';
    end
    root_system = RootSystem(Type,RootSystemRank,MatrixSize);
    root_list = root_system.RootList;

    LongRoots = {};
    MediumRoots = {};
    index_long = 1;
    index_medium = 1;
    for i=1:length(root_list)
        alpha = root_list{i};
        if dot(alpha,alpha)==4
            LongRoots{index_long} = alpha;
            index_long = index_long + 1;
        elseif dot(alpha,alpha)==2
            MediumRoots{index_medium} = alpha;
            index_medium = index_medium + 1;
        end
    end

    % Building the FormMatrix
    if n-2*q > 0
        vec_C = sym('c',[1,n-2*q]);
        FormMatrix = GetH(n,q,vec_C);
    else
        FormMatrix = GetH(n,q,[]);
    end
    
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
    
    SU_n_q = PinnedGroup(NameString,MatrixSize,root_system,FormMatrix,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        CommutatorCoefficientMap,WeylGroupCoefficientMap)
    LongRoots
    MediumRoots

    % Two long roots case
    % This case is complete
%     for i=1:length(LongRoots)
%         alpha = LongRoots{i};
%         assert(RootSpaceDimensionSU(root_system,alpha)==1);
%         w_alpha_1 = W_SU(n,root_system,FormMatrix,alpha,1);
%         w_alpha_1_inverse = w_alpha_1^(-1);
%         for j=1:length(LongRoots)
%             beta = LongRoots{j};
%             assert(RootSpaceDimensionSU(root_system,beta)==1);
%             syms v;
%             X_beta_v = X_SU(n,root_system,FormMatrix,beta,v);
%             conjugation = w_alpha_1*X_beta_v*w_alpha_1_inverse;
% 
%             reflected_root = RootSystem.ReflectRoot(alpha,beta);
%             assert(RootSpaceDimensionSU(root_system,beta)==RootSpaceDimensionSU(root_system,reflected_root));
%             syms phi_v;
%             RHS = X_SU(n,root_system,FormMatrix,reflected_root,phi_v);
% 
%             position = find(RHS==phi_v);
% 
%             alpha
%             beta
%             conjugation(position)
%         end
%     end

    % Alpha is medium, beta is long
    % This case is complete
%     for i=1:length(MediumRoots)
%         alpha = MediumRoots{i};
%         assert(RootSpaceDimensionSU(root_system,alpha)==2);
%         w_alpha_1 = W_SU(n,root_system,FormMatrix,alpha,[1,0]);
%         w_alpha_1_inverse = w_alpha_1^(-1);
%         for j=1:length(LongRoots)
%             beta = LongRoots{j};
%             assert(RootSpaceDimensionSU(root_system,beta)==1);
%             syms v;
%             X_beta_v = X_SU(n,root_system,FormMatrix,beta,v);
%             conjugation = w_alpha_1*X_beta_v*w_alpha_1_inverse;
% 
%             reflected_root = RootSystem.ReflectRoot(alpha,beta);
%             assert(RootSpaceDimensionSU(root_system,beta)==RootSpaceDimensionSU(root_system,reflected_root));
%             syms phi_v;
%             RHS = X_SU(n,root_system,FormMatrix,reflected_root,phi_v);
% 
%             position = find(RHS==phi_v);
% 
%             alpha
%             beta
%             conjugation(position)
%         end
%     end

    % The case where alpha is long and beta is medium is also complete
   
    % Both alpha and beta are medium
    for i=1:length(MediumRoots)
        alpha = MediumRoots{i};
        assert(RootSpaceDimensionSU(root_system,alpha)==2);
        w_alpha_1 = W_SU(n,root_system,FormMatrix,alpha,[1,0]);
        w_alpha_1_inverse = w_alpha_1^(-1);
        for j=1:length(MediumRoots)
            beta = MediumRoots{j};
            assert(RootSpaceDimensionSU(root_system,beta)==2);
            syms v1;
            syms v2;
            v = [v1,v2];
            X_beta_v = X_SU(n,root_system,FormMatrix,beta,v);
            conjugation = w_alpha_1*X_beta_v*w_alpha_1_inverse;

            reflected_root = RootSystem.ReflectRoot(alpha,beta);
            assert(RootSpaceDimensionSU(root_system,beta)==RootSpaceDimensionSU(root_system,reflected_root));
            syms phi_v_1;
            syms phi_v_2;
            phi_v = [phi_v_1,phi_v_2];
            RHS = X_SU(n,root_system,FormMatrix,reflected_root,phi_v);

            position = find(RHS==phi_v_1+phi_v_2*1i);
            
            alpha
            beta
            conjugation(position)
        end
    end

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
function bool = IsIn_little_su(MatrixSize,MatrixToTest,FormMatrix)
    bool = (length(MatrixToTest)==MatrixSize && ...
        trace(MatrixToTest)==0 && ...
        SymbolicIsEqual(ctranspose(MatrixToTest)*FormMatrix,-FormMatrix*MatrixToTest));
end
function bool = IsInSU(MatrixSize,MatrixToTest,FormMatrix)
    bool = (length(MatrixToTest)==MatrixSize && ...
        det(MatrixToTest) == 1 && ...
        SymbolicIsEqual(ctranspose(MatrixToTest)*FormMatrix*MatrixToTest,FormMatrix));
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

    mat = SymbolicZeros(MatrixSize);
    for i=1:RootSystemRank
        mat(i,i) = DiagonalValues(i);
        mat(RootSystemRank+i,RootSystemRank+i) = DiagonalValues(i)^(-1);
    end
    for i=2*RootSystemRank+1:MatrixSize
        mat(i,i) = 1;
    end
end


