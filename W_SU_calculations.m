function W_SU_calculations()
    n = 10;
    q = 5;
    my_directory = 'C:\Users\Joshua\Documents\Math\Research\Grinnell - Map 499 Spring 2023\Matlab\';
    output_filename = [my_directory,'SU_10_5_weyl'];
    
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
        CommutatorCoefficientMap,WeylGroupCoefficientMap);
    obj = SU_n_q;

    % create a blank table
    % Each row of the table decribes a root alpha, root beta, and the
    %   resulting conjugation coefficient
    % alpha is stored in the first q columns, beta in the q+1 to 2q
    % columns, and the coefficient in the 2q+1 column
    table = cell(length(root_list)^2,2*q+1);
    row_number = 1;

    for i=1:length(root_list)
        alpha = root_list{i};
        dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
        vec1 = zeros(1,dim_V_alpha);
        vec1(1) = 1;
        w_alpha_1 = W_SU(n,root_system,FormMatrix,alpha,vec1);
        w_alpha_1_inverse = w_alpha_1^(-1);

        for j=1:length(root_list)

            beta = root_list{j};
            dim_V_beta = obj.RootSpaceDimension(obj.Root_System,beta);
            v = sym('v',[dim_V_beta,1]);
            X_beta_v = X_SU(n,root_system,FormMatrix,beta,v);
            conjugation = simplify(w_alpha_1*X_beta_v*w_alpha_1_inverse);

            reflected_root = RootSystem.ReflectRoot(alpha,beta);
            assert(RootSpaceDimensionSU(root_system,reflected_root)==dim_V_beta);

            if dim_V_beta == 1
                syms phi_v
                RHS = X_SU(n,root_system,FormMatrix,reflected_root,phi_v);
                position = find(RHS==phi_v);
                coeff = conjugation(position);

            elseif dim_V_beta == 2
                syms phi_v_1
                syms phi_v_2
                phi_v = [phi_v_1,phi_v_2];
                RHS = X_SU(n,root_system,FormMatrix,reflected_root,phi_v);
                position = find(RHS==phi_v_1+phi_v_2*1i);
                coeff = conjugation(position);

            else
                assert(false)
            end

            % the first q entries of the row describe alpha
            % the q+1 to 2q entries of the row describe beta
            % the last (2q+1) entry gives the coefficient
            for k=1:q
                table{row_number,k} = alpha(k);
                table{row_number,q+k} = beta(k);            
            end
            table{row_number,2*q+1} = string(coeff);
            row_number = row_number+1;
        end
    end

    writecell(table,output_filename);

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


