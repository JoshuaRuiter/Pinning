function W_SO_calculations()
    % SU_{4,1}
    n = 4;
    q = 1;
    
   NameString = strcat('special orthogonal group of size',{' '},num2str(n),{' '},'with Witt index',{' '},num2str(q));
    MatrixSize = n;

    % Setting up the root system
    root_system = RootSystem('B',q,n);
    root_list = root_system.RootList;

    % Building the FormMatrix
    if n-2*q > 0
        vec_C = sym('c',[1,n-2*q]);
        FormMatrix = GetB(n,q,vec_C);
    else
        FormMatrix = GetB(n,q,[]);
    end

    RootSpaceDimension = @RootSpaceDimensionSO;
    RootSpaceMap = @LieX_SO;
    RootSubgroupMap = @X_SO;
    WeylGroupMap = @W_SO;
    GenericTorusElementMap = @GenericTorusElementSO;
    IsGroupElement = @IsInSO;
    IsTorusElement = @IsTorusElementSO;
    IsLieAlgebraElement = @IsIn_little_so;
    CommutatorCoefficientMap = @CommutatorCoefficientSO;
    WeylGroupCoefficientMap = @WeylGroupConjugationCoefficientSO;

    SO_n_q = PinnedGroup(NameString,MatrixSize,root_system,FormMatrix,...
        RootSpaceDimension,RootSpaceMap,RootSubgroupMap,WeylGroupMap,GenericTorusElementMap,...
        IsGroupElement,IsTorusElement,IsLieAlgebraElement,...
        CommutatorCoefficientMap,WeylGroupCoefficientMap)


    % SO_{4,1} calculations
    % The associated tests for these calculations now pass
%     for i=1:length(root_list)
%         alpha = root_list{i};
%         assert(RootSpaceDimensionSO(root_system,alpha)==2);
%         w_alpha_1 = W_SO(n,root_system,FormMatrix,alpha,[1,0]);
%         w_alpha_1_inverse = w_alpha_1^(-1);
%         for j=1:length(root_list)
%             beta = root_list{j};
%             assert(RootSpaceDimensionSO(root_system,beta)==2);
%             syms v1;
%             syms v2;
%             v = [v1, v2];
%             X_beta_v = X_SO(n,root_system,FormMatrix,beta,v);
%             conjugation = w_alpha_1*X_beta_v*w_alpha_1_inverse;
% 
%             reflected_root = RootSystem.ReflectRoot(alpha,beta);
%             assert(RootSpaceDimensionSO(root_system,beta)==RootSpaceDimensionSO(root_system,reflected_root));
%             syms phi_v_1;
%             syms phi_v_2;
%             phi_v = [phi_v_1,phi_v_2];
%             RHS = X_SO(n,root_system,FormMatrix,reflected_root,phi_v);
% 
%             position_1 = find(RHS==phi_v_1);
%             position_2 = find(RHS==phi_v_2);
% 
%             alpha
%             beta
%             phi_v_1 = conjugation(position_1)
%             phi_v_2 = conjugation(position_2)
%         end
%     end

end

function dim = RootSpaceDimensionSO(Root_System,alpha) %#ok<INUSD>
    dim = Root_System.VectorLength - 2*Root_System.Rank;
end