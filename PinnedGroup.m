classdef PinnedGroup

    properties
        % To do: add type requirements
        NameString
        MatrixSize
        Root_System
        RootList
        RootSystemRank
        FormMatrix
        RootSpaceDimension
        RootSpaceMap
        RootSubgroupMap
        WeylGroupMap
        GenericTorusElementMap
        IsGroupElement
        IsTorusElement
        IsLieAlgebraElement
        CommutatorCoefficientMap
        WeylGroupCoefficientMap
    end

    methods

        % Constructor
        function obj = PinnedGroup(NameString, MatrixSize, ...
            Root_System, FormMatrix, RootSpaceDimension, ...
            RootSpaceMap, RootSubgroupMap, WeylGroupMap, GenericTorusElementMap, ...
            IsGroupElement, IsTorusElement, IsLieAlgebraElement, ...
            CommutatorCoefficientMap, WeylGroupCoefficientMap)
            obj.NameString = NameString;
            obj.MatrixSize = MatrixSize;

            obj.Root_System = Root_System;
            obj.RootList = Root_System.RootList;
            obj.RootSystemRank = Root_System.Rank;
            assert(Root_System.VectorLength == MatrixSize);

            obj.FormMatrix = FormMatrix;
            obj.RootSpaceDimension = RootSpaceDimension;
            obj.RootSpaceMap = RootSpaceMap;
            obj.RootSubgroupMap = RootSubgroupMap;
            obj.WeylGroupMap = WeylGroupMap;
            obj.GenericTorusElementMap = GenericTorusElementMap;
            obj.IsGroupElement = IsGroupElement;
            obj.IsTorusElement = IsTorusElement;
            obj.IsLieAlgebraElement = IsLieAlgebraElement;
            obj.CommutatorCoefficientMap = CommutatorCoefficientMap;
            obj.WeylGroupCoefficientMap = WeylGroupCoefficientMap;
        end

        % Tests
        function RunTests(obj)
            fprintf("Running tests to verify a pinning of the " + obj.NameString + "...\n")
            TestBasics(obj);
            TestRootSubgroupMapsAreHomomorphisms(obj);
            TestTorusConjugationFormula(obj);
            TestCommutatorFormula(obj);
            TestWeylGroupElements(obj);
            TestWeylGroupConjugationFormula(obj);
            fprintf("\n\nAll tests passed.\n\n")
        end
        function TestBasics(obj)
            fprintf("\n\tChecking basic properties...");
            fprintf("\n\t\tChecking root spaces belong to Lie algebra...")
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
                LieX_alpha_u = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                assert(obj.IsLieAlgebraElement(obj.MatrixSize,LieX_alpha_u,obj.FormMatrix))
            end
            fprintf("passed.")
        
            fprintf("\n\t\tChecking root subgroups belong to the group...")
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                assert(obj.IsGroupElement(obj.MatrixSize,X_alpha_u,obj.FormMatrix));
            end
            fprintf("passed.")
        
            fprintf("\n\tBasic tests passed.\n")
        end
        function TestRootSubgroupMapsAreHomomorphisms(obj)
            % Run tests to confirm that RootSubgroupMap is a homomorphism
            % That is, RootSubgroupMap(alpha,u+v) =
            % RootSubgroupMap(alpha,u)*RootSubgroupMap(alpha,v)
            % for symbolic variables u and v
        
            fprintf("\n\tChecking root subgroup maps are homomorphisms...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
                v = sym('v',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                X_alpha_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,v);
                product = X_alpha_u*X_alpha_v;
                X_alpha_u_plus_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u+v);
                assert(SymbolicIsEqual(product,X_alpha_u_plus_v));
            end
            fprintf("passed.")
        end
        function TestTorusConjugationFormula(obj)
        
            fprintf("\n\tChecking torus conjugation formula...");
            vec_t = sym('t',[obj.RootSystemRank,1]);
            t = obj.GenericTorusElementMap(obj.MatrixSize, obj.RootSystemRank, vec_t);
        
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
        
                alpha_of_t = PinnedGroup.CharacterEval(obj.MatrixSize,alpha,t);
        
                LieX_alpha_u = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                LHS1 = t*LieX_alpha_u*t^(-1);
                RHS1 = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,alpha_of_t*u);
                assert(SymbolicIsEqual(LHS1,RHS1));
        
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                LHS2 = t*X_alpha_u*t^(-1);
                RHS2 = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,alpha_of_t*u);
                assert(SymbolicIsEqual(LHS2,RHS2));
            end
            
            fprintf("passed.")
        end
        function TestCommutatorFormula(obj)
            fprintf("\n\tChecking commutator formula...");
        
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
        
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
        
                    if obj.Root_System.IsRoot(alpha+beta) && ~RootSystem.IsProportionalRoot(alpha,beta)
                        dim_V_beta = obj.RootSpaceDimension(obj.Root_System,beta);
                        v = sym('v',[dim_V_beta,1]);
                        X_beta_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,beta,v);
                        LHS = Commutator(X_alpha_u, X_beta_v);
                        
                        RHS = eye(length(LHS));
                        combos = obj.Root_System.LinearCombos(alpha,beta);
                        for k=1:length(combos)
                            root = combos{k}{1};
                            p = combos{k}{2};
                            q = combos{k}{3};
                            assert(isequal(root,p*alpha+q*beta))
                            N = obj.CommutatorCoefficientMap(obj.MatrixSize,obj.Root_System,alpha,beta,p,q,u,v);
                            RHS = RHS * obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,p*alpha+q*beta,N);
                        end
        
                        assert(SymbolicIsEqual(LHS,RHS));
                    end
                end
            end
            fprintf("passed.")
        end
        function TestWeylGroupElements(obj)
            fprintf("\n\tChecking Weyl group elements normalize the torus...")
            vec_t = sym('t',[obj.RootSystemRank,1]);
            t = obj.GenericTorusElementMap(obj.MatrixSize, obj.RootSystemRank, vec_t);
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);
                u = sym('u',[dim_V_alpha,1]);
                w_alpha_u = obj.WeylGroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,u);
                conjugation = simplify(w_alpha_u*t*w_alpha_u^(-1));
                assert(obj.IsTorusElement(obj.MatrixSize,obj.RootSystemRank,conjugation));
            end
            fprintf("passed.")
        end
        function TestWeylGroupConjugationFormula(obj)    
        
            fprintf("\n\tChecking Weyl group conjugation formula...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.Root_System,alpha);

                % OLD VERSION
                %vec1 = ones(1,dim_V_alpha);

                % NEW VERSION
                vec1 = zeros(1,dim_V_alpha);
                vec1(1) = 1;

                w_alpha_1 = simplify(obj.WeylGroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,vec1));
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
                    dim_V_beta = obj.RootSpaceDimension(obj.Root_System,beta);
                    v = sym('v',[dim_V_beta,1]);
        
                    X_beta_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,beta,v);
                    LHS = simplify(w_alpha_1*X_beta_v*w_alpha_1^(-1));
        
                    reflected_root = RootSystem.ReflectRoot(alpha,beta);
                    coeff = obj.WeylGroupCoefficientMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,alpha,beta,v);
                    dim_V_reflected_root = obj.RootSpaceDimension(obj.Root_System,reflected_root);
                    assert(dim_V_reflected_root == dim_V_beta)
                    assert(length(coeff)==dim_V_reflected_root)
                    RHS = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.FormMatrix,reflected_root,coeff);
                    
                    assert(SymbolicIsEqual(LHS,RHS))

                end
            end
            fprintf("passed.");
        end

    end

    methods (Static)

        function val = CharacterEval(MatrixSize, alpha, DiagonalMatrix)
            % Given a positive integer MatrixSize,
            % a character alpha (a vector of length n, which maybe a root but not necessarily)
            % and an (n by n) diagonal matrix t,
            % compute alpha(t)
            assert(SymbolicIsDiag(DiagonalMatrix));
            assert(length(DiagonalMatrix) == MatrixSize);
            assert(length(alpha) == MatrixSize);
        
            % EXAMPLE CALCULATION FOR THE CASE N=3
            %val = t(1,1)^(alpha(1))*t(2,2)^(alpha(2))*t(3,3)^(alpha(3));
            
            val = 1;
            for i = 1:MatrixSize
                val = val*DiagonalMatrix(i,i)^(alpha(i));
            end
        end
    end

end