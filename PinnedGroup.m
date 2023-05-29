classdef PinnedGroup

    properties
        % To do: add type requirements
        NameString
        MatrixSize
        Root_System
        RootList
        RootSystemRank
        Form
        RootSpaceDimension
        RootSpaceMap
        RootSubgroupMap
        WeylGroupMap
        GenericTorusElementMap
        IsGroupElement
        IsTorusElement
        IsLieAlgebraElement
        HomDefectCoefficientMap
        CommutatorCoefficientMap
        WeylGroupCoefficientMap
    end

    methods

        % Constructor
        function obj = PinnedGroup(NameString, MatrixSize, ...
            Root_System, Form, RootSpaceDimension, ...
            RootSpaceMap, RootSubgroupMap, WeylGroupMap, GenericTorusElementMap, ...
            IsGroupElement, IsTorusElement, IsLieAlgebraElement, ...
            HomDefectCoefficientMap, CommutatorCoefficientMap, WeylGroupCoefficientMap)
            obj.NameString = NameString;
            obj.MatrixSize = MatrixSize;

            obj.Root_System = Root_System;
            obj.RootList = Root_System.RootList;
            obj.RootSystemRank = Root_System.Rank;
            assert(Root_System.VectorLength == MatrixSize);

            obj.Form = Form;
            obj.RootSpaceDimension = RootSpaceDimension;
            obj.RootSpaceMap = RootSpaceMap;
            obj.RootSubgroupMap = RootSubgroupMap;
            obj.WeylGroupMap = WeylGroupMap;
            obj.GenericTorusElementMap = GenericTorusElementMap;
            obj.IsGroupElement = IsGroupElement;
            obj.IsTorusElement = IsTorusElement;
            obj.IsLieAlgebraElement = IsLieAlgebraElement;
            obj.HomDefectCoefficientMap = HomDefectCoefficientMap;
            obj.CommutatorCoefficientMap = CommutatorCoefficientMap;
            obj.WeylGroupCoefficientMap = WeylGroupCoefficientMap;
        end

        % Tests
        function RunTests(obj)
%             obj.Root_System.VerifyProperties();
            fprintf("Running tests to verify a pinning of the " + obj.NameString + "...\n")
%             TestBasics(obj);
%             TestRootSpaceMapsAreHomomorphisms(obj);
%             TestRootSubgroupMapsAreAlmostHomomorphisms(obj);
%             TestTorusConjugationFormula(obj);
            TestCommutatorFormula(obj);
%             TestWeylGroupElements(obj);
%             TestWeylGroupConjugationFormula(obj);
            fprintf("\n\nAll tests passed.\n\n")
        end
        function TestBasics(obj)
            fprintf("\n\tChecking basic properties...");

            fprintf("\n\t\tChecking root spaces belong to Lie algebra...")
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                LieX_alpha_u = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                assert(obj.IsLieAlgebraElement(obj.MatrixSize,LieX_alpha_u,obj.Form))
            end
            fprintf("passed.")
        
            fprintf("\n\t\tChecking root subgroups belong to the group...")
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                assert(obj.IsGroupElement(obj.MatrixSize,X_alpha_u,obj.Form));
            end
            fprintf("passed.")
        
            fprintf("\n\tBasic tests passed.\n")
        end
        function TestRootSpaceMapsAreHomomorphisms(obj)
            % Run tests to confirm that RootSpaceMap is a homomorphism
            % That is, RootSpaceMap(alpha,u+v) =
            % RootSpaceMap(alpha,u)+RootSubgroupMap(alpha,v)
            % for symbolic variables u and v

            fprintf("\n\tChecking root space maps are homomorphisms...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                v = sym('v',[dim_V_alpha,1]);
                L_alpha_u = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                L_alpha_v = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,v);
                product = L_alpha_u+L_alpha_v;
                L_alpha_u_plus_v = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u+v);
                assert(SymbolicIsEqual(product,L_alpha_u_plus_v));
            end
            fprintf("passed.")
        end
        function TestRootSubgroupMapsAreAlmostHomomorphisms(obj)
            % Run tests to confirm that RootSubgroupMap is almost a homomorphism
            % That is, RootSubgroupMap(alpha,u)*RootSubgroupMap(alpha,v) =
            % RootSubgroupMap(alpha,u+v)*RootSubgroupMap(2alpha,q_alpha(u,v))
            % in theory this can continue if higher multiples of alpha are roots, 
            % but that can't ever happen because of the structure of root systems.

            % The function q_alpha is called the "Hom(omomorphism)DefectCoefficient map"
            % because it captures exactly how much RootSubgroupMap(alpha,-) 
            % fails to be a homomorphism.

            % Often, such as in the special linear group, q_alpha(u,v)=0 and so these
            % are bonafide homomorphisms.

            fprintf("\n\tChecking root subgroup maps are almost homomorphisms...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                v = sym('v',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,...
                    obj.Root_System,obj.Form,alpha,u);
                X_alpha_v = obj.RootSubgroupMap(obj.MatrixSize,...
                    obj.Root_System,obj.Form,alpha,v);
                product = X_alpha_u*X_alpha_v;

                X_alpha_u_plus_v = obj.RootSubgroupMap(obj.MatrixSize,...
                    obj.Root_System,obj.Form,alpha,u+v);
                RHS = X_alpha_u_plus_v;
                
                if obj.Root_System.IsRoot(2*alpha)
                    q_uv = obj.HomDefectCoefficientMap(obj.MatrixSize,...
                        obj.Root_System,obj.Form,alpha,u,v);
                    X_2alpha_quv = obj.RootSubgroupMap(obj.MatrixSize,...
                        obj.Root_System,obj.Form,2*alpha,q_uv);
                    RHS = RHS * X_2alpha_quv;
                end
                assert(SymbolicIsEqual(product,RHS));
            end
            fprintf("passed.")
        end
        function TestTorusConjugationFormula(obj)
        
            fprintf("\n\tChecking torus conjugation formula...");
            vec_t = sym('t',[obj.RootSystemRank,1]);
            t = obj.GenericTorusElementMap(obj.MatrixSize, obj.RootSystemRank, vec_t);
        
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
        
                alpha_of_t = PinnedGroup.CharacterEval(obj.MatrixSize,alpha,t);
        
                LieX_alpha_u = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                LHS1 = t*LieX_alpha_u*t^(-1);
                RHS1 = obj.RootSpaceMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,alpha_of_t*u);
                assert(SymbolicIsEqual(LHS1,RHS1));
        
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                LHS2 = t*X_alpha_u*t^(-1);
                RHS2 = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,alpha_of_t*u);
                assert(SymbolicIsEqual(LHS2,RHS2));
            end
            
            fprintf("passed.")
        end
        function TestCommutatorFormula(obj)

            warning('off','all') % TEMPORARY

            fprintf("\n\tChecking commutator formula...");

            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                X_alpha_u = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
        
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
        
                    if obj.Root_System.IsRoot(alpha+beta) && ~RootSystem.IsProportionalRoot(alpha,beta)
                        dim_V_beta = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, beta);
                        v = sym('v',[dim_V_beta,1]);
                        X_beta_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,beta,v);
                        LHS = Commutator(X_alpha_u, X_beta_v);

                        RHS = eye(length(LHS));
                        combos = obj.Root_System.LinearCombos(alpha,beta);
                        for k=1:length(combos)
                            root = combos{k}{1};
                            p = combos{k}{2};
                            q = combos{k}{3};
                            assert(isequal(root,p*alpha+q*beta))
                            N = obj.CommutatorCoefficientMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,beta,p,q,u,v);
                            RHS = RHS * obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,p*alpha+q*beta,N);
                        end

                        % ORIGINAL TEST
%                         assert(SymbolicIsEqual(LHS,RHS))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin segmented tests for special unitary group commutator coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         % alpha medium, beta long
%                         if dot(alpha,alpha) == 2 && dot(beta,beta) == 4
%                             % Passes for all special unitary groups
%                             assert(SymbolicIsEqual(LHS,RHS));
%                         end
% 
%                         % alpha long, beta medium
%                         if dot(alpha,alpha) == 4 && dot(beta,beta) == 2
%                             % Passes for all special unitary groups
%                             assert(SymbolicIsEqual(LHS,RHS));
%                         end
% 
%                         % alpha medium, beta medium, alpha+beta long
%                         if dot(alpha,alpha)==2 && dot(beta,beta)==2 && dot(alpha+beta,alpha+beta)==4
%                             % Passes for all special unitary groups
%                             assert(SymbolicIsEqual(LHS,RHS));
%                         end
% 
%                         % alpha medium, beta medium, alpha+beta medium
%                         if dot(alpha,alpha)==2 && dot(beta,beta)==2 && dot(alpha+beta,alpha+beta)==2
%                             % These tests pass for all special unitary groups
%                                
%                             mn = find(alpha);
%                             assert(length(mn)==2)
%                             m = mn(1);
%                             n = mn(2);
%                             assert(m<n)
%                             pq = find(beta);
%                             assert(length(pq)==2)
%                             p = pq(1);
%                             q = pq(2);
%                             assert(p<q);
% 
%                             alpha_m = zeros(1,obj.Root_System.VectorLength);
%                             alpha_m(m) = 1;
%                             alpha_n = zeros(1,obj.Root_System.VectorLength);
%                             alpha_n(n) = 1;
%                             alpha_p = zeros(1,obj.Root_System.VectorLength);
%                             alpha_p(p) = 1;
%                             alpha_q = zeros(1,obj.Root_System.VectorLength);
%                             alpha_q(q) = 1;
%                             eps_m = alpha(m);
%                             eps_n = alpha(n);
%                             eps_p = beta(p);
%                             eps_q = beta(q);
%                             assert(isequal(alpha,eps_m*alpha_m+eps_n*alpha_n));
%                             assert(isequal(beta,eps_p*alpha_p+eps_q*alpha_q));
% 
%                             if m == p || m == q
%                                 r = m;
%                                 s = n;
%                                 if m == p
%                                     t = q;
%                                 else
%                                     t = p;
%                                 end
%                             elseif n == p || n == q
%                                 r = n;
%                                 s = m;
%                                 if n == p
%                                     t = q;
%                                 else
%                                     t = p;
%                                 end
%                             end
% 
%                             alpha_r = zeros(1,obj.Root_System.VectorLength);
%                             alpha_r(r) = 1;
%                             alpha_s = zeros(1,obj.Root_System.VectorLength);
%                             alpha_s(s) = 1;
%                             alpha_t = zeros(1,obj.Root_System.VectorLength);
%                             alpha_t(t) = 1;
%                             eps_r = alpha(r);
%                             eps_s = alpha(s);
%                             eps_t = beta(t);
%                             assert(beta(r)==-eps_r)
%                             assert(isequal(alpha,eps_r*alpha_r+eps_s*alpha_s))
%                             assert(isequal(beta,-eps_r*alpha_r+eps_t*alpha_t))
% 
%                             if alpha(r)==-alpha(s) && beta(r)==-beta(t)
%                                 % The entries of each root have opposite signs
%                                 % This behaves essentially like a copy of SL_n
%                                 % Tests pass for (6,3,1) and (6,3,-1)
%                                 fprintf("A")
% %                                 alpha
% %                                 beta
% %                                 LHS
% %                                 simplify(rdivide(LHS,RHS))
%                                 assert(SymbolicIsEqual(LHS,RHS));
%                                 
%                             elseif r < s && r < t
%                                 % r is smaller than both s and t
%                                 if alpha(r)==alpha(s) && beta(r)==beta(t)
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("B1")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 elseif (t-s)*beta(r)*beta(t) > 0
%                                     % Entries with the same sign are spaced farther apart than
%                                     % entries with the opposite sign.
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("B2")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 else
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("B3")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 end
%                         
%                             elseif r > s && r > t
%                                 % r is bigger than both s and t
%                                 if alpha(r)==alpha(s) && beta(r)==beta(t)
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("C1")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 elseif s < t
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("C2")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 else
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("C3")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 end
%                                 
%                             else
%                                 % r is between s and t
%                                 if (t-s)*beta(r)*beta(t)>0 && alpha(s)==beta(t)
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("D1")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 elseif s < t
%                                     % Tests pass for (6,3,1) and (6,3,-1)
%                                     fprintf("D2")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
% 
%                                 else
%                                     fprintf("D3")
% %                                     alpha
% %                                     beta
% %                                     LHS
% %                                     simplify(rdivide(LHS,RHS))
%                                     assert(SymbolicIsEqual(LHS,RHS));
%                                 end
%                             end
%                         end

                        % alpha and beta are short
                        % This only occurs for non-quasisplit groups
                        if dot(alpha,alpha)==1 && dot(beta,beta)==1
%                             % Passes for all special unitary groups
%                             alpha
%                             beta
%                             simplify(LHS)
%                             simplify(rdivide(LHS,RHS))
                            assert(SymbolicIsEqual(LHS,RHS));
                        end

                        % alpha is short, beta is medium
                        % This only occurs for non-quasisplit groups
                        if dot(alpha,alpha)==1 && dot(beta,beta)==2
%                             alpha
%                             beta
%                             LHS
%                             simplify(rdivide(LHS,RHS))
%                             assert(SymbolicIsEqual(LHS,RHS));
                        end

                        % alpha is medium, beta is short
                        % This only occurs for non-quasisplit groups
                        if dot(alpha,alpha)==2 && dot(beta,beta)==1
%                             alpha
%                             beta
%                             LHS
%                             simplify(rdivide(LHS,RHS))
%                             assert(SymbolicIsEqual(LHS,RHS));
                        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End segmented tests for special unitary group commutator coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    end
                end
            end
            fprintf("passed.")

            warning('on','all') % TEMPORARY

        end
        function TestWeylGroupElements(obj)
            fprintf("\n\tChecking Weyl group elements normalize the torus...")
            vec_t = sym('t',[obj.RootSystemRank,1]);
            t = obj.GenericTorusElementMap(obj.MatrixSize, obj.RootSystemRank, vec_t);
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                u = sym('u',[dim_V_alpha,1]);
                w_alpha_u = obj.WeylGroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,u);
                conjugation = simplify(w_alpha_u*t*w_alpha_u^(-1));
                assert(obj.IsTorusElement(obj.MatrixSize,obj.RootSystemRank,conjugation));
            end
            fprintf("passed.")
        end
        function TestWeylGroupConjugationFormula(obj)    
        
            fprintf("\n\tChecking Weyl group conjugation formula...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                dim_V_alpha = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, alpha);
                vec1 = zeros(1,dim_V_alpha);
                vec1(1) = 1;
                w_alpha_1 = simplify(obj.WeylGroupMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,vec1));

                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
                    dim_V_beta = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, beta);
                    v = sym('v',[dim_V_beta,1]);
        
                    X_beta_v = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,beta,v);
                    LHS = simplify(w_alpha_1*X_beta_v*w_alpha_1^(-1));
        
                    reflected_root = RootSystem.ReflectRoot(alpha,beta);
                    coeff = obj.WeylGroupCoefficientMap(obj.MatrixSize,obj.Root_System,obj.Form,alpha,beta,v);
                    dim_V_reflected_root = obj.RootSpaceDimension(obj.MatrixSize, obj.Root_System, reflected_root);
                    assert(dim_V_reflected_root == dim_V_beta)
                    assert(length(coeff)==dim_V_reflected_root)
                    RHS = obj.RootSubgroupMap(obj.MatrixSize,obj.Root_System,obj.Form,reflected_root,coeff);
                    
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