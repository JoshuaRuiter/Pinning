classdef SCM

    % Symbolic Coordinates Matrix
    % A class to represent a matrix with various entries
    % in various symbolically specified entries
    
    % This is accomplished by simply having a list of SymbolicSingleEntryMatrix objects

    properties
        Size            % An integer, the size of the matrix
        Concrete_base   % A numeric matrix of size obj.Size
        List_of_terms   % A cell array of SymbolicSingleEntryMatrix objects
    end

    methods

        % Constructor
        function obj = SCM(Size,Concrete_base,List_of_terms)
            obj.Size = Size;

            % Default concrete base is the zero matrix
            if nargin == 1
                obj.Concrete_base = zeros(obj.Size);
            else
                assert(length(Concrete_base)==Size)
                obj.Concrete_base = Concrete_base;
            end

            % Default list of terms is empty
            if nargin <= 2
                obj.List_of_terms = {};
            else
                obj.List_of_terms = List_of_terms;
            end
        end

        function bool = isElementary(obj)
            bool = isequal(obj.Concrete_base,eye(obj.Size)) && length(obj.List_of_terms)==1;
        end

        % Overloading unary negative operator
        function negative = uminus(obj)
            num_terms = length(obj.List_of_terms);
            List_of_negative_terms = cell(1,num_terms);
            for t=1:num_terms
                List_of_negative_terms{t} = -ob.List_of_terms{t};
            end
            negative = SCM(obj.Size,-obj.Concrete_base,List_of_negative_terms);
        end

        function inverse = invert(obj)
            % Compute the inverse in elementary matrix case, i.e. when
            % the concrete_base is the identity, and there is only one term
            assert(obj.isElementary())
            inverse = SCM(obj.Size,obj.Concrete_base,{-obj.List_of_terms{1}});
        end

        function output_table = tableForm(obj, with_size)
            % Create a tabular output version of the SCM object
            % The table is formatted as follows:
            %   There is one row for each term in the List_of_terms 
            %   In that row is (Row, Column, Entry)
            %   or (Row, Column, Entry, Size) if with_size == 1
            % with the entry converted to a symbolic format
            % if it is an SKD object

            if nargin == 1
                with_size = 0;
            end

            if with_size
                num_columns = 4;
            else
                num_columns = 3;
            end
            
            num_terms = length(obj.List_of_terms);
            output_table = sym(zeros(num_terms, num_columns));    
            for t=1:length(obj.List_of_terms)
                term_t = obj.List_of_terms{t};
                output_table(t,1:num_columns) = term_t.rowForm();
            end
        end

        function simplified_copy = cancelPairs(obj)
            % If there are any SSEM's in the List_of_terms
            % that are negatives of each other, cancel them
            % out in pairs.

            simplified_copy = obj;
            i=1;
            j=2;

            while i<=length(simplified_copy.List_of_terms)-1
                assert(i<j)
                term_i = simplified_copy.List_of_terms{i};
                while j<=length(simplified_copy.List_of_terms)
                    term_j = simplified_copy.List_of_terms{j};
                    if isequal(term_i,-term_j)
                        % Remove both terms, starting with the later one
                        simplified_copy.List_of_terms(j) = [];
                        simplified_copy.List_of_terms(i) = [];
                       
                        % Adjust i to compensate for now missing entry
                        % and break out of the j loop
                        % Note this only breaks out of the j loop,
                        % not the i loop
                        i=i-1;
                        break
                    end
                    j=j+1;
                end
                i=i+1;
                j=i+1;
            end
        end

        function simplified_copy = simplifyWithAssumption(obj, assumption)
            % Given an assumption, zero out any terms which are
            % zero under that assumption
            simplified_copy = obj;
            t = 1;
            while t <= length(simplified_copy.List_of_terms)
                term_t = simplified_copy.List_of_terms{t}; % an SSEM object
                entry_t = term_t.Entry;                    % could be an SKD object
                if isequal(entry_t,0)
                    % Zero entry, can remove term_t
                    simplified_copy.List_of_terms(t) = [];
                    % Adjust t to compensate for the removed term
                    t = t-1;
                elseif isa(entry_t, 'SKD')
                    conditions_t = entry_t.Conditions; % a list of equations
                    % Iterate over all conditions, check if any single 
                    % condition contradicts the assumption
                    for c=1:length(conditions_t)
                        if isequal(conditions_t(c), ~assumption)
                            % Remove the entry
                            simplified_copy.List_of_terms(t) = [];
                            t = t-1;
                            break; % This only exits the for loop, not the while loop
                        end
                    end
                end
                t = t+1;
            end
        end

        function new_obj = addDerivedConditions(obj)
            new_terms = cell(size(obj.List_of_terms));
            for t=1:length(obj.List_of_terms)
                term_t = obj.List_of_terms{t};
                new_term_t = term_t.addDerivedConditions();
                new_terms{t} = new_term_t;
            end
            new_obj = SCM(obj.Size,obj.Concrete_base,new_terms);
        end

    end

    methods (Static)

        function sum = add(mat1, mat2)
            assert(mat1.Size == mat2.Size)
            sum = SCM(mat1.Size,...
                mat1.Concrete_base+mat2.Concrete_base,...
                [mat1.List_of_terms,mat2.List_of_terms]);
        end

        function product = multiply(obj1, obj2)
            assert(obj1.Size == obj2.Size)
            size = obj1.Size;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WARNING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Currently, this is only implemented in the cases where
            % the concrete bases for each matrix are the zero or identity matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            assert(isequal(obj1.Concrete_base,zeros(size)) ...
                || isequal(obj1.Concrete_base,eye(size)))
            assert(isequal(obj2.Concrete_base,zeros(size)) ...
                || isequal(obj2.Concrete_base,eye(size)))

            length1 = length(obj1.List_of_terms);
            length2 = length(obj2.List_of_terms);

            % Build the list of SSEM objects
            % by doing a big "distribute" on the two lists
            product_list_of_terms = cell(1,length1*length2);
            index = 1;
            for s=1:length1
                ssem1 = obj1.List_of_terms{s};
                for t=1:length2
                    ssem2 = obj2.List_of_terms{t};

                    % Compute the product of two SSEM objects
                    prod = SSEM.multiply(ssem1,ssem2);
                    product_list_of_terms{index} = prod;
                    index = index + 1;
                end
            end

            % If the concrete base for mat1 is the identity, just concatenate on the 
            % list of terms for obj2
            if isequal(obj1.Concrete_base,eye(size))
                product_list_of_terms = [product_list_of_terms,obj2.List_of_terms];
            end
            
            % If the concrete base for mat2 is the identity, just concatenate on the 
            % list of terms for obj1
            if isequal(obj2.Concrete_base,eye(size))
                product_list_of_terms = [product_list_of_terms,obj1.List_of_terms];
            end

            if isequal(obj2.Concrete_base,eye(size)) && isequal(obj1.Concrete_base,eye(size))
                product_concrete_base = eye(size);
            else
                product_concrete_base = zeros(size);
            end

            product = SCM(size,product_concrete_base,product_list_of_terms);
        
        end

        function output = commutator(obj1, obj2)
            % Assuming obj1 and obj1 are elementary, 
            % compute their commutator
            assert(obj1.isElementary())
            assert(obj2.isElementary())

            prod1 = SCM.multiply(obj1,obj2);
            prod2 = SCM.multiply(invert(obj1),invert(obj2));

            output = SCM.multiply(prod1,prod2);

        end

        function runTests()
            fprintf("Running SCM tests...");
                
            syms a
            syms b
            syms c
            syms d
            syms u
            syms v
            n = 5;

            E_ab_u = SSEM(a,b,n,u);
            e_ab_u = SCM(n,eye(n),{E_ab_u});
            assert(e_ab_u.isElementary())
            e_ab_u_inverse = e_ab_u.invert();
            e_ab_u_manual_inverse = SCM(n,eye(n),{-E_ab_u});
            assert(isequal(e_ab_u_inverse,e_ab_u_manual_inverse))

            E_cd_v = SSEM(c,d,n,v);
            e_cd_v = SCM(n,eye(n),{E_cd_v});
            assert(e_cd_v.isElementary())

            E_ab_u_SCM = E_ab_u.convertToSCM();
            E_cd_v_SCM = E_cd_v.convertToSCM();
            assert(isequal(E_ab_u_SCM.Concrete_base,zeros(n)));
            assert(isequal(E_cd_v_SCM.Concrete_base,zeros(n)));
            assert(length(E_ab_u_SCM.List_of_terms)==1);
            assert(length(E_cd_v_SCM.List_of_terms)==1);
            sum = SCM.add(E_ab_u_SCM,E_cd_v_SCM);
            assert(sum.Size == n)
            assert(isequal(sum.Concrete_base,zeros(n)))
            assert(length(sum.List_of_terms)==2)
            
            prod_1 = SCM.multiply(e_ab_u,e_cd_v);
            assert(isequal(prod_1.Concrete_base,eye(n)))
            assert(length(prod_1.List_of_terms)==3)

            I_SCM_v1 = SCM(n,eye(n),{});
            prod_2 = SCM.multiply(e_ab_u,I_SCM_v1);
            assert(isequal(prod_2.Concrete_base,eye(n)))
            assert(length(prod_2.List_of_terms)==1)
            assert(isequal(prod_2,e_ab_u));

            com = SCM.commutator(e_ab_u,e_cd_v);
            assert(isequal(com.Concrete_base,eye(n)))
            assert(length(com.List_of_terms)==15)

            % Add tests for cancelPairs method
            com_after_cancel = com.cancelPairs();

            fprintf(" all tests passed.\n")
        end

    end

end