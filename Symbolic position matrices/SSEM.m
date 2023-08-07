classdef SSEM

    % Symbolic Single Entry Matrix
    % A class to represent a square matrix with a single nonzero entry
    % The row, column, and entry properties can be symbolic variables

    % For example, we can model the symbolic n-by-n  matrix E_ij(u) using 
    % SymbolicSingleEntryMatrix(i,j,n,u)

    % Two such matrices can be multiplied together
    % even in this symbolic state using the multiply() method
   
    % To get an actual matrix representation, you must specify
    % the row, column, and size variables to concrete integers
    % This is accomplished using the MatrixForm() method

    properties
       Row
       Column
       Entry
       Size
    end

    methods

        % Constructor
        function obj = SSEM(Row,Column,Size,Entry)
            obj.Row = Row;          % Can be symbolic
            obj.Column = Column;    % Can be symbolic
            obj.Entry = Entry;      % Can be symbolic, usually is
            obj.Size = Size;        % Must be a positive integer, NOT symbolic
        end

        function mat = MatrixForm(obj,row_sub,column_sub)
            if nargin == 1
                assert(isa(obj.Row,'numeric'))
                assert(isa(obj.Column,'numeric'))
                row_sub = obj.Row;
                column_sub = obj.Column;
            end
            assert(row_sub <= obj.Size)
            assert(column_sub <= obj.Size)
            mat = sym(zeros(obj.Size));
            mat(row_sub,column_sub) = obj.Entry;
        end

        function output = readableForm(obj)
            % Output a slightly more readable version
            % where the entry is converted to a string if it is an SKD
            output = obj;
            if isa(obj.Entry,'SKD')
                output.Entry = output.Entry.symForm();
            end
        end

        function negative = uminus(obj)
            negative = SSEM(obj.Row,obj.Column,obj.Size,-obj.Entry);
        end

        function SCMversion = convertToSCM(obj)
            SCMversion = SCM(obj.Size,zeros(obj.Size),{obj});
        end

        function output_array = rowForm(obj, with_size)
            % Create an array version of the SSEM object
            % This is an array of length 3, which contains the four entries:
            %   Row, Column, Entry (in this order)
            % with the entry converted to "symbolic form"
            % if it is an SKD object, which it very often is

            % If with_size is 0 or omitted, do the above. If with_size is 1,
            % then include a column for matrix size.
            if nargin == 1
                with_size = 0;
            end


            if with_size
                num_columns = 4;
            else
                num_columns = 3;
            end
            output_array = sym(zeros(1,num_columns));

            output_array(1) = obj.Row;
            output_array(2) = obj.Column;

            if isa(obj.Entry,'SKD')
                output_array(3) = obj.Entry.symForm();
            else
                output_array(3) = obj.Entry;
            end

            if with_size
                output_array(4) = obj.Size;
            end
        end

        function new_obj = addDerivedConditions(obj)
            if isa(obj.Entry,'SKD')
                new_entry = obj.Entry.addDerivedConditions();
            else
                new_entry = obj.Entry;
            end
            new_obj = SSEM(obj.Row,obj.Column,obj.Size,new_entry);
        end

    end

    methods (Static)

        function product = multiply(mat1, mat2)
            % Multiply two symbolic single entry matrices
            assert(mat1.Size == mat2.Size)
            new_condition = (mat1.Column == mat2.Row);
            if ~isa(mat1.Entry,'SKD') && ~isa(mat2.Entry,'SKD')
                if isAlways(new_condition,Unknown="false")
                    % If the column and row are always equal, then
                    % the new entry does not need to be an SKD object
                    new_entry = mat1.Entry*mat2.Entry;
                else
                    new_entry = SKD(new_condition,mat1.Entry*mat2.Entry);
                end
            else
                new_entry = SKD.multiply(mat1.Entry,mat2.Entry);
                new_entry = new_entry.addCondition(new_condition);
            end
            product = SSEM(mat1.Row,mat2.Column,mat1.Size,new_entry);
        end

        function runTests()
            fprintf("Running SSEM tests...");

            syms a
            syms b
            syms c
            syms d
            syms u
            syms v

            n = 5;
            E_ab_u = SSEM(a,b,n,u);
            E_ab_u_mat_form = E_ab_u.MatrixForm(1,2);
            E_ab_u_mat = sym(zeros(n));
            E_ab_u_mat(1,2) = u;
            assert(isequal(E_ab_u_mat_form,E_ab_u_mat))

            negative_E_ab_u = -E_ab_u;
            negative_E_ab_u_mat_form = negative_E_ab_u.MatrixForm(1,2);
            assert(isequal(negative_E_ab_u_mat_form,-E_ab_u_mat))

            E_cd_v = SSEM(c,d,n,v);
            prod = SSEM.multiply(E_ab_u,E_cd_v);
            del = SKD(b==c,u*v);
            assert(isequal(prod.Entry,del))

            fprintf(" all tests passed.\n")
        end

    end

end