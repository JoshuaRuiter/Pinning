classdef SKD
    
    % Symbolic Kronecker Delta
    % A class to model the Kronecker delta function del_ij 
    % which is 1 when i==j and 0 otherwise
    % We allow slightly more general, to have a variable value when true
    % but still always be zero when false.

    properties
        Conditions
        ValueIfTrue
    end

    methods

        % Constructor
        function obj = SKD(conditions,val_if_true)
            obj.Conditions = conditions;
            obj.ValueIfTrue = val_if_true;
        end

        function negative = negate(obj)
            negative = SKD(obj.Conditions,-obj.ValueIfTrue);
        end

        function my_string = strForm(obj)
            % Output a string version of the SKD
            condition_string = "";
            for i=1:length(obj.Conditions)
                condition_string = strcat(condition_string,"_",string(obj.Conditions(i)));
            end
            my_string = strcat("del", condition_string,"_",string(obj.ValueIfTrue));

            % Remove = signs, spaces, and *
            my_string = strrep(my_string,"=","");
            my_string = strrep(my_string," ","");
            my_string = strrep(my_string,"*","");

            % Replace negative signs with "neg"
            my_string = strrep(my_string,"-","neg_");

            % Replace ^ with "pow"
            my_string = strrep(my_string,"^","pow");
        end

        function my_sym = symForm(obj)
            % Output a symbolic version of the SKD
            my_sym = sym(obj.strForm());
        end

    end

    methods (Static)

        function prod = multiply(obj1,obj2)
            % Mutliply two SymbolicKroneckerDelta objects together
            % or a SymbolicKroneckerDelta object with another variable type
            % resulting in a new SymbolicKroneckerDelta object

            if isa(obj1,'SKD') && isa(obj2,'SKD')
                % Note that this is dependent on the assumption that 
                % the value when false is zero
                new_conditions = [obj1.Conditions,obj2.Conditions];
                new_value_if_true = obj1.ValueIfTrue*obj2.ValueIfTrue;
            elseif isa(obj1,'SKD') && ~isa(obj2,'SKD')
                new_conditions = obj1.Conditions;
                new_value_if_true = obj1.ValueIfTrue*obj2;
            elseif ~isa(obj1,'SKD') && isa(obj2,'SKD')
                new_conditions = obj2.Conditions;
                new_value_if_true = obj2.ValueIfTrue*obj1;
            else 
                % This should be impossible
                assert(false,['Trying to use the SKD...' ...
                    ' multiplication with two non-SKD types.'])
            end
            
            prod = SKD(new_conditions,new_value_if_true);
        end

        function runTests()
            fprintf("Running SKD tests...");

            syms a
            syms b
            syms c
            syms d
            syms u
            syms v
            
            del_ab = SKD(a==b,1);
            del_cd = SKD(c==d,1);
            prod_1 = SKD.multiply(del_ab,del_cd);
            del_3 = SKD([a==b,c==d],1);
            assert(isequal(prod_1,del_3))

            prod_2 = SKD.multiply(del_ab,u);
            u_del_ab = SKD(a==b,u);
            assert(isequal(prod_2.Conditions,u_del_ab.Conditions))
            assert(isequal(prod_2.ValueIfTrue,u_del_ab.ValueIfTrue))
            assert(isequal(prod_2,u_del_ab))

            prod_3 = SKD.multiply(u,del_ab);
            assert(isequal(prod_3.Conditions,u_del_ab.Conditions))
            assert(isequal(prod_3.ValueIfTrue,u_del_ab.ValueIfTrue))
            assert(isequal(prod_3,u_del_ab))

            fprintf(" all tests passed.\n")
        end

    end
end