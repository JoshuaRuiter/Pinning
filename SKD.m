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

        % Overloading the unary - operator
        function negative = uminus(obj)
            negative = SKD(obj.Conditions,-obj.ValueIfTrue);
        end

        function my_string = conditionString(obj)
            % Output a string version of the conditions list
            my_string = string(obj.Conditions(1));
            for i=2:length(obj.Conditions)
                my_string = strcat(my_string,"_",string(obj.Conditions(i)));
            end
            my_string = strrep(my_string,"=","");
            my_string = strrep(my_string," ","");
            my_string = strrep(my_string,"*","");
        end

        function my_string = strForm(obj)
            % Output a string version of the SKD
            condition_string = obj.conditionString();
            my_string = strcat("del_", condition_string,"_",string(obj.ValueIfTrue));

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

        function new_obj = addCondition(obj, new_condition)
            new_obj = obj;
            if ~isAlways(new_condition,Unknown="false") && ~ismember(new_condition, new_obj.Conditions)
                new_obj.Conditions(end+1) = new_condition;
            end
        end

        function new_obj = addFirstOrderDerivedConditions(obj)
            % Build a new copy of the SKD object, which includes
            % derived conditions from pairs of existing conditions
            % For example, if the original object has conditions a==b and b==c,
            % then the new object should also have these conditions, with the
            % addition derived condition a==c
            new_obj = obj;

            % Iterate over all pairs, adding derived conditions as needed
            for i=1:length(obj.Conditions)
                eq_i = obj.Conditions(i);
                for j=i+1:length(obj.Conditions)
                    eq_j = obj.Conditions(j);
                    if isAlways(lhs(eq_i)==lhs(eq_j),Unknown="false")
                        new_obj = new_obj.addCondition(rhs(eq_i)==rhs(eq_j));
                    elseif isAlways(rhs(eq_i)==lhs(eq_j),Unknown="false")
                        new_obj = new_obj.addCondition(lhs(eq_i)==rhs(eq_j));
                    elseif isAlways(lhs(eq_i)==rhs(eq_j),Unknown="false")
                        new_obj = new_obj.addCondition(rhs(eq_i)==lhs(eq_j));
                    elseif isAlways(rhs(eq_i)==rhs(eq_j),Unknown="false")
                        new_obj = new_obj.addCondition(lhs(eq_i)==lhs(eq_j));
                    end
                end
            end
        end

        function new_obj = addDerivedConditions(obj)
            % Build a new copy of the SKD object, which includes
            % derived conditions from pairs of existing conditions,
            % and derived conditions from newly added conditions, etc.

            new_obj = obj;

            % Iterate adding first order derived conditions
            % until no more are added. At that point, stop.
            while true
                starting_conditions = length(new_obj.Conditions);
                new_obj = new_obj.addFirstOrderDerivedConditions();
                ending_conditions = length(new_obj.Conditions);
                if starting_conditions == ending_conditions
                    % If no new conditions were added, then break the loop.
                    break;
                end
            end
        end

    end

    methods (Static)

        function new_conditions = combineConditions(obj1,obj2)
            % Combine two lists of conditions, without repeats
            new_conditions = obj1.Conditions;
            for i=1:length(obj2.Conditions)
                condition_i = obj2.Conditions(i);
                if ~isAlways(condition_i,Unknown="false") && ~ismember(condition_i,obj1.Conditions)
                    new_conditions(end+1) = obj2.Conditions(i);
                end
            end
        end

        function prod = multiply(obj1,obj2)
            % Mutliply two SymbolicKroneckerDelta objects together
            % or a SymbolicKroneckerDelta object with another variable type
            % resulting in a new SymbolicKroneckerDelta object

            if isa(obj1,'SKD') && isa(obj2,'SKD')
                % Note that this is dependent on the assumption that 
                % the value when false is zero for at least one,
                % but we have built in this assumption so that the
                % value when false is zero for both
                new_conditions = SKD.combineConditions(obj1,obj2);
                new_value_if_true = obj1.ValueIfTrue*obj2.ValueIfTrue;
            elseif isa(obj1,'SKD') && ~isa(obj2,'SKD')
                new_conditions = obj1.Conditions;
                new_value_if_true = obj1.ValueIfTrue*obj2;
            elseif ~isa(obj1,'SKD') && isa(obj2,'SKD')
                new_conditions = obj2.Conditions;
                new_value_if_true = obj2.ValueIfTrue*obj1;
            else 
                % This should be impossible
                assert(false,'Trying to use the SKD ' + ...
                    'multiplication with two non-SKD types.')
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
            
            del_ab = SKD(a==b,1);
            del_cd = SKD(c==d,1);
            prod_1 = SKD.multiply(del_ab,del_cd);
            del_3 = SKD([a==b,c==d],1);
            assert(isequal(prod_1,del_3))

            neg_del_ab = -del_ab;
            neg_del_ab_2 = SKD(a==b,-1);
            assert(isequal(neg_del_ab,neg_del_ab_2))

            prod_2 = SKD.multiply(del_ab,u);
            u_del_ab = SKD(a==b,u);
            assert(isequal(prod_2.Conditions,u_del_ab.Conditions))
            assert(isequal(prod_2.ValueIfTrue,u_del_ab.ValueIfTrue))
            assert(isequal(prod_2,u_del_ab))

            prod_3 = SKD.multiply(u,del_ab);
            assert(isequal(prod_3.Conditions,u_del_ab.Conditions))
            assert(isequal(prod_3.ValueIfTrue,u_del_ab.ValueIfTrue))
            assert(isequal(prod_3,u_del_ab))

            del_ab = SKD(a==b,1);
            del_bc = SKD(b==c,1);
            del_ab_bc = SKD([a==b,b==c],1);
            prod_4 = SKD.multiply(del_ab,del_bc);
            assert(length(del_ab_bc.Conditions)==2)
            assert(isequal(del_ab_bc,prod_4))
            del_abc = del_ab_bc.addFirstOrderDerivedConditions();
            assert(length(del_abc.Conditions)==3)
            del_abc_2 = del_ab_bc.addDerivedConditions();
            assert(length(del_abc_2.Conditions)==3)

            del_1 = SKD([a==b,b==c,c==d],1);
            assert(length(del_1.Conditions)==3)
            del_2 = del_1.addFirstOrderDerivedConditions();
            assert(length(del_2.Conditions)==5)
            del_3 = del_2.addFirstOrderDerivedConditions();
            assert(length(del_3.Conditions)==6)
            del_4 = del_3.addFirstOrderDerivedConditions();
            assert(length(del_4.Conditions)==6)
            del_5 = del_4.addDerivedConditions();
            assert(length(del_5.Conditions)==6)

            fprintf(" all tests passed.\n")
        end

    end

end