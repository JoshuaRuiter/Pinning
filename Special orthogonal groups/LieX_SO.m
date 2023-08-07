function mat = LieX_SO(MatrixSize, Root_System, Form, alpha, v)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    n = MatrixSize;
    q = Root_System.Rank;
    diff = n-2*q;
    mat = sym(zeros(MatrixSize));
    c = Form.AnisotropicPartVector;

    assert(Root_System.IsRoot(alpha));
    assert(MatrixSize > 2*q);
    assert(length(v)==RootSpaceDimensionSO(MatrixSize, Root_System, alpha))
    
    %sum alpha to find out the type, possible sums are -2,-1,0,1,2
    type = sum(alpha);
 
    % Short root
    % has root space of dimension diff=n-2q
    if type == 1 || type == -1
        % alpha is a short root, of the form +/- alpha_i
        assert(length(v)==diff);
        for i = 1:q
            if alpha(i) == 1
                for s=1:diff
                    mat(i,2*q+s) = -c(s)*v(s);
                    mat(2*q+s,q+i) = v(s);
                end
                break
            elseif alpha(i) == -1
                for s=1:diff
                    mat(q+i,2*q+s) = -c(s)*v(s);
                    mat(2*q+s,i) = v(s);
                end
                break
            end
         end

    % Long root
    % has root space with dimension 1
    elseif type == 0
        % alpha is a long root of the form alpha_i-alpha_j
        assert(length(v)==1);
        a = 0;
        b = 0;
        for i = 1:q
            if alpha(i) == 1
                a = i;
            elseif alpha(i) == -1
                b = i;
            end
        end
        mat(a,b) = v(1);
        mat(q+b, q+a) = -v(1);

    elseif type == 2
        assert(length(v)==1);

        % This accomplishes the code below in one line
        a = find(alpha==1);
%         a = [0,0];
%         index = 1;
%         for i = 1:q
%             if alpha(i) == 1
%                 a(index) = i;
%                 index = index + 1;
%             end
%         end

        mat(a(1),q+a(2)) = v(1);
        mat(a(2),q+a(1)) = -v(1);

    elseif type == -2
        assert(length(v)==1);

        % This accomplishes the code below in one line
        a = find(alpha==-1);
%         a = [0,0];
%         index = 1;
%         for i = 1:q
%             if alpha(i) == -1
%                 a(index) = i;
%                 index = index + 1;
%             end
%         end

        mat(q+a(1),a(2)) = v(1);
        mat(q+a(2),a(1)) = -v(1);
    end    
end