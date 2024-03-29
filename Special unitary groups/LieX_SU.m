function mat = LieX_SU(MatrixSize, Root_System, NIForm, alpha, u)
    % Takes inputs alpha (a root), root_system the root system it 
    % comes from, and u (a vector, possibly symbolic)
    % Output the associated element of the Lie algebra

    n = MatrixSize;
    q = Root_System.Rank;
    H = NIForm.Matrix;
    C = NIForm.AnisotropicMatrix;
    eps = NIForm.Epsilon;
    P = NIForm.PrimitiveElement;

    if eps == 1
        P_eps = P;
    else % eps == -1
        P_eps = 1;
    end

    % validate inputs - alpha is a root, and u is the right length
    assert(Root_System.IsRoot(alpha))
    assert(length(u) == RootSpaceDimensionSU(n,Root_System,alpha))
    assert(length(C) == n-2*q)

    % create a nxn matrix of zeros, that allows symbolic
    mat = sym(zeros(n));
    
    if IsLong(alpha)
        % In this case, the root alpha is of the form 
        % +/- 2alpha_i for some i
        i = find(alpha~=0);
        
        % In this case, the input should be a scalar 
        % (a single element of the base field k)
        assert(RootSpaceDimensionSU(n,Root_System,alpha)==1)
        assert(length(u)==1);
        
        if sum(alpha) == 2
            % alpha = 2alpha_i
            mat(i,i+q) = u*P_eps;
        elseif sum(alpha) == -2
            % alpha = -2alpha_i
            mat(i+q,i) = u*P_eps;
        else
            % something is wrong, throw an error
            printf("A long root has the wrong form.");
            assert(false);
        end

    elseif IsShort(alpha)
        % In this case, the root alpha is of the form 
        % +/- alpha_i for some i
        i = find(alpha~=0);

        % In this case, the input should have length 2*(n-2*q)
        assert(RootSpaceDimensionSU(n,Root_System,alpha)==2*(n-2*q))
        assert(length(u)==2*(n-2*q))

        % Convert u to a vector of length n-2*q, with complex entries
        u_complex = sym(zeros(1,n-2*q));
        for j=1:n-2*q
            u_complex(j) = u(2*j-1) + u(2*j)*P;
        end

        if sum(alpha) == 1
            % alpha = alpha_i
            for j=1:n-2*q
                mat(2*q+j,q+i) = u_complex(j);
                mat(i,2*q+j) = -eps*C(j,j)*conjugate(u_complex(j),P);
            end

        elseif sum(alpha) == -1
            % alpha = -alpha_i
            for j=1:n-2*q
                mat(2*q+j,i) = u_complex(j);
                mat(q+i,2*q+j) = -C(j,j)*conjugate(u_complex(j),P);
            end

        else
            % something is wrong, throw an error
            printf("A short root has the wrong form.");
            assert(false);
        end

    elseif IsMedium(alpha)

        % In this case, the input should vector of length 2,
        % which carries the "real" part in the first entry 
        % and the "imaginary" part in the second entry
        assert(length(u)==2);
        
        % switch to viewing u as a single quadratic extension element
        % P is the primitive element
        u_complex = complexify(u,P);
        u_conjugate = conjugate(u_complex,P);

        if sum(alpha)==0
            % alpha = alpha_i - alpha_j
            i = find(alpha==1);
            j = find(alpha==-1);

            mat(i,j) = u_complex; 
            mat(j+q,i+q) = -u_conjugate;
            
        elseif sum(alpha)==2
            % alpha = alpha_i + alpha_j, with i and j distinct
            % convention: i is less than j

            % F finds the first two indices of alpha where the value is 1
            % and puts them into a 1x2 array
            F = find(alpha==1,2);
            i = F(1);
            j = F(2);
            assert(i<j);

            mat(i,j+q) = u_complex;
            mat(j,i+q) = -eps*u_conjugate;
            
        elseif sum(alpha)==-2
            % alpha = -alpha_i - alpha_j, with i and j distinct
            % convention: i is less than j

            % F finds the first two indices of alpha where the value is -1
            % and puts them into a 1x2 array
            F = find(alpha==-1,2); 
            i = F(1);
            j = F(2);
            assert(i<j);

            mat(j+q,i) = u_complex;
            mat(i+q,j) = -eps*u_conjugate;

        else
            % something is wrong, throw an error
            printf("A normal root isn't of the right form.");
            assert(false);
        end
    else
        % Something is wrong, throw an error
        printf("A root is neither long, short, nor medium.");
        assert(false);
    end

    assert(isequal(size(mat),[n,n]));
end

function bool = IsShort(alpha)
    bool = (dot(alpha,alpha)==1);
end
function bool = IsMedium(alpha)
    bool = (dot(alpha,alpha)==2);
end
function bool = IsLong(alpha)
    bool = (dot(alpha,alpha)==4);
end