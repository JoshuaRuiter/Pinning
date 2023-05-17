classdef NIForm

    % A class representing a nondegenerate isotropic form
    %   Can be hermitian, skew-hermitian, or symmetric bilinear

    properties
        Dimension % often called n
        Index % Witt index, often called q
        NameString % String, 'hermitian' or 'skew-hermitian' or 'symmetric bilinear'
        Epsilon 
            % Epsilon = 1 for hermitian or symmetric bilinear
            % Epsilon = -1 for skew-hermitian
        Matrix % often called H
        AnisotropicPartVector % often called vec_C or Vector_C
        AnisotropicMatrix % often called C
        PrimitiveElement 
            % Primitive element of associated quadratic extension
            % Only relevant for hermitian and skew-hermitian forms
    end

    methods

        % Constructor
        function obj = NIForm(n,q,eps,vec_C,PrimitiveElement,NameString)
            % Construct a form associated to a (skew)-hermitian form on a
            % vector space of dimension n, with Witt index q, using eps=1
            % for hermitian and eps=-1 for skew-hermitian, with vec_C
            % listing the diagonal values of the anisotropic part

            % Validate inputs
            assert(n >= 2*q) % Witt index cannot exceed half the dimension
            assert(length(vec_C)==n-2*q) % vec_C must have a certain length
            assert(abs(eps)==1) % epsilon is 1 or -1
            assert(strcmpi(NameString,'hermitian') || ...
                    strcmpi(NameString,'skew-hermitian') || ...
                    strcmpi(NameString,'symmetric bilinear'));

            obj.Dimension = n;
            obj.Index = q;
            obj.Epsilon = eps;
            obj.AnisotropicPartVector = vec_C;
            obj.PrimitiveElement = PrimitiveElement;
            obj.NameString = NameString;

            % Now we build the associated matrix  
            % The constructed matrix is a 3x3 block matrix with blocks of 
            %       (q x q) identity in the (1,2) block, 
            %       (q x q) Epsilon*identity in the (2,1) block
            %       diag(AnisotropicPartVector) in the (3,3) block
            % Create a matrix of zeros of the right size
            if strcmpi(NameString,'hermitian') || strcmpi(NameString,'skew-hermitian')
                obj.Matrix = QuadMat.zeros(n,PrimitiveElement);
%                 one_quad = QuadraticExtensionElement(1,0,PrimitiveElement);
%                 eps_quad = QuadraticExtensionElement(eps,0,PrimitiveElement);

                % C must satisfy vec_C = eps*conj(vec_C)
                for i=1:n-2*q
                    assert(eq(vec_C(i),eps_quad*conj(vec_C(i))));
                end
                
                % create the two identity blocks
                for i=1:q
                    % this makes the (q x q) identity block in the (1,2) block
                    obj.Matrix(i,q+i) = one_quad;
            
                    % this makes the (q x q) Epsilon*identity block in the (2,1) block
                    obj.Matrix(q+i,i) = eps_quad;
                end
            
                % this creates the diag(vec_C) block in the (3,3) block
                if n > 2*q
                    for i=1:n-2*q
                        obj.Matrix(2*q+i,2*q+i) = vec_C(i);
                    end
                end

            elseif strcmpi(NameString,'symmetric bilinear')
                obj.Matrix = sym(zeros(n));
                
                % create the two identity blocks
                for i=1:q
                    % this makes the (q x q) identity block in the (1,2) block
                    obj.Matrix(i,q+i) = 1;
            
                    % this makes the (q x q) Epsilon*identity block in the (2,1) block
                    obj.Matrix(q+i,i) = eps;
                end
            
                % this creates the diag(vec_C) block in the (3,3) block
                if n > 2*q
                    for i=1:n-2*q
                        obj.Matrix(2*q+i,2*q+i) = vec_C(i);
                    end
                end

            else
                assert(false,['The label for a nondegenerate isotropic form is invalid. '...
                    'Valid labels are hermitian, skew-hermitian, and symmetric bilinear.'])
            end
        end

        function bool = IsSymmetric(obj)
            bool = isequal(obj.Matrix,transpose(obj.Matrix));
        end
        function bool = IsHermitian(obj)
            bool = isequal(obj.Matrix,ctranspose(obj.Matrix));
        end
        function bool = IsSkewHermitian(obj)
            bool = isequal(obj.Matrix,-ctranspose(obj.Matrix));
        end

    end

end