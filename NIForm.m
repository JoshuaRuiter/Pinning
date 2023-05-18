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
            % Should be a symbolic variable, or perhaps the imaginary unit 1i
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

            if strcmpi(NameString,'symmetric bilinear')
                assert(eps == 1)
            elseif not(isempty(vec_C))
                % (Skew)-hermitian, non-quasisplit case
                % C must satisfy vec_C = eps*conj(vec_C)
                assert(isequal(conjugate(vec_C,PrimitiveElement),eps*vec_C))
            end

            obj.Dimension = n;
            obj.Index = q;
            obj.Epsilon = eps;
            obj.AnisotropicPartVector = vec_C;
            obj.AnisotropicMatrix = diag(vec_C);
            obj.PrimitiveElement = PrimitiveElement;
            obj.NameString = NameString;

            % Now we build the associated matrix  
            % The constructed matrix is a 3x3 block matrix with blocks of 
            %       (q x q) identity in the (1,2) block, 
            %       (q x q) Epsilon*identity in the (2,1) block
            %       diag(AnisotropicPartVector) in the (3,3) block
            % Create a matrix of zeros of the right size
            obj.Matrix = sym(zeros(n));
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

            if strcmpi(NameString,'hermitian')
                assert(obj.IsHermitian());
            elseif strcmpi(NameString,'skew-hermitian')
                assert(obj.IsSkewHermitian());
            elseif strcmpi(NameString,'symmetric bilinear')
                assert(obj.IsSymmetric());
            end

        end

        function bool = IsSymmetric(obj)
            bool = isequal(obj.Matrix,transpose(obj.Matrix));
        end
        function bool = IsHermitian(obj)
            bool = isequal(obj.Matrix,transpose(conjugate(obj.Matrix,obj.PrimitiveElement)));
        end
        function bool = IsSkewHermitian(obj)
            bool = isequal(obj.Matrix,-transpose(conjugate(obj.Matrix,obj.PrimitiveElement)));
        end

    end

end