classdef QuadMat

    % A class to represent a matrix of elements from quadratic field extension

    properties
        RealMat
        ImagMat
        PrimitiveElement
    end

   methods
        
        % Constructor
        function obj = QuadMat(RealMat, ImagMat, PrimitiveElement)
            assert(isequal(size(RealMat),size(ImagMat)));
            obj.PrimitiveElement = PrimitiveElement;
            obj.RealMat = RealMat;
            obj.ImagMat = ImagMat;
        end

        % Overloading the + operator
        function sum = plus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            sum = QuadMat(a.RealMat + b.RealMat, a.ImagMat + b.ImagMat, a.PrimitiveElement);
        end

        % Overloading the - operator
        function difference = minus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            difference = QuadMat(a.RealMat - b.RealMat, a.ImagMat - b.ImagMat, a.PrimitiveElement);
        end

        % Overloading the .* operator
        % (element-wise multiplication)
        function product = times(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            P = a.PrimitiveElement;

            assert(isequal(size(a),size(b)));
            dimensions = size(a);
            rows = dimensions(1);
            cols = dimensions(2);
            real_prod = zeros(dimensions);
            imag_prod = zeros(dimensions);
            for i=1:rows
                for j=1:cols
                    real_prod(i,j) = a.RealMat(i,j)*b.RealMat(i,j) + a.ImagMat(i,j)*b.ImagMat(i,j)*P^2;
                    imag_prod(i,j) = a.RealMat(i,j)*b.ImagMat(i,j) + a.ImagMat(i,j)*b.RealMat(i,j);
                end
            end
            product = QuadMat(real_prod,imag_prod,P);
        end

        % Overloading the * operator
        % (matrix multiplication)
        function product = mtimes(a,b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            P = a.PrimitiveElement;
            
            size_a = size(a);
            size_b = size(b);
            assert(size_a(2) == size_b(1))
            common_size = size_a(2);

            rows = size_a(1);
            cols = size_b(2);
            dimensions = [rows,cols];

            real_prod = zeros(dimensions);
            imag_prod = zeros(dimensions);
            for i=1:rows
                for j=1:cols
                    for k=1:common_size
                        real_prod(i,j) = real_prod(i,j) + ...
                            a.RealMat(i,k)*b.RealMat(k,j) + a.ImagMat(i,k)*b.ImagMat(k,j)*P^2;
                        imag_prod(i,j) = imag_prod(i,j) + ...
                            a.RealMat(i,k)*b.ImagMat(k,j) + a.ImagMat(i,k)*b.RealMat(k,j);
                    end
                end
            end
            product = QuadMat(real_prod,imag_prod,P);
        end

        % Conjugation operation
        function conjugation = conj(obj)
            conjugation = QuadMat(obj.RealMat,-obj.ImagMat,obj.PrimitiveElement);
        end

        % Equality
        function bool = isequal(a,b)
            bool = (isequal(a.PrimitiveElement,b.PrimitiveElement) ...
                    && isequal(a.Real,b.Real) ...
                    && isequal(a.Imag,b.Imag));
        end

        % Output a complex matrix for easier visualization
        function mat = ComplexForm(a)
            mat = a.RealMat + a.ImagMat*a.PrimitiveElement;
        end

        % Output a vector [rows,columns]
        function dimensions = size(obj)
            dimensions = size(obj.RealMat);
        end

        % Overloading subscripted reference
        function quad = subsref(i,j)
            quad = QuadMat(obj.RealMat(i,j),obj.ImagMat(i,j),obj.PrimitiveElement);
        end

   end

   methods (Static)

       function zero_quad_mat = zeros(n,PrimitiveElement)
            zero_quad_mat = QuadMat(zeros(n),zeros(n),PrimitiveElement);
       end

   end

end