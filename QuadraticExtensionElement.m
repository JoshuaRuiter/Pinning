classdef QuadraticExtensionElement

    % A class to represent an element of a quadratic field extension
    % If a = Coord1 and b = Coord2, and c = PrimitiveElement, then this
    % object represents a+b*c

    properties
        Real
        Imag
        PrimitiveElement
    end

    methods
        
        % Constructor
        function obj = QuadraticExtensionElement(Real, Imag, PrimitiveElement)
            obj.PrimitiveElement = PrimitiveElement;
            obj.Real = Real;
            obj.Imag = Imag;
        end

        % Overloading the + operator
        function sum = plus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            sum.PrimitiveElement = a.PrimitiveElement;
            sum.Real = a.Real + b.Real;
            sum.Imag = a.Imag + b.Imag;
        end

        % Overloading the - operator
        function difference = minus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            difference.PrimitiveElement = a.PrimitiveElement;
            difference.Real = a.Real - b.Real;
            difference.Imag = a.Imag - b.Imag;
        end

        % Overloading the * operator
        function product = times(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            product.PrimitiveElement = a.PrimitiveElement;
            product.Real = a.Real*b.Real + a.Imag*b.Imag*(product.PrimitiveElement)^2;
            product.Imag = a.Real*b.Imag + a.Imag*b.Real;
        end

        % Overloading the .* operator        
        function product = mtimes(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            product.PrimitiveElement = a.PrimitiveElement;
            product.Real = a.Real*b.Real + a.Imag*b.Imag*(product.PrimitiveElement)^2;
            product.Imag = a.Real*b.Imag + a.Imag*b.Real;
        end

        % I think this overloads the basic conj operation
        function conjugation = conj(a)
            conjugation.PrimitiveElement = a.PrimitiveElement;
            conjugation.Real = a.Real;
            conjugation.Imag = -a.Imag;
        end

        % Overloading the == operator
        function bool = eq(a,b)
            bool = (a.PrimitiveElement == b.PrimitiveElement ...
                    && a.Real == b.Real ...
                    && a.Imag == b.Imag);
        end

        function bool = IsReal(a)
            bool = (a.Imag==0);
        end

        function bool = IsImag(a)
            bool = (a.Real==0);
        end

        function num = ComplexForm(a)
            num = a.Real + a.Imag*a.PrimitiveElement;
        end

        function mat = Display(quad_mat)
            row_col = size(quad_mat);
            rows = row_col(1);
            cols = row_col(2);
            mat = sym(zeros(rows,cols));
            for i=1:rows
                for j=1:cols
                    mat(i,j) = ComplexForm(quad_mat(i,j));
                end
            end
        end

    end

    methods (Static)

        function mat = zeros(n,PrimitiveElement)
            % Create a square n-by-n matrix whose entries have the
            % QuadraticExtensionElement type
            z = QuadraticExtensionElement(0,0,PrimitiveElement);
            mat = repelem(z,n,n);
        end


    end
end
