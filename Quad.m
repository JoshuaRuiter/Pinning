classdef Quad

    % A class to represent an element of a quadratic field extension

    properties
        Real
        Imag
        PrimitiveElement
    end

    methods
        
        % Constructor
        function obj = Quad(Real, Imag, PrimitiveElement)
            obj.PrimitiveElement = PrimitiveElement;
            obj.Real = Real;
            obj.Imag = Imag;
        end

        % Overloading the + operator
        function sum = plus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            sum = Quad(a.Real+b.Real,a.Imag+b.Imag,a.PrimitiveElement);
        end

        % Overloading the - operator
        function difference = minus(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            difference = Quad(a.Real-b.Real,a.Imag-b.Imag,a.PrimitiveElement);
        end

        % Overloading the * operator
        % (entry-wise multiplication)
        function product = times(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            product = Quad(a.Real*b.Real + a.Imag*b.Imag*(a.PrimitiveElement)^2,...
                a.Real*b.Imag + a.Imag*b.Real,a.PrimitiveElement);
        end

        % Overloading the .* operator
        % (same as * for Quad type)
        function product = mtimes(a, b)
            assert(isequal(a.PrimitiveElement,b.PrimitiveElement))
            product = Quad(a.Real*b.Real + a.Imag*b.Imag*(a.PrimitiveElement)^2,...
                a.Real*b.Imag + a.Imag*b.Real,a.PrimitiveElement);
        end

        % I think this overloads the basic conj operation
        function conjugation = conj(obj)
            conjugation = Quad(obj.Real,-obj.Imag,obj.PrimitiveElement);
        end

        % Overloading the == operator
        function bool = isequal(a,b)
            bool = (isequal(a.PrimitiveElement,b.PrimitiveElement) ...
                    && isequal(a.Real,b.Real) ...
                    && isequal(a.Imag,b.Imag));
        end

        function num = ComplexForm(obj)
            num = obj.Real + obj.Imag*obj.PrimitiveElement;
        end

    end

end
