function conjugated_element = conjugate(thing_to_conjugate, PrimitiveElement)
    % Compute a conjugate of a quadratic extension element
    conjugated_element = subs(thing_to_conjugate, PrimitiveElement, -PrimitiveElement);
end