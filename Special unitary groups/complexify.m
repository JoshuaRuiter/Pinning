function my_scalar = complexify(vector,PrimitiveElement)
    assert(length(vector)==2);
    my_scalar = vector(1) + vector(2)*PrimitiveElement;
end