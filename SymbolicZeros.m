function mat = SymbolicZeros(n)
    % Build a symbolic matrix whose entries are all zero

    % This does seem to be necessary in some cases, rather than just
    % using the built-in Matlab function zeros(n), because 
    % when you use zeros(n) and set an entry to be a symbolic variable
    % and then multiply symbolic matrices together, Matlab
    % seems to do unexpected things when comparing such matrices
    % to check if they are equal.

    mat = sym('mat',[n n]);
    for i=1:n
        for j=1:n
            mat(i,j) = 0;
        end
    end
end