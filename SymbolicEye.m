function mat = SymbolicEye(n)
    % Build a symbolic identity matrix whose entries are all zero

    % This does seem to be necessary in some cases, rather than just
    % using the built-in Matlab function zeros(n), because 
    % when you use zeros(n) and set an entry to be a symbolic variable
    % and then multiply symbolic matrices together, Matlab
    % seems to do unexpected things when comparing such matrices
    % to check if they are equal.
    
    mat = sym('mat',[n n]);
    for i=1:n
        for j=1:n
            if i==j
                mat(i,j) = 1;
            else
                mat(i,j) = 0;
            end
        end
    end
end