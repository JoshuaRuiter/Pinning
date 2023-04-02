function bool = SymbolicIsDiag(A)
    % Take an input matrix A
    % Return true if A is a diagonal matrix
    bool = true;
    for i=1:length(A)
        for j=1:length(A)
            if i~=j && A(i,j)~=0
                % If an off-diagonal (i,j) entry is nonzero, then the
                % matrix is not diagonal
                bool = false;

                % If we get a single false result, we can just escape the
                % function altogether, we don't need to complete the for
                % loop
                return
            end
        end
    end
end