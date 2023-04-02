function bool = SymbolicIsEqual(A,B)
    % Compare matrices A and B, 
    % where potentially one or both has symbolic values

    % default result is true
    bool = true;

    % return false if A and B have different sizes
    if length(A) ~= length(B)
        bool = false;
        return
    end

    for i=1:length(A)
        for j=1:length(A)
            if ~isAlways(A(i,j)==B(i,j))
                % if the (i,j) entry of A and the (i,j) entry of B are
                % not equal, return false
                bool = false;

                % If we get a single false result, we can just escape the
                % function altogether, we don't need to complete the for
                % loop
                return
            end
        end
    end
end