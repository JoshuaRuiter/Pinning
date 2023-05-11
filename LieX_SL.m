function mat = LieX_SL(MatrixSize, Root_System, FormMatrix, alpha, u)
    % Inputs: an integer n, a root alpha in the A_(n-1) root system (vector of length n)
    % and u (a scalar, possibly a symbolic variable)
    % Output the associated element of the Lie algebra sl_n
    
    % rowNumber is the row number of where the nonzero entry goes
    % columnNumber is the column of where the nonzero entry goes
    rowNumber = alpha==1;
    columnNumber = alpha==-1;

    mat = SymbolicZeros(MatrixSize);
    mat(rowNumber,columnNumber) = u;
end