function mat = Commutator(A,B)
    % Given two square matrices of the same size, 
    % return the commutator matrix [A,B] = A*B*A^(-1)*B^(-1)
    mat = A*B*A^(-1)*B^(-1);
end