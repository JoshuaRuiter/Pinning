function my_matrix = GetB(n,q,vec_C)
    % Given positive integers n and q
    % and a vector vec_C = [c1, c2, c3, ... ]
    % output the matrix B of a nondegenerate isotropic
    % symmetric bilinear form on a vector space of dimension n
    % B is a square matrix with 3x3 block structure:
    %       (1,3) block -- Top right corner is (q x q) identity
    %       (2,2) block -- Middle block is diag(vec_C)
    %       (3,1) block -- Bottom left block is (q x q) identity
    %       All other blocks are zero

    assert(length(vec_C)==n-2*q);

    my_matrix = SymbolicZeros(n);
    % Put ones in for the two identity blocks
    for i=1:q
        my_matrix(i,n-q+i) = 1; % (1,3) identity block
        my_matrix(n-q+i,i) = 1; % (3,1) identity block
    end

    % Put vec_C along the diagonal ofthe (2,2) block
    for i=1:(n-2*q)
        my_matrix(i+q,i+q) = vec_C(i);
    end
end