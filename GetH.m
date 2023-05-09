function myMatrix = GetH(n,q,vec_C)
    % return the (n x n) matrix H used to define SU_{n,q}(L,h)
    % H is a 3x3 block matrix with blocks of 
    %       (q x q) identity in the (1.2) block, 
    %       (q x q) negative identity in the (2,1) block
    %       diag(vec_C) in the (3,3) block

    myMatrix = sym(zeros(n));
    assert(length(vec_C)==n-2*q);

    % create the two identity blocks
    for i=1:q
        % this makes the (q x q) identity block in the (1,2) block
        myMatrix(i,q+i) = 1;

        % this makes the (q x q) negative identity block in the (2,1) block
        myMatrix(q+i,i) = -1;
    end

    % this creates the diag(vec_C) block in the (3,3) block
    if n > 2*q
        for i=1:n-2*q
            myMatrix(2*q+i,2*q+i) = vec_C(i);
        end
    end
end