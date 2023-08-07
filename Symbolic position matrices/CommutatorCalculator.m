% This script is designed to help calculate
% commutators of elementary or near-elementary 
% matrices

function CommutatorCalculator()
    % Do the calculation for SL3
    % The 3 doesn't really matter, the calculation
    % will be exactly the same for any n that you put
    % into this. The 3 only matters when you try and 
    % visualize the resulting matrices.
    CommutatorCalculatorSL(3);

    % Do the calculation for SU_{6,2}
    % INCOMPLETE
%     CommutatorCalculatorSU(6,2);
end

function CommutatorCalculatorSL(n)
    % Calculate the commutator between two 
    % elementary matrices, i.e. pinning matrices
    % for SL_n

    fprintf("Computing the commutator of two aribtrary pinning \n" + ...
        "matrices for SL_n, namely e_ab(u) and e_cd(v). \n" + ...
        "A row in the table below denotes a single entry \n" + ...
        "matrix (additive) term, with the first two columns \n" + ...
        "giving the position and the third column giving \n" + ...
        "the entry. Not depicted in the table is an implicit \n" + ...
        "identity matrix additive term.");

    syms a
    syms b
    syms c
    syms d
    syms u
    syms v
    
    % Create two arbitrary symbolic single-entry matrices
    E_ab_u = SSEM(a,b,n,u);
    E_cd_v = SSEM(c,d,n,v); 

    % Create two arbitrary elementary matrices
    e_ab_u = SCM(n,eye(n),{E_ab_u});
    e_cd_v = SCM(n,eye(n),{E_cd_v});

    % Compute their commutator
    com = SCM.commutator(e_ab_u,e_cd_v);

    % Cancel out pairs of matrices that add to zero
    sim_com_1 = com.cancelPairs();

    % Remove terms dependent on a==b or c==d
    sim_com_2 = sim_com_1.simplifyWithAssumption(a~=b);
    sim_com_3 = sim_com_2.simplifyWithAssumption(c~=d);

    % This builds a table which can be exported
    output_table = sim_com_3.tableForm()
end

function CommutatorCalculatorSU(n,q)
    % Calculate the commutator between two 
    % pinning matrices for SL_{n,q}
    % INCOMPLETE
end
