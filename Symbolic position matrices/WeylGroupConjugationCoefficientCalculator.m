% This script is designed to help calculate
% conjugations of Weyl group elements

function WeylGroupConjugationCoefficientCalculator()
    % Do the calculation for SL3
    % The 3 doesn't really matter, the calculation
    % will be exactly the same for any n that you put
    % into this. The 3 only matters when you try and 
    % visualize the resulting matrices.
    WeylGroupConjugationCoefficientCalculatorSL(3);

end

function WeylGroupConjugationCoefficientCalculatorSL(n)    
    % Calculate the result of a Weyl group conjugation
    % associated to any root in A_(n-1)
    % of any root subgroup in SL_n

    % In other terms, calculate w_alpha(1) * X_beta(v) * w_alpha(1)^(-1)
    % where alpha and beta are any two roots
    
    % Let alpha = alpha_a - alpha_b
    % and beta = alpha_c - alpha_d

    % Create w_alpha_1  
    syms a
    syms b
        
        % Create X_alpha_1
        E_ab_1 = SSEM(a,b,n,1);
        e_ab_1 = SCM(n,eye(n),{E_ab_1});
        X_alpha_1 = e_ab_1;

        % Create X_m_alpha_m_1
        E_ba_m_1 = SSEM(b,a,n,-1);
        e_ba_m_1 = SCM(n,eye(n),{E_ba_m_1});
        X_m_alpha_m_1 = e_ba_m_1;

        % Now we can build w_alpha_1
        temp1 = SCM.multiply(X_alpha_1,X_m_alpha_m_1);
        w_alpha_1 = SCM.multiply(temp1, X_alpha_1);

    % Create w_alpha_1_inverse
    
        % Create X_alpha_m_1
        E_ab_m_1 = SSEM(a,b,n,-1);
        e_ab_m_1 = SCM(n,eye(n),{E_ab_m_1});
        X_alpha_m_1 = e_ab_m_1;

        % Create X_m_alpha_1
        E_ba_1 = SSEM(b,a,n,1);
        e_ba_1 = SCM(n,eye(n),{E_ba_1});
        X_m_alpha_1 = e_ba_1;

        % Now we can build w_alpha_1_inverse
        temp2 = SCM.multiply(X_alpha_m_1, X_m_alpha_1);
        w_alpha_1_inverse = SCM.multiply(temp2,X_alpha_m_1);

    % Create X_alpha(v)
    syms c
    syms d
    syms v
    E_cd_v = SSEM(c,d,n,v); 
    e_cd_v = SCM(n,eye(n),{E_cd_v});
    X_alpha_v = e_cd_v;

    % Compute the conjugate
    temp3 = SCM.multiply(w_alpha_1,X_alpha_v);
    conjugation = SCM.multiply(temp3,w_alpha_1_inverse);

    % Cancel out negative pairs
    sim_con_1 = conjugation.cancelPairs();

    % Remove terms that are zero when a~=b or c~=d, 
    % after first adding in derived conditions
    sim_con_2 = sim_con_1.addDerivedConditions();
    sim_con_3 = sim_con_2.simplifyWithAssumption(a~=b);
    sim_con_4 = sim_con_3.simplifyWithAssumption(c~=d);

    % Output a readable table
    output_table = sim_con_4.tableForm()

end