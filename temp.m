% Trying to use the SymbolicSingleEntryMatrix class
% to work out commutator coefficients and
% Weyl group conjugation coefficients
% for special linear groups, as a proof of concpet

% Symbolic Kronecker Delta tests
SKD.runTests()

% Symbolic Single Entry Matrix tests
SSEM.runTests()

% Symbolic Coordinates Matrix tests
SCM.runTests()

% General commutator in SL_5
n = 5;
syms a
syms b
syms c
syms d
syms u
syms v

E_ab_u = SSEM(a,b,n,u);
E_cd_v = SSEM(c,d,n,v); 
e_ab_u = SCM(n,eye(n),{E_ab_u});
e_cd_v = SCM(n,eye(n),{E_cd_v});
com = SCM.commutator(e_ab_u,e_cd_v);

% com
% com.Concrete_base
% length_before_simplifying = length(com.List_of_terms)
simplified_commutator = com.cancelPairs();
% length_after_simplifying = length(simplified_commutator.List_of_terms)

for t=1:length(simplified_commutator.List_of_terms)
    t;
    term_t = simplified_commutator.List_of_terms{t};
    entry_t = simplified_commutator.List_of_terms{t}.Entry;
    str_entry_t = simplified_commutator.List_of_terms{t}.Entry.strForm();
    sym_entry_t = simplified_commutator.List_of_terms{t}.Entry.symForm();
    modified_term_t = term_t;
    modified_term_t.Entry = sym_entry_t;
    modified_term_t
end
