%%%%%%%%%%%%%%%%
%% THESE PASS %%
%%%%%%%%%%%%%%%%
RunSLTests(3)
RunSLTests(4)
RunSLTests(5)
RunSOTests(4,1)

% Root system axioms
for n=2:10
    A_n = RootSystem('A',n,n);
    A_n.VerifyProperties();

    B_n = RootSystem('B',n,n);
    B_n.VerifyProperties();

    C_n = RootSystem('C',n,n);
    C_n.VerifyProperties();
    
    BC_n = RootSystem('BC',n,n);
    BC_n.VerifyProperties();
end

%%%%%%%%%%%%%%%%%%%%%%%
%% THESE DO NOT PASS %%
%%%%%%%%%%%%%%%%%%%%%%%

% This passes up until the Weyl group conjugation tests
% Not sure what is wrong here
%RunSUTests(4,2) 

%RunSUTests(6,2)
%RunSOTests(6,2)

