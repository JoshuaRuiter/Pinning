%%%%%%%%%%%%%%%%
%% THESE PASS %%
%%%%%%%%%%%%%%%%
% RunSLTests(3)
% RunSLTests(4)
% RunSLTests(5)
% RunSOTests(4,1)

% Root system axioms
% for n=2:10
%     A_n = RootSystem('A',n,n+1);
%     A_n.VerifyProperties();
% 
%     B_n = RootSystem('B',n,n);
%     B_n.VerifyProperties();
% 
%     C_n = RootSystem('C',n,n);
%     C_n.VerifyProperties();
%     
%     BC_n = RootSystem('BC',n,n);
%     BC_n.VerifyProperties();
% end

% G_2 = RootSystem('G',2,3);
% G_2.VerifyProperties();
% 
% F_4 = RootSystem('F',4,4);
% F_4.VerifyProperties();

%%%%%%%%%%%%%%%%%%%%%%%
%% THESE DO NOT PASS %%
%%%%%%%%%%%%%%%%%%%%%%%
% for n=4:10
%     D_n = RootSystem('D',n,n+2);
%     D_n.VerifyProperties();
% end
% 
% for n=6:8
%     E_n = RootSystem('E',n,8);
%     E_n.VerifyProperties();
% end

% This passes up until the Weyl group conjugation tests
RunSUTests(4,2)

% Goals to work towards
%RunSUTests(6,3)
%RunSUTests(6,2)
%RunSOTests(6,2)

