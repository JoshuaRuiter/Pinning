%%%%%%%%%%%%%%%%
%% THESE PASS %%
%%%%%%%%%%%%%%%%

%%% Group pinning tests
% RunSLTests(3)
% RunSLTests(4)
% RunSLTests(5)
% RunSOTests(4,1)
% RunSUTests(4,2)

%%% Root system axioms
% for n=2:10
%     A_n = RootSystem('A',n);
%     A_n.VerifyProperties();
% 
%     B_n = RootSystem('B',n);
%     B_n.VerifyProperties();
% 
%     C_n = RootSystem('C',n);
%     C_n.VerifyProperties();
%     
%     BC_n = RootSystem('BC',n);
%     BC_n.VerifyProperties();
% 
%     if n >= 4
%         D_n = RootSystem('D',n);
%         D_n.VerifyProperties();
%     end
% end
% 
% for n=6:8
%     E_n = RootSystem('E',n,8);
%     E_n.VerifyProperties();
% end
% 
% G_2 = RootSystem('G',2,3);
% G_2.VerifyProperties();
% 
% F_4 = RootSystem('F',4,4);
% F_4.VerifyProperties();

%%%%%%%%%%%%%%%%%%%%%%%
%% THESE DO NOT PASS %%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Everything passes except Weyl group conjugation coefficients
%%% Cases where the Weyl group conjugation formula passes tests:
%%%   - two long roots
%%% Remaining cases:
%%%   - two medium roots
%%%   - beta long, alpha medium
%%%   - beta medium, alpha long
% RunSUTests(6,3)
% RunSUTests(8,4)
% RunSUTests(10,5)

%%% Goals to work towards
% RunSOTests(6,2)
RunSUTests(5,2)
RunSUTests(6,2)