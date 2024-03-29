% This script runs a battery of tests on code describing root systems,
% special linear groups, special orthogonal groups, and special unitary groups.

% Tests for the Symbolic Kronecker Delta,
% Symbolic Single Entry Matrix, and
% Symbolic Coordinates Matrix classes
% SKD.runTests()
% SSEM.runTests()
% SCM.runTests()

% % Root system tests
% % All of these pass
% upperBound = 10; % Increase this to run tests on higher rank root systems
% for n=2:upperBound
%     A_n = RootSystem('A',n);
%     A_n.VerifyProperties();
%     B_n = RootSystem('B',n);
%     B_n.VerifyProperties();
%     C_n = RootSystem('C',n);
%     C_n.VerifyProperties();
%     BC_n = RootSystem('BC',n);
%     BC_n.VerifyProperties();
% end
% for n=4:upperBound
%     D_n = RootSystem('D',n);
%     D_n.VerifyProperties();
% end
% for n=6:8
%     E_n = RootSystem('E',n,8);
%     E_n.VerifyProperties();
% end
% G_2 = RootSystem('G',2,3);
% G_2.VerifyProperties();
% F_4 = RootSystem('F',4,4);
% F_4.VerifyProperties();

% % Special linear group tests
% % All of these pass
% RunSLTests(3)
% RunSLTests(4)
% RunSLTests(5)

% % % Special orthogonal group tests
% % % Quasisplit special orthogonal groups (n=2q+2)
RunSOTests(4,1)
RunSOTests(6,2)
RunSOTests(8,3)
RunSOTests(10,4)
% % % Non-quasisplit special orthogonal groups (n>2q+2)
RunSOTests(5,1)
RunSOTests(6,1)
RunSOTests(7,2)
RunSOTests(8,2)
RunSOTests(9,3)
RunSOTests(10,3)

% % Quasisplit special unitary group tests (n=2q)
% % These pass tests up to and including commutator coefficients
% RunSUTests(4,2,1)
% RunSUTests(4,2,-1)
% RunSUTests(6,3,1)
% RunSUTests(6,3,-1)
% RunSUTests(8,4,1)
% RunSUTests(8,4,-1)
% RunSUTests(10,5,1)
% RunSUTests(10,5,-1)

% Non-quasisplit special unitary groups (n>2q)
% % These pass tests before commutator coefficients, 
% % but do not pass all commutator coefficient tests yet.
% RunSUTests(5,2,1)
% RunSUTests(5,2,-1)
% RunSUTests(6,2,1)
% RunSUTests(6,2,1)
% RunSUTests(6,2,-1)
% RunSUTests(8,3,1)
% RunSUTests(8,3,-1)
% RunSUTests(9,3,1)
% RunSUTests(9,3,-1)


