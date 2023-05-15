% This script runs a battery of tests on code describing root systems,
% special linear groups, special orthogonal groups, and special unitary groups.

% Root system tests
% All of these pass
upperBound = 6; % Increase this to run tests on higher rank root systems
for n=2:upperBound
    A_n = RootSystem('A',n);
    A_n.VerifyProperties();
    B_n = RootSystem('B',n);
    B_n.VerifyProperties();
    C_n = RootSystem('C',n);
    C_n.VerifyProperties();
    BC_n = RootSystem('BC',n);
    BC_n.VerifyProperties();
end
for n=4:upperBound
    D_n = RootSystem('D',n);
    D_n.VerifyProperties();
end
for n=6:8
    E_n = RootSystem('E',n,8);
    E_n.VerifyProperties();
end
G_2 = RootSystem('G',2,3);
G_2.VerifyProperties();
F_4 = RootSystem('F',4,4);
F_4.VerifyProperties();

% Special linear group tests
% All of these pass
RunSLTests(3)
RunSLTests(4)
RunSLTests(5)

% Special orthogonal group tests
% Quasisplit special orthogonal groups (n=2q+2)
RunSOTests(4,1)
RunSOTests(6,2)
RunSOTestS(8,3)
RunSOTestS(10,4)
% Non-quasisplit special orthogonal groups (n>2q+2)
RunSoTests(7,2)
RunSOTests(8,2)
RunSOTests(9,3)
RunSOTests(10,3)

% Special unitary group tests
% Quasisplit special unitary groups (n=2q)
RunSUTests(4,2)
RunSUTests(6,3)
RunSUTests(8,4)
RunSUTests(10,5)
% Non-quasisplit special unitary groups (n>2q)
RunSUTests(5,2)
RunSUTests(6,2)
RunSUTests(8,3)
RunSUTests(9,3)


