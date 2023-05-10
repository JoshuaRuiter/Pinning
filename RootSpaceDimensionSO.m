function dim = RootSpaceDimensionSO(MatrixSize,Root_System,alpha) %#ok<INUSD>
    switch abs(sum(alpha))
        case 0
            dim = MatrixSize - 2*RootSystemRank;
        case 1
            dim = 1;
        case 2
            dim = MatrixSize - 2*RootSystemRank;
    end
end