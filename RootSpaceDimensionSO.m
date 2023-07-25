function dim = RootSpaceDimensionSO(MatrixSize,Root_System,alpha)
    switch abs(sum(alpha))
        case 0
            dim = 1;
        case 1
            dim = MatrixSize - 2*Root_System.Rank;
        case 2
            dim = 1;
    end
end