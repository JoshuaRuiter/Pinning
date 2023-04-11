classdef RootSystem

    properties
        Type % string
        Rank % positive integer
        VectorLength % positive integer
        RootList % cell array
    end

    methods

        % Constructor
        function obj = RootSystem(Type,Rank,VectorLength)

            assert(Rank > 0)

            obj.Type = Type;
            obj.Rank = Rank;
            switch upper(Type)
                case 'A'
                    defaultVectorLength = Rank+1;
                case 'E'
                    defaultVectorLength = 8;
                case 'G'
                    defaultVectorLength = 3;
                otherwise
                    defaultVectorLength = Rank;
            end
            if nargin == 2
                VectorLength = defaultVectorLength;
            else
                assert(VectorLength >= defaultVectorLength)
            end
            obj.VectorLength = VectorLength;

            % Construct the list of roots
            switch upper(Type)

                case 'A'
                    defaultVectorLength = Rank + 1;
                    obj.RootList = RootSystem.BuildTypeA(Rank,VectorLength);
                    assert(length(obj.RootList)==Rank*(Rank+1))
                
                case 'B'
                    obj.RootList = RootSystem.BuildTypeB(Rank,VectorLength);
                    assert(length(obj.RootList)==2*Rank^2)

                case 'C'
                    obj.RootList = RootSystem.BuildTypeC(Rank,VectorLength);
                    assert(length(obj.RootList)==2*Rank^2)
                    
                case 'BC'
                    obj.RootList = RootSystem.BuildTypeBC(Rank,VectorLength);
                    assert(length(obj.RootList)==2*Rank*(Rank+1))
                
                case 'D'
                    obj.RootList = RootSystem.BuildTypeD(Rank,VectorLength);
                    assert(length(obj.RootList)==2*Rank*(Rank-1))

                case 'E'
                    defaultVectorLength = 8;
                    assert(5 < Rank & Rank < 9)
                    obj.RootList = RootSystem.BuildTypeE(Rank,VectorLength);
                    switch Rank
                        case 6
                            assert(length(obj.RootList)==72)
                        case 7
                            assert(length(obj.RootList)==126)
                        case 8
                            assert(length(obj.RootList)==240)
                        otherwise
                            assert(false,'For root system E, rank must be 6, 7, or 8')
                    end

                case 'F'
                    assert(Rank == 4)
                    obj.RootList = RootSystem.BuildTypeF(Rank,VectorLength);
                    assert(length(obj.RootList)==48)

                case 'G'
                    defaultVectorLength = 3;
                    assert(Rank == 2)
                    obj.RootList = RootSystem.BuildTypeG(Rank,VectorLength);
                    assert(length(obj.RootList)==12)

                otherwise
                    assert(false,'Invalid root system type. The valid types are A, B, C, BC, D, E, F, and G')

            end
        end

        function bool = IsRoot(obj,alpha)
            bool = any(cellfun(@isequal, obj.RootList, repmat({alpha}, size(obj.RootList))));
        end

        function combos = LinearCombos(obj,alpha,beta)
            % Return all positive integral linear combinations of alpha and
            % beta which are roots, along with data describing what
            % integers coefficients are used

            assert(obj.IsRoot(alpha))
            assert(obj.IsRoot(beta))

            if not(obj.IsRoot(alpha+beta))
                % If alpha + beta is not a root, we're done
                combos = {};
                return;
            else
                % If alpha + beta is a root, add it to the list of combos
                % also tracking the data of (i,j) where the newly added
                % combo is i*alpha+j*beta
                combos = {{alpha+beta,1,1}};
            end

            % Run a loop where each iteration, we try adding alpha and beta
            % to each existing combo
            while true
                newCombos = obj.IncrementCombos(alpha,beta,combos);
                if length(combos) == length(newCombos)
                    % If at any point, no new combos are added, then we are done
                    break;
                end
                combos = newCombos;
            end
        end
        function newCombos = IncrementCombos(obj,alpha,beta,oldCombos)
            newCombos = oldCombos;
            for k=1:length(oldCombos)
                root = oldCombos{k}{1};
                i = oldCombos{k}{2};
                j = oldCombos{k}{3};

                % If root+alpha is a root, then add it to the list of
                % linear combinations
                if obj.IsRoot(root+alpha)
                    % Determine if root+alpha is new
                    root_plus_alpha_new = true;
                    for p=1:length(oldCombos)
                        old_root = oldCombos{p}{1};
                        if isequal(root+alpha,old_root)
                            root_plus_alpha_new = false;
                            break;
                        end
                    end

                    if root_plus_alpha_new
                        % root+alpha is a new root, add it to the list
                        % along with the (i,j) data
                        newCombos{end+1} = {root+alpha,i+1,j};
                    end
                end

                % If root+beta is a root, then add it to the list of
                % linear combinations
                if obj.IsRoot(root+beta)
                    % Determine if root+alpha is new
                    root_plus_beta_new = true;
                    for p=1:length(oldCombos)
                        old_root = oldCombos{p}{1};
                        if isequal(root+beta,old_root)
                            root_plus_beta_new = false;
                            break;
                        end
                    end

                    if root_plus_beta_new
                        % root+beta is a new root, add it to the list
                        % along with the (i,j) data
                        newCombos{end+1} = {root+beta,i,j+1};
                    end
                end
            end
        end

        function VerifyProperties(obj)
            fprintf("Running tests to verify root system axioms for the type " ...
                + obj.Type + num2str(obj.Rank) + " root system.\n");

            % 0 is not a root
            assert(~obj.IsRoot(zeros(1,obj.VectorLength)));

            % Any reflection of any root should be a root
            fprintf("\tChecking that any reflection of a root is another root...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
                    reflected_root = RootSystem.ReflectRoot(alpha,beta);
                    assert(obj.IsRoot(reflected_root));
                end
            end
            fprintf("passed.\n");

            % Angle brackets are integers
            fprintf("\tChecking that angle brackets are integers...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};

                    % rem(a,b) is remainder of integer division a/b
                    % so rem(a,1)=0 if and only if x is an integer
                    angle_bracket = 2*dot(alpha,beta)/dot(beta,beta);
                    assert(rem(angle_bracket,1)==0);
                end
            end

            fprintf("passed.\n");
            % Only contains +/- and +/-2 multiples
            fprintf("\tChecking that the only multiples of a root are +/-1 and +/-2 (or +/-0.5)...");
            for i=1:length(obj.RootList)
                alpha = obj.RootList{i};
                for j=1:length(obj.RootList)
                    beta = obj.RootList{j};
                    
                    % rdivide does entry-wise division of matrices/vectors
                    % rdivide puts NaN for anything divide by zero
                    entrywise_division = rdivide(alpha,beta);
                    alpha_beta_are_proportional = (max(entrywise_division)== min(entrywise_division));
                    if alpha_beta_are_proportional
                        % alpha and beta are proportional

                        ratio = NaN;
                        for k=1:length(entrywise_division)
                            if ~isnan(entrywise_division(k))
                                ratio = entrywise_division(k);
                                break;
                            end
                        end
                        assert(~isnan(ratio));
                        assert(abs(ratio)==1 || abs(ratio)==2 || abs(ratio)==0.5);
                    end
                end
            end
            fprintf("passed.\n");

            fprintf("All tests passed.\n\n")
        end
    end

    methods (Static)

        function myCellArray = BuildTypeA(Rank,VectorLength)
            assert(VectorLength>Rank)
            myCellArray = cell(1, Rank*(Rank+1));
            index = 1;
            for i=1:Rank+1
                for j=1:Rank+1
                    if i ~= j 
                        x = zeros(1,VectorLength);
                        x(i) = 1;
                        x(j) = -1;
                        myCellArray{index} = x;
                        index = index + 1;
                    end
                end
            end
        end
        function myCellArray = BuildTypeB(Rank,VectorLength)
            % Build the type B root system
            assert(VectorLength >= Rank)
            q = Rank;
            myCellArray = cell(1,2*q^2);
            index = 1;
            
            for i=1:q

                % Make the short roots of the form +alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 1;
                index = index + 1;

                % Make the short roots of the form -alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = -1;
                index = index+ 1;

                % Make the long roots
                for j=i+1:q

                    % Make roots of the form alpha_i + alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i + alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form alpha_i - alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i - alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;
                end
            end
        end
        function myCellArray = BuildTypeC(Rank,VectorLength)
            % Build the type C root system
            assert(VectorLength >= Rank)
            q = Rank;
            myCellArray = cell(1,2*q^2);
            index = 1;
            
            for i=1:q

                % Make the long roots of the form 2*alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 2;
                index = index + 1;

                % Make the long roots of the form -2*alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = -2;
                index = index + 1;

                % Make the short length roots
                for j=i+1:q

                    % Make roots of the form alpha_i + alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i + alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form alpha_i - alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i - alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;
                end
            end
        end
        function myCellArray = BuildTypeBC(Rank,VectorLength)
            % Build the type BC root system
            assert(VectorLength >= Rank)
            q = Rank;
            myCellArray = cell(1,2*q*(q+1));
            index = 1;
            
            for i=1:q

                % Make the short roots of the form +alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 1;
                index = index + 1;

                % Make the short roots of the form -alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = -1;
                index = index+ 1;

                % Make the long roots of the form 2*alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 2;
                index = index + 1;

                % Make the long roots of the form -2*alpha_i
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = -2;
                index = index + 1;

                % Make the medium length roots
                for j=i+1:q

                    % Make roots of the form alpha_i + alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i + alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form alpha_i - alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i - alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;
                end
            end
        end
        function myCellArray = BuildTypeD(Rank,VectorLength) 
            % Technically, we could construct D3, but this is isomorphic to A3.
            assert(Rank > 3)
            assert(VectorLength >= Rank)
            q = Rank;
            myCellArray = cell(1,2*q*(q-1));
            % The roots have 2 entries that are +-1, either sign.
            index = 1;
            
            for i=1:q
                for j=i+1:q

                    % Make roots of the form alpha_i + alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i + alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form alpha_i - alpha_j
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 1;
                    myCellArray{index}(j) = -1;
                    index = index + 1;

                    % Make roots of the form -(alpha_i - alpha_j)
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -1;
                    myCellArray{index}(j) = 1;
                    index = index + 1;
                end
            end

        end
        function myCellArray = BuildTypeE(Rank,VectorLength)
            assert(5 < Rank && Rank < 9)
            assert(VectorLength > 7)
            q = Rank;
            rootCount = 2^(q-1) + 2*(q-1)*(q-2);
            if q==7
                rootCount = rootCount + 2;
            end
            myCellArray = cell(1,rootCount);

            % The roots are pairs of +-2 among the first q-1 entries
                % or flip an even number of the first q-1 +-1s
            % For E7, we also have (0,0,0,0,0,0,2,2) and the negative
            index = 1;
            % First, add the roots of 8 +-1s
            for i = 0:2^(q-1)-1
                signFlipCount = 0;
                myCellArray{index} = zeros(1,VectorLength);
                for j = 1:8
                    if mod(i,2^j) < 2^(j-1)
                        myCellArray{index}(j) = 1;
                        signFlipCount = signFlipCount + 1;
                    else
                        myCellArray{index}(j) = -1;
                    end
                end
                % If an odd number of bits were flipped, overwrite with the next root.
                if mod(signFlipCount, 2) == 0
                    index = index+1;

                    myCellArray{index} = zeros(1,VectorLength);
                    for j = 1:8
                        if mod(i,2^j) < 2^(j-1)
                            myCellArray{index}(j) = -1;
                        else
                            myCellArray{index}(j) = 1;
                        end
                    end
                    index = index+1;
                end
            end
            % Next, add roots of 2 +-2s
            % Here, we essentially construct D_(q-1)
            if q==8
                % Technically we are still working with E8, but the
                    % computation is simpler if we treat it as E9
                % We did not do this earlier because we would get duplicate
                    % roots of 8 +-1s
                q = 9;
            end
            for i=1:q-1
                for j=i+1:q-1
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 2;
                    myCellArray{index}(j) = 2;
                    index = index + 1;

                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -2;
                    myCellArray{index}(j) = -2;
                    index = index + 1;
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = 2;
                    myCellArray{index}(j) = -2;
                    index = index + 1;

                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index}(i) = -2;
                    myCellArray{index}(j) = 2;
                    index = index + 1;
                end
            end
            % Finally, add the 2 additional roots for E7
            if q==7
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(7) = 2;
                myCellArray{index}(8) = 2;
                index = index + 1;

                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(7) = -2;
                myCellArray{index}(8) = -2;
            end
        end
        function myCellArray = BuildTypeF(Rank,VectorLength)
            assert(Rank == 4)
            assert(VectorLength > 3)
            myCellArray = cell(1,48);
            % 1, 2, or 4 entries are +-1 (+-2 in the case of 1 entry)
            index = 1;
            for i = 1:4
                for j = 1:4
                    myCellArray{index} = zeros(1,VectorLength);
                    myCellArray{index+1} = zeros(1,VectorLength);
                    if i<j
                        myCellArray{index}(i) = 1;
                        myCellArray{index}(j) = 1;
                        myCellArray{index+1}(i) = -1;
                        myCellArray{index+1}(j) = -1;
                    elseif i>j
                        myCellArray{index}(i) = 1;
                        myCellArray{index}(j) = -1;
                        myCellArray{index+1}(i) = -1;
                        myCellArray{index+1}(j) = 1;
                    else
                        % Here, i==j
                        myCellArray{index}(i) = 2;
                        myCellArray{index+1}(j) = -2;
                    end
                    index = index + 2;
                end
            end
            for i = 0:15
                myCellArray{index} = zeros(1,VectorLength);
                if mod(i,2) == 0
                    myCellArray{index}(1) = 1;
                else
                    myCellArray{index}(1) = -1;
                end
                if mod(i,4) < 2
                    myCellArray{index}(2) = 1;
                else
                    myCellArray{index}(2) = -1;
                end
                if mod(i,8) < 4
                    myCellArray{index}(3) = 1;
                else
                    myCellArray{index}(3) = -1;
                end
                if i < 8
                    myCellArray{index}(4) = 1;
                else
                    myCellArray{index}(4) = -1;
                end
                index = index + 1;
            end

        end
        function myCellArray = BuildTypeG(Rank,VectorLength)
            assert(Rank == 2)
            assert(VectorLength > 2)
            myCellArray = cell(1,12);
            % 2 entries are +-1, while the third is 0 or +-2.
            % In addition, the entries must sum to 0.
            index = 1;
            for i = 1:3
                j = 1;
                k = 2;
                if i == 1
                    j = 3;
                elseif i == 2
                    k = 3;
                end
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 0;
                myCellArray{index}(j) = 1;
                myCellArray{index}(k) = -1;
                index = index + 1;
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 0;
                myCellArray{index}(j) = -1;
                myCellArray{index}(k) = 1;
                index = index + 1;
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = -2;
                myCellArray{index}(j) = 1;
                myCellArray{index}(k) = 1;
                index = index + 1;
                myCellArray{index} = zeros(1,VectorLength);
                myCellArray{index}(i) = 2;
                myCellArray{index}(j) = -1;
                myCellArray{index}(k) = -1;
                index = index + 1;
            end

        end

        function bool = IsProportionalRoot(alpha,beta)
            entrywise_division = rdivide(alpha,beta);
            bool = (max(entrywise_division)== min(entrywise_division));
        end

        function reflectedRoot = ReflectRoot(alpha, beta)
            % Given two roots alpha and beta, compute the reflection
            % of beta across the hyperplane perpendicular to alpha. 
            
            % In the notation of the notes document, compute sig_alpha(beta)
            reflectedRoot = beta - 2*dot(alpha,beta)/dot(alpha,alpha)*alpha;
        end
    end

end