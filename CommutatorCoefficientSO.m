function coeff = CommutatorCoefficientSO(MatrixSize,Root_System,Form,alpha,beta,i,j,u,v)
    % Given two roots alpha and beta, output the coefficient
    % N_{ij}^{alpha beta} (u,v)
    % that arises in the commutator of X(alpha,u) and X(beta,v)
    
    assert(Root_System.VectorLength == MatrixSize)
    assert(Root_System.IsRoot(alpha))
    assert(Root_System.IsRoot(beta))

    if ~Root_System.IsRoot(alpha+beta)
        % If alpha+beta is not a root, then the coefficient is zero
        % which is the same as saying that the commutator of X(alpha,u)
        % and X(beta,v) is X(alpha+beta,0)=I
        coeff = 0;

    else
        n = MatrixSize;
        q = Root_System.Rank;
        diff = n - 2*q;

        %c = sym('c',[1 diff]); % OLD VERSION
        c = Form.AnisotropicPartVector; % FIXED VERSION

        if IsShort(alpha) && IsShort(beta)
            % Two short roots case
            % add some assertion about the lengths of u and v

        elseif IsShort(alpha) && IsLong(beta)
            % One short and one long root case
            % add some assertion about the lengths of u and v

        elseif IsLong(alpha) && IsShort(beta)
            % One short and one long root case
            % add some assertion about the lengths of u and v

            % In this case, just reverse the roles of alpha and beta,
            % then use the previous case just with an additional negative sign

        elseif IsLong(alpha) && IsLong(beta)
            % Two long roots case
            % add some assertion about the lengths of u and v

        else
            % This should be impossible
            assert(False,"Two root lengths in the type B root system have the wrong length.");
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD VERSION
% Translate/transfer these cases to the above set of if-elseif statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if (i==1 && j==1 && length(u) == diff && length(v) == diff && abs(sum(alpha))== 1 && sum(beta)==1)
%             % Placeholder
%             coeff = -compare(alpha,beta)*(sum(c*u*v));
%         
%         elseif (i==1 && j==1 && sum(alpha) == 1 && sum(beta) == 0 && length(u) == diff && length(v) == 1)
%             coeff = sum(-u*v);
%         elseif (i==1 && j==1 && sum(alpha) == 0 && sum(beta) == 1 && length(u) == 1 && length(v) == diff)
%             coeff = sum(u*v);
%         elseif (i==2 && j==1 && sum(alpha) == 1 && sum(beta) == 0 && length(u) == diff && length(v) == 1)
%             coeff = sum(-v^2*c/2*u);
%         elseif (i==1 && j==2 && sum(alpha) == 0 && sum(beta) == 1 && length(u) == 1 && length(v) == diff)
%             coeff = sum(v^2*c/2*u);
%         end
%         % [Xai, Xaj] = Xai+aj(N11) -- N11(u,v) =  e( c1u1v1 + ...
%         % e is a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

function d = compare(alpha,beta)
    d = find(alpha)<find(beta);
end
function bool = IsShort(alpha)
    bool = (dot(alpha,alpha)==1);
end
function bool = IsLong(alpha)
    bool = (dot(alpha,alpha)==2);
end