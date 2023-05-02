 
% % Type B
% for q=2:8
%     fprintf(['B',num2str(q),'\n'])
%     B_q = RootSystem('B',q);
%     for i=1:length(B_q.RootList)
%         alpha = B_q.RootList{i};
%         for j=1:length(B_q.RootList)
%             beta = B_q.RootList{j};
%             combos = B_q.LinearCombos(alpha,beta);
%             num_combos = length(combos);
%             
%             if num_combos > 0
%                 if dot(alpha,alpha)==1 && dot(beta,beta)==1
%                     assert(num_combos==1)
%                 elseif dot(alpha,alpha)==1 && dot(beta,beta)==2
%                     assert(num_combos==2)
%                 elseif dot(alpha,alpha)==2 && dot(beta,beta)==1
%                     assert(num_combos==2)
%                 elseif dot(alpha,alpha)==2 && dot(beta,beta)==2
%                     assert(num_combos==1)
%                 end
%             end
%         end
%     end
% end

% % Type C
% for q=2:8
%     fprintf(['C',num2str(q),'\n'])
%     C_q = RootSystem('C',q);
%     for i=1:length(C_q.RootList)
%         alpha = C_q.RootList{i};
%         for j=1:length(C_q.RootList)
%             beta = C_q.RootList{j};
%             combos = C_q.LinearCombos(alpha,beta);
%             num_combos = length(combos);
%             
%             if dot(alpha,alpha)==4 && dot(beta,beta)==4
%                 assert(num_combos==0)
%             end
% 
%             if num_combos > 0
%                 if dot(alpha,alpha)==2 && dot(beta,beta)==2
%                     assert(num_combos==1)
%                 elseif dot(alpha,alpha)==2 && dot(beta,beta)==4
%                     assert(num_combos==2)
%                 elseif dot(alpha,alpha)==4 && dot(beta,beta)==2
%                     assert(num_combos==2)
%                 end
%             end
%         end
%     end
% end

% Type BC
for q=2:8
    fprintf(['BC',num2str(q),'\n'])
    BC_q = RootSystem('BC',q);
    for i=1:length(BC_q.RootList)
        alpha = BC_q.RootList{i};
        for j=1:length(BC_q.RootList)
            beta = BC_q.RootList{j};
            
            if not(RootSystem.IsProportionalRoot(alpha,beta))
                
                combos = BC_q.LinearCombos(alpha,beta);
                num_combos = length(combos);
                
    
                if dot(alpha,alpha)==4 && dot(beta,beta)==4
                    assert(num_combos==0)
                end
    
                if num_combos > 0
                    if dot(alpha,alpha)==1 && dot(beta,beta)==1
                        assert(num_combos==1)
                    elseif dot(alpha,alpha)==1 && dot(beta,beta)==2
                        assert(num_combos==3)
                    elseif dot(alpha,alpha)==1 && dot(beta,beta)==4
                        assert(num_combos==0)
                    elseif dot(alpha,alpha)==2 && dot(beta,beta)==1
                        assert(num_combos==3)
                    elseif dot(alpha,alpha)==2 && dot(beta,beta)==2
                        assert(num_combos==1)
                    elseif dot(alpha,alpha)==2 && dot(beta,beta)==4
                        assert(num_combos==2)
                    elseif dot(alpha,alpha)==4 && dot(beta,beta)==1
                        assert(num_combos==0)
                    elseif dot(alpha,alpha)==4 && dot(beta,beta)==2
                        assert(num_combos==2)
                    elseif dot(alpha,alpha)==4 && dot(beta,beta)==4
                        assert(num_combos==0)
                    end
                end
            end
        end
    end
end