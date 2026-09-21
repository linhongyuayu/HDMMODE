function [Penalty_index] = Penalty_Choose(Pop,Archive,rr)
    Pop_dec = Pop.decs;
    Archive_dec = Archive.decs;
    [n1,dim]=size(Pop_dec);
    [n2,~]=size(Archive_dec);
    D=zeros(n1,n2);
    for i=1:dim
        D=D+(abs(Pop_dec(:,i)-repmat(Archive_dec(:,i)',n1,1))<=rr(i));
    end
    [Penalty_index,~] = find(D==dim);

end