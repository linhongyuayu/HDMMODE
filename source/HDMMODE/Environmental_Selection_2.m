function Population=Environmental_Selection_2(Pop,N,xu,xl,delta,fk,Penalty_index)
PopDec=Pop.decs;
PopObj=Pop.objs;


PopObj(Penalty_index,:) = repmat(max(PopObj,[],1),length(Penalty_index),1);
NNN = size(PopDec,1);
Normalize_PopDec = (PopDec- repmat(xl,NNN,1))./repmat(xu-xl,NNN,1);
label_3=false;
range=max(PopDec,[],1)-min(PopDec,[],1);
r=(range)*0.2;
C=NCM(PopDec,r);
K=length(C);
local_A=[];

for i=1:K
    cluster=C(i).p;
    if length(cluster)<=delta
        continue;
    end

    [FrontNo1,~]=NDSort(Pop(cluster).objs,length(cluster));
    Inter_index = cluster(FrontNo1==1);
    Inter_Pop = Pop(Inter_index);

    if length(Inter_Pop)~=1
        Normalize_inter_PopDec = Normalize_PopDec(Inter_index',:);

        inter_PopObj = PopObj(Inter_index,:);
        inter_Distance = pdist2(inter_PopObj,inter_PopObj);
        inter_Distance(logical(eye(length(inter_Distance)))) = inf;
        inter_Distance = sort(inter_Distance,2);
        inter_CrowdDis_obj = (inter_Distance(:,floor(sqrt(size(inter_PopObj,1))))+2);
        Max_inter_CrowdDis_obj = max(inter_CrowdDis_obj);
        Min_inter_CrowdDis_obj = min(inter_CrowdDis_obj);
        Nor_inter_CrowdDis_obj = (inter_CrowdDis_obj - Min_inter_CrowdDis_obj)./(Max_inter_CrowdDis_obj - Min_inter_CrowdDis_obj);

        for j=1:length(Inter_Pop)
            Self_inter_dec = Normalize_inter_PopDec(j,:);
            Niche_inter_Dec = Niche_induce(Self_inter_dec,Normalize_PopDec,fk);
            distance_inter_Dec = sum(  abs((ones(size(Niche_inter_Dec,1),1)*Self_inter_dec-Niche_inter_Dec)).^(fk),2).^(1/fk) ;
            Mean_inter_distance_Dec(j,:) = sum(distance_inter_Dec,1)./length(distance_inter_Dec);
            Min_inter_distance_Dec(j,:) = min(distance_inter_Dec);
            inter_MEDX(j,:) = Mean_inter_distance_Dec(j).*Min_inter_distance_Dec(j);
            clear Niche_inter_Dec distance_inter_Dec
        end
        Max_inter_MEDX = max(inter_MEDX);
        Min_inter_MEDX = min(inter_MEDX);
        Nor_inter_MEDX = (inter_MEDX - Min_inter_MEDX)./(Max_inter_MEDX - Min_inter_MEDX);

        inter_F_score = 2./(1./(Nor_inter_CrowdDis_obj) + 1./Nor_inter_MEDX);
        [~,inter_F_index] = sort(inter_F_score,1,'descend');
        inter_min_index_local = Inter_index(inter_F_index(1));
        local_A=[local_A inter_min_index_local];
    else
        local_A=[local_A cluster(FrontNo1==1)];
    end
    clear FrontNo1 inter_min_index_local inter_F_score  nter_F_index Min_inter_distance_Dec Mean_inter_distance_Dec
    clear inter_CrowdDis_obj Inter_index Inter_Pop Normalize_inter_PopDec inter_Distance inter_PopObj inter_MEDX inter_F_index
    clear Max_inter_CrowdDis_obj Min_inter_CrowdDis_obj Nor_inter_CrowdDis_obj Max_inter_MEDX Min_inter_MEDX Nor_inter_MEDX
end


if length(local_A)>N
    P = local_A;
else
    label_3=true;
    [FrontNo,~]=Ndsort(PopObj,length(Pop));
    P = find(FrontNo==1);
    temp=setdiff(P,local_A);
    P=[local_A temp];
    Front_number=1;
end

%%
while  length(P)<=N
    Front_number=Front_number+1;
    temp1=find(FrontNo==Front_number);
    temp2=setdiff(temp1, local_A);
    P=[P temp2];
end
%%
if label_3==true && Front_number==1
    temp2=temp;
end
interPop = setdiff(P, temp2);


if length(P)>N
    Add_Number=N-length(interPop);
    Last_index = temp2';
    Last_PopDec = Normalize_PopDec(Last_index,:);
    Last_PopObj = PopObj(Last_index,:);
    Last_number = size(Last_PopDec,1);
    Last_Distance = pdist2(Last_PopObj,Last_PopObj);
    Last_Distance(logical(eye(length(Last_Distance)))) = inf;
    Last_Distance = sort(Last_Distance,2);
    Last_CrowdDis_obj = (Last_Distance(:,floor(sqrt(size(Last_PopObj,1))))+2);
    Max_Last_CrowdDis_obj = max(Last_CrowdDis_obj);
    Min_Last_CrowdDis_obj = min(Last_CrowdDis_obj);
    Nor_Last_CrowdDis_obj = (Last_CrowdDis_obj - Min_Last_CrowdDis_obj)./(Max_Last_CrowdDis_obj - Min_Last_CrowdDis_obj);

    for i=1:Last_number
        Self_dec = Last_PopDec(i,:);
        Niche_Dec = Niche_induce(Self_dec,Normalize_PopDec(P',:),fk);
        distance_Dec = sum( abs((ones(size(Niche_Dec,1),1)*Self_dec-Niche_Dec)).^(fk),2).^(1/fk);
        Mean_distance_Dec(i,:) = sum(distance_Dec,1)./length(distance_Dec);
        Min_distance_Dec(i,:) = min(distance_Dec);
        Last_MEDX(i,:) = Mean_distance_Dec(i).*Min_distance_Dec(i);
        clear Niche_Dec distance_Dec
    end
    Max_Last_MEDX = max(Last_MEDX);
    Min_Last_MEDX = min(Last_MEDX);
    Nor_Last_MEDX = (Last_MEDX - Min_Last_MEDX)./(Max_Last_MEDX - Min_Last_MEDX);

    F_score = 2./(1./(Nor_Last_CrowdDis_obj) + 1./Nor_Last_MEDX);
    [~,F_index] = sort(F_score,1,'descend');
    Add_index = Last_index(F_index(1:Add_Number));
end
NewP = [interPop Add_index'];
Population = Pop(NewP);

clear Mean_distance_Dec Min_distance_Dec MEDX CrowdDis_obj F_score F_index Add_index
clear Last_PopObj Last_PopDec Last_index FrontNo temp temp1 temp2 interPop
clear Max_Last_MEDX Min_Last_MEDX Nor_Last_MEDX
end

function Niche_Dec = Niche_induce(Self_dec,PopDec,fk)
nn=4;
Nn=size(PopDec,1);
[~,Neighbour_k]=sort( sum( abs((ones(Nn,1)*Self_dec-PopDec)).^(fk),2).^(1/fk) );
Go=false;
while Go==false
    Neighbour_dec = PopDec(Neighbour_k(2:nn),:);
    n=nn-1;
    for j=1:size(Neighbour_dec,2)
        Mean_dec(:,j) = sum(Neighbour_dec(:,j),1)/n;
        std_dec(:,j) = sqrt( sum( (Neighbour_dec(:,j) - repmat(Mean_dec(:,j),n,1)).^2  ,1)/n );
        range_l(:,j) = Mean_dec(:,j) - 3.*std_dec(:,j);
        range_u(:,j) = Mean_dec(:,j) + 3.*std_dec(:,j);
    end
    if nn<Nn
        Neighbour_next_dec = PopDec(Neighbour_k(nn+1),:);
        for jj=1:size(Neighbour_dec,2)
            index(:,jj) = Neighbour_next_dec(:,jj)>=range_l(:,jj) && Neighbour_next_dec(:,jj)<=range_u(:,jj);
        end
        Sum_index = sum(index,2);
        if Sum_index~=size(Neighbour_dec,2)
            Go=true;
        else
            nn=nn+1;
            Go=false;
        end
        clear Neighbour_dec
    else
        Go =true;
    end
    Niche_Dec = PopDec(Neighbour_k(2:nn),:);
end
end