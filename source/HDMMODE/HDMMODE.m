classdef HDMMODE < ALGORITHM
    % <multi> <real/binary/permutation> <multimodal>
    % Coevolutionary constrained multi-objective optimization framework
    % type --- 1 --- Type of operator (1. GA 2. DE)
    %%
    %------------------------------- Reference --------------------------------
    % Y. Tian, T. Zhang, J. Xiao, X. Zhang, and Y. Jin, A coevolutionary
    % framework for constrained multi-objective optimization problems, IEEE
    % Transactions on Evolutionary Computation, 2020.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------



    methods
        function main(Algorithm,Problem)

            delta=3;
            fk=1;
            change_rate = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);

            Population = Problem.Initialization();
            gen=1;
            K = 2;
            last_gen = 100;
            change_threshold = 0.01;
            Flag = 0;
            NND = floor(Problem.N/K);

            Populations = cell(1,K);
            index       = randperm(floor(Problem.N/K)*K);
            temp        = reshape(index,K,floor(Problem.N/K));
            MaxGen = Problem.maxFE./ Problem.N;

            for i = 1 : K
                Populations{i} = Population(temp(i,:));
            end
            Archive = cell(1,K);
            rank = 1:K;

            for i=1:K
                [Front_number,~]= NDSort(Populations{rank(i)}.objs,inf);
                inter_Archive = [Populations{rank(i)}(Front_number==1),Archive{rank(i)}];
                [a_Front_number,~]= NDSort(inter_Archive.objs,inf);
                Archive{rank(i)} = inter_Archive(a_Front_number==1);
                clear a_Front_number Front_number
                

                while length(Archive{rank(i)})>1*length(Populations{rank(i)})
                    Del  = Truncation2( Archive{rank(i)}.objs, length(Archive{rank(i)}) - 1*length(Populations{rank(i)}) );
                    Archive{rank(i)} = Archive{rank(i)}(Del');
                    clear Del
                end
                            
            end


            while Algorithm.NotTerminated(Population)
                for i = 1 : K
                    if i==1

                        change_rate = Normalization2(Population,change_rate,gen);
                        CrowdDis{rank(i)}=Crowding(Populations{rank(i)}.decs,fk);
                        MatingPool = TournamentSelection(2,2*length(Populations{rank(i)}),-CrowdDis{rank(i)});
                        Offspring  = OperatorDE(Problem, Populations{rank(i)},Populations{rank(i)}(MatingPool(1:end/2)),Populations{rank(i)}(MatingPool(end/2+1:end)));


                        if Flag == 1
                            Flag = 1;
                        else
                            if gen > last_gen && Convertion2(change_rate,gen,last_gen,change_threshold)
                                Flag = 1;
                            end
                        end

                        if Flag==1
                            if NND == floor(Problem.N/K)/2 
                               NND = floor(Problem.N/K)/2;
                            else
                                NND = NND-1;
                            end
                            Populations{rank(i)} = [Populations{rank(i)},Offspring];
                            [Populations{rank(i)}] = Environmental_Selection(Populations{rank(i)},NND,Problem.upper,Problem.lower,delta,fk);
                            [Front_number,~]= NDSort(Populations{rank(i)}.objs,inf);
                            inter_Archive = [Populations{rank(i)}(Front_number==1),Archive{rank(i)}];
                            [a_Front_number,~]= NDSort(inter_Archive.objs,inf);
                            Archive{rank(i)} = inter_Archive(a_Front_number==1);
                            clear a_Front_number Front_number Offspring
                            

                            while length(Archive{rank(i)})>1*length(Populations{rank(i)})
                                Del  = Truncation2( Archive{rank(i)}.objs, length(Archive{rank(i)}) - 1*length(Populations{rank(i)}) );
                                Archive{rank(i)} = Archive{rank(i)}(Del');
                                clear Del
                            end
                            
                        else

                            Populations{rank(i)} = [Populations{rank(i)},Offspring];
                            [Populations{rank(i)}] = Environmental_Selection(Populations{rank(i)},floor(Problem.N/K),Problem.upper,Problem.lower,delta,fk);
                            [Front_number,~]= NDSort(Populations{rank(i)}.objs,inf);
                            inter_Archive = [Populations{rank(i)}(Front_number==1),Archive{rank(i)}];
                            [a_Front_number,~]= NDSort(inter_Archive.objs,inf);
                            Archive{rank(i)} = inter_Archive(a_Front_number==1);
                            clear a_Front_number Front_number Offspring
                            

                            while length(Archive{rank(i)})>1*length(Populations{rank(i)})
                                Del  = Truncation2( Archive{rank(i)}.objs, length(Archive{rank(i)}) - 1*length(Populations{rank(i)}) );
                                Archive{rank(i)} = Archive{rank(i)}(Del');
                                clear Del
                            end
                            
                        end

                    else
                        
                        CrowdDis{rank(i)}=Crowding(Populations{rank(i)}.decs,fk);
                        MatingPool = TournamentSelection(2,2*length(Populations{rank(i)}),-CrowdDis{rank(i)});
                        Offspring  = OperatorDE(Problem,Populations{rank(i)},Populations{rank(i)}(MatingPool(1:end/2)),Populations{rank(i)}(MatingPool(end/2+1:end)));
                        Populations{rank(i)} = [Populations{rank(i)},Offspring];

                    
                        [All_front,~] = NDSort(Archive{rank(1:i-1)}.objs,inf);
                        Archive_for_Compare = Archive{rank(1:i-1)}(All_front==1);

                        Range=max([Populations{rank(1)}.decs;Populations{rank(K)}.decs],[],1)-...
                            min([Populations{rank(1)}.decs;Populations{rank(K)}.decs],[],1); 

                        rr=(Range).*(0.2); 
                        Penalty_index = Penalty_Choose(Populations{rank(i)},Archive_for_Compare,rr); 

                        [Populations{rank(i)}] = Environmental_Selection_2(Populations{rank(i)},(2.*floor(Problem.N/K) - NND),...
                            Problem.upper,Problem.lower,delta,fk,Penalty_index);
                      
                            
                    end
                end
                Population = [Populations{:}];
                gen=gen+1;

            end
        end
    end

end
