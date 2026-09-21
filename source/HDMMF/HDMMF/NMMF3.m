classdef NMMF3 < PROBLEM
    % <multi> <real> <multimodal>
    % Multi-modal multi-objective test function
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    properties  
        Sx; 
        h;  
        THETA_;
        dd
    end
    methods



        function Setting(obj)
            if isempty(obj.D); obj.D = 10; end
            obj.M = 2;
            obj.dd=obj.D - obj.M - 3;
            d = obj.M+obj.dd;

            obj.lower(:,1)= -1 + 1e-10;
            obj.upper(:,1)= 1 - 1e-10;
            obj.lower(:,2:obj.M)    = zeros(1,obj.M-1) + 1e-10;
            obj.upper(:,2:obj.M)    = ones(1,obj.M-1) - 1e-10;
            obj.lower(:,obj.M+1:d)= -5.*ones(1,d-obj.M);
            obj.upper(:,obj.M+1:d)= 5.*ones(1,d-obj.M);
            obj.lower(:,d+1:obj.D)= zeros(1,obj.D-d);
            obj.upper(:,d+1:obj.D)= ones(1,obj.D-d);
            
            obj.encoding = 'real';
        end

            function PopObj = CalObj(obj,PopDec)
                OptX = 0.2;               
                [N,D]  = size(PopDec);
                Pop=PopDec;
                M      = obj.M;

                obj.THETA_=zeros(N,1);
                for i=1:N
                    if Pop(i,1)>0
                        obj.THETA_(i) = 2/pi*atan(abs(Pop(i,2))./Pop(i,1));
                    elseif Pop(i,1)==0
                        obj.THETA_(i) = 1;
                    elseif Pop(i,1)<0
                        obj.THETA_(i) = 2/pi*atan(abs(Pop(i,2))./(Pop(i,1)+ 1));
                    end
                end 

                t=cumsum(PopDec(:,obj.M+1 : obj.M+obj.dd),2);
                num_of_peaks=2;
                temp1=(sin( (1/4).*num_of_peaks*pi.*t(:,end))).^5;
                Q1 = 1-temp1;
                Q=50.*Q1;               
           
                obj.h = 20 - 20 * exp(-0.2 * sqrt(sum((PopDec(:,obj.M+obj.dd+1:end) - OptX).^2,2)/(D-obj.M-obj.dd))) + exp (1)...  
                            - exp(sum(cos(2 * pi .*(PopDec(:,obj.M+obj.dd+1:end) - OptX)),2)/(D-obj.M-obj.dd));  

                T_=zeros(N,1);
                G_=zeros(N,M);
                for i=1:N
                    if Pop(i,1)>=0
                        T_(i) = (0.97 - Pop(i,1) - Pop(i,2)).^2 + Q(i) + obj.h(i);
                        G_(i,1:M) = 1-[ones(1,1) cumprod(sin(pi/2*obj.THETA_(i)),2)] .* [cos(pi/2*obj.THETA_(i)) ones(1,1)];
                    elseif Pop(i,1)<0
                        T_(i) = (1/4 - (Pop(i,1)+1).^2-6*Pop(i,2).^2).^2 + + Q(i) + obj.h(i);
                        G_(i,1:M) = 1-[ones(1,1) cumprod(sin(pi/2*obj.THETA_(i)),2)] .* [cos(pi/2*obj.THETA_(i)) ones(1,1)];
                    end
                end
                PopObj = G_ .* repmat((1+T_),1,M) ;          
            end

        function R = GetOptimum(obj,N)
            R = load('NMMF3_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                 R = load('NMMF3_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF3_Reference_PSPF_data');
            obj.POS = PS;
            obj.optimum = PF;
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
     end
 end
