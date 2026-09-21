classdef NMMF4 < PROBLEM
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
        num_of_peaks
        alpha
    end
    methods

        

        function Setting(obj)
            obj.M = 2;
            obj.dd=7;
            d=obj.M+obj.dd;
            if isempty(obj.D); obj.D = 10; end
            obj.lower(:,1:obj.M)    = zeros(1,obj.M) + 1e-10;
            obj.lower(:,obj.M+1 : obj.M+ obj.dd)    = -10.*ones(1,d-obj.M) + 1e-10;
            obj.lower(:,d+1:obj.D)= zeros(1,obj.D-d);

            obj.upper    = ones(1,obj.M) - 1e-10;
            obj.upper(:,obj.M+1:d)= 10.*ones(1,d-obj.M);
            obj.upper(:,d+1:obj.D)= ones(1,obj.D-d);
            obj.encoding = 'real';
            obj.num_of_peaks=8;
            obj.alpha = 0.96;
            
        end

        function PopObj = CalObj(obj,PopDec)
            OptX = 0.2;
            [N,~]  = size(PopDec);
            Pop=PopDec;
            M      = obj.M;


            obj.THETA_=zeros(N,1);
            for i=1:N
                obj.THETA_(i) = 2/pi*atan(Pop(i,2)./Pop(i,1));
                if Pop(i,1)==0
                    obj.THETA_(i) = 1;
                end
            end


            y=PopDec(:,obj.M + 1 : obj.M + obj.dd);
            T=cumsum(y,2);
            n=obj.M + 2 : obj.M + obj.dd;
            for k=1:N
                Q(k,:)=  100.*(1- sin( (1/8).*obj.num_of_peaks*pi*(T(k,1)-2)  )).^2 + 50.*sum(n.*(obj.alpha*T(k,2:end) + T(k,1:end-1)).^2) +...
                    20.*(1-(PopDec(k,1).^2 + PopDec(k,2).^2)).^2; 
            end


            obj.h = sum((PopDec(:,obj.M+obj.dd+1:obj.D) - OptX).^2,2);


            T_=zeros(N,1);
            G_=zeros(N,M);
            for i=1:N
                T_(i) = (1 - (Pop(i,1)^2 + Pop(i,2)^2)).^2  + Q(i) + obj.h(i) ;
                G_(i,1:M) = [ones(1,1) cumprod(sin(pi/2*obj.THETA_(i)),2)] .* [cos(pi/2*obj.THETA_(i)) ones(1,1)];
            end
            PopObj = G_ .* repmat((1+T_),1,M) ; 
        end


        function R = GetOptimum(~,~)
            R = load('NMMF4_v3_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end


        function R = GetPF(obj)
            if obj.M == 2
                R = load('NMMF4_v3_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF4_Reference_PSPF_data');
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