classdef NMMF2 < PROBLEM
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
        h;  
        dd
        num_of_peaks
    end
    methods


        %% 
        function Setting(obj)
            obj.M=2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower(:,1)    = 0.1;
            obj.upper(:,1)    = 1.1;
            obj.lower(:,2:obj.D)= -5.*ones(1,obj.D-1);
            obj.upper(:,2:obj.D)= 5.*ones(1,obj.D-1);
            obj.encoding = 'real';

            obj.num_of_peaks=2;
        end
        %% 
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            y=PopDec(:,2:D);
            T=cumsum(y,2);
            n=3:D;
            for k=1:N
                G(k,:)=  (1 + 10.*(1-sin( (1/8).*obj.num_of_peaks*pi*(T(k,1)-2))).^2 + 5.*sum(n.*(2*T(k,2:end) - T(k,1:end-1)).^2));
            end
            PopObj(:,1) = PopDec(:,1);
            PopObj(:,2) = G./PopDec(:,1);
        end

        %%
        function R = GetOptimum(obj,N)
            R = load('NMMF2_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end
        %% 
        function R = GetPF(obj)
            if obj.M == 2
                R = load('NMMF2_Reference_PSPF_data.mat','draw_pf');
                R=R.draw_pf;
            else
                R = [];
            end
        end
        %% 
        function score = CalMetric(obj,metName,Population)
            load('NMMF2_Reference_PSPF_data');
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