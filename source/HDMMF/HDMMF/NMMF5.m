classdef NMMF5 < PROBLEM
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
        THETA_
    end
    methods




        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 16; end
            obj.lower    = zeros(1,obj.M) + 1e-10;
            obj.lower(:,obj.M+1:obj.D)= -15.*ones(1,obj.D-obj.M);
            obj.upper    = 15.*ones(1,obj.M) - 1e-10;
            obj.upper(:,obj.M+1:obj.D)= 15.*ones(1,obj.D-obj.M);
            obj.encoding = 'real';
        end

        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            Pop=PopDec;
            M      = obj.M;
            OptX=0.2;

            obj.THETA_=zeros(N,1);
            for i=1:N
                obj.THETA_(i) = 2/pi*atan(Pop(i,2)./Pop(i,1));
                if Pop(i,1)==0
                    obj.THETA_(i) = 1;
                end
            end

            dd=2;
            x1=PopDec(:,M+1:M+dd);
            x2=PopDec(:,M+dd+1:M+2*dd);
            x3=PopDec(:,M+2*dd+1:M+3*dd);
            x4=PopDec(:,M+3*dd+1:M+4*dd);
            x5=PopDec(:,M+4*dd+1:M+5*dd);
            x6=PopDec(:,M+5*dd+1:M+6*dd);

            r1=r_cal(x1);
            r2=r_cal(x2);
            for i=1:N
                t1(i,:)=100.*((x1(i,2)-3)-0.01.*x1(i,1).^2+1).^2 + 0.01.*r1 + 1;
                t2(i,:)=100.*((x2(i,2)-3)-0.01.*x2(i,1).^2+1).^2 + 0.01.*r2 + 10;
                t3(i,:)=100.*((x3(i,2)-3)-0.01.*x3(i,1).^2+1).^2 + 0.01.*(x3(i,1)+10).^2  + 1;
                t4(i,:)=100.*((x4(i,2)-3)-0.01.*x4(i,1).^2+1).^2 + 0.01.*(x4(i,1)+10).^2  + 5;
                t5(i,:)=100.*((x5(i,2)-3)-0.01.*x5(i,1).^2+1).^2 + 0.01.*(x5(i,1)+10).^2  + 4;
                t6(i,:)=100.*((x6(i,2)-3)-0.01.*x6(i,1).^2+1).^2 + 0.01.*(x6(i,1)+10).^2  + 3;
            end
            for i=1:13
                t(i)=0.1*i;
                y(i)=exp(-t(i)) - 5.*exp(10.*t(i)) + 3.*exp(-4.*t(i));
                q(:,i)= ( t3.* exp(-t(i).* t1) - t4.* exp(t(i).* t2) +  t6.* exp(-t(i).* t5) - y(i)).^2;
            end
            Q=sum(q(:,1:end),2);

            obj.h = sum((PopDec(:,obj.M+6*dd+1:D) - OptX).^2,2);

            T_=zeros(N,1);
            G_=zeros(N,M);
            for i=1:N
                T_(i) = (1 - Pop(i,1)^2-Pop(i,2)^2).^2  + Q(i) + obj.h(i);
                G_(i,1:M) = [ones(1,1) obj.THETA_(i)] .* [1-obj.THETA_(i) ones(1,1)];
            end
            PopObj = G_ .* repmat((1+T_),1,M) ;
        end



        function R = GetOptimum(~,~)
            R = load('NMMF5_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                 R = load('NMMF5_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF5_Reference_PSPF_data');
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
   
            function r=r_cal(x)
                n=size(x,1);
                for jj=1:n
                    if x(jj,1)>5
                        r=25.*(x(jj,1)-6).^2;
                    elseif x(jj,1)<-5
                        r=(x(jj,1)+10).^2;
                    else
                        r=(24/25).*x(jj,1).^2+1;
                    end
                end
             end