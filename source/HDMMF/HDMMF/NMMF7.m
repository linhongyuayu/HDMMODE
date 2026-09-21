classdef NMMF7 < PROBLEM
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
    end
    methods

        function Setting(obj)
            obj.M=2;
            if isempty(obj.D); obj.D = 16; end
            obj.lower    = zeros(1,obj.M) + 1e-10;
            obj.lower(:,obj.M+1:obj.D)= -15.*ones(1,obj.D-obj.M);
            obj.upper    = 15.*ones(1,obj.M) - 1e-10;
            obj.upper(:,obj.M+1:obj.D) = 15.*ones(1,obj.D-obj.M);
            obj.encoding = 'real';
        end

        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            Pop = PopDec;
            M      = obj.M;
            dd=2;
            OptX=0.2;
            d = obj.M + 4*dd;

            obj.THETA_=zeros(N,1);
            for i=1:N
                obj.THETA_(i) = 2/pi*atan(Pop(i,2)./Pop(i,1));
                if Pop(i,1)==0
                    obj.THETA_(i) = 1;
                end
            end

            xx=PopDec(:,1:M);
            x1=PopDec(:,M+1:M+dd);
            x2=PopDec(:,M+dd+1:M+2*dd);
            x3=PopDec(:,M+2*dd+1:M+3*dd);
            x4=PopDec(:,M+3*dd+1:M+4*dd);

            r1=r_cal(x1);
            r2=r_cal(x2);
            for i=1:N
                t1(i,:)=100.*((x1(i,2)-3) - 0.01.*x1(i,1).^2 + 1).^2 + 0.01.*r1(i) + 20.*(1 - (xx(i,1).^2 + xx(i,2).^2)).^2;
                t2(i,:)=100.*((x2(i,2)-3) - 0.01.*x2(i,1).^2 + 1).^2 + 0.01.*r2(i) + 1 + 20.*(1- (xx(i,1).^2 + xx(i,2).^2)).^2;
                t3(i,:)=(100.*((x3(i,2)-3) - 0.01.*x3(i,1).^2 + 1).^2 + 0.01.*(x3(i,1)+10).^2) +1 + 20.*(1-(xx(i,1).^2 + xx(i,2).^2)).^2;
                t4(i,:)=(100.*((x4(i,2)-3) - 0.01.*x4(i,1).^2 + 1).^2 + 0.01.*(x4(i,1)+10).^2) +1 + 20.*(1- (xx(i,1).^2 + xx(i,2).^2)).^2;
            end
            Q=(exp(-1.*t1)-t2).^4 + 100.*(t2-t3).^6 + (tan( 5000.*(t3-t4))).^4 + t1.^8;

            obj.h = sum((PopDec(:,d+1:D) - OptX).^2,2);

            T_=zeros(N,1);
            G_=zeros(N,M);
            for i=1:N
                T_(i) = (1 - Pop(i,1)^2-Pop(i,2)^2).^2  + Q(i) + obj.h(i);
                A = 2; 
                G_ = 1-[ones(N,1) cumprod(sin(pi/2*obj.THETA_),2)] .* [cos(pi/2*obj.THETA_) ones(N,1)];% P2
                G_(:,1) = obj.THETA_(:,1) - cos(2*pi*A*obj.THETA_(:,1)+pi/2)/2/A/pi;
            end
            PopObj = G_ .* repmat((1+T_),1,M) ;
            
    
            function r=r_cal(x)
                n=size(x,1);
                for jj=1:n
                    if x(jj,1)>5
                        r(jj,:)=25.*(x(jj,1)-6).^2;
                    elseif x(jj,1)<-5
                        r(jj,:)=(x(jj,1)+10).^2;
                    else
                        r(jj,:)=(24/25).*x(jj,1).^2+1;
                    end
                end
            end
        end


        function R = GetOptimum(~,~)
            R = load('NMMF7_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                R = load('NMMF7_Reference_PSPF_data.mat','draw_pf');
                R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF7_Reference_PSPF_data');
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