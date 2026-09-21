classdef NMMF12 < PROBLEM
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
        THETA_;
        K   
        q1;  
        q2   
        Q1   
        Q2   
        Q
    end
    methods

        function Setting(obj)
            obj.M = 2;
            obj.K=3;
            obj.q1=6;
            obj.q2=6;
            if isempty(obj.D); obj.D = 18; end

            obj.lower(:,1:obj.K)    = zeros(1,obj.K);
            obj.lower(:,obj.K+1 : obj.K+obj.q1) = -10.*ones(1,obj.q1);
            obj.lower(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2) = -10.*ones(1,obj.q2);

            obj.upper    = ones(1,obj.K);
            obj.upper(:,obj.K+1 : obj.K+obj.q1)= 10.*ones(1,obj.q1);
            obj.upper(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2)= 10.*ones(1,obj.q2);

            obj.lower(:,obj.K+obj.q1+obj.q2+1:obj.D)=zeros(1,obj.D-obj.K-obj.q1-obj.q2);
            obj.upper(:,obj.K+obj.q1+obj.q2+1:obj.D)=ones(1,obj.D-obj.K-obj.q1-obj.q2);
            obj.encoding = 'real';
        end

        function PopObj = CalObj(obj,PopDec)
            OptX = 0.2;
            [N,~]  = size(PopDec);
            M      = obj.M;
            Pop=PopDec;

            T(:,1)=PopDec(:,3);
            T(:,2)=sqrt( PopDec(:,1).^2 + PopDec(:,2).^2 );
            obj.THETA_=zeros(N,1);
            for i=1:N
                obj.THETA_(i) = 2/pi*atan(T(i,2)./T(i,1));
                if T(i,1)==0
                    obj.THETA_(i) = 1;
                end
            end

            for i=obj.K+1 : 2 : obj.K+obj.q1
                PopDec(:,i)=(1+cos(pi./2.*(i./(obj.K+obj.q1)))).*((1.8).*Pop(:,i)-obj.lower(i))...
                    -obj.THETA_.*(obj.upper(i)-obj.lower(i));
            end

            for i=obj.K+2 : 2 : obj.K+obj.q1
                PopDec(:,i)=(1+i./(obj.K+obj.q1)).*((1.8).*Pop(:,i)-obj.lower(i))...
                    -Pop(:,1).*(obj.upper(i)-obj.lower(i));
            end

            for i=obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2
                PopDec(:,i)=(1+cos(pi./2.*(i./(obj.K+obj.q1+obj.q2)))).*((1.8).*Pop(:,i)-obj.lower(i))...
                    -( (Pop(:,obj.K+1)-obj.lower(obj.K+1))./(obj.upper(obj.K+1)-obj.lower(obj.K+1))).*(obj.upper(i)-obj.lower(i));
            end

            X_Q1=PopDec(:,obj.K+1:obj.K+obj.q1);
            obj.Q1=Pinter(X_Q1,obj);
            X_Q2=PopDec(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2);
            obj.Q2=Pathological(X_Q2);
            obj.Q=obj.Q1+obj.Q2;
            obj.h = sum((PopDec(:,obj.K + obj.q1 + obj.q2 +1 : obj.D) - OptX).^2,2);

            T_=zeros(N,1);
            G_=zeros(N,M);
            for i=1:N
                T_(i) = (1 - (T(i,1)^2 + T(i,2)^2)).^2  + obj.Q(i) + obj.h(i);
                G_(i,1:M) = [ones(1,1) cumprod(obj.THETA_(i),2)] .* [1-obj.THETA_(i) ones(1,1)];
            end
            PopObj = G_ .* repmat((1+T_),1,M) ;
        end


        function R = GetOptimum(~,~)
            R = load('NMMF12_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                R = load('NMMF12_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF12_Reference_PSPF_data');
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

function y=Pinter(x,obj)
[N1,q1]=size(x);
y=zeros(N1,1);

for i=2:q1
    j=obj.K+i;
    if i==q1
        x(:,i+1)=x(:,1);
    end

    A=10.*x(:,i-1).*sin(x(:,i)-x(:,i-1)) + 5.*sin(x(:,i+1)-x(:,i));
    B= (x(:,i)-x(:,i-1)).^2 - j.*x(:,i) + j.*(x(:,i+1)) -5.*cos(x(:,i+1)-x(:,i)) +5;

    y=y+ 10.*(x(:,1)-1).^2 + j.*(x(:,i)-x(:,i-1)).^2 + 20.*j.*sin(A).^2 + j.*log10(1+i.*B.^2);
end
end

function y=Pathological(x)
[N1,q1]=size(x);
y=zeros(N1,1);
r=r_cal(x(:,1));
for i=2:q1
    if i==q1
        x(:,i+1)=x(:,1);
    end
    y1= 50.*sin(sqrt(100.* (x(:,i+1)-x(:,i)).^2 + (x(:,i)-x(:,i-1)).^2)).^2 + 10*r - 5;
    y2= 1+1.*( (x(:,i+1)-x(:,i)).^2 - 2.*(x(:,i+1)-x(:,i)).*(x(:,i)-x(:,i-1)) + (x(:,i)-x(:,i-1)).^2).^2;

    y=y+ 5+ y1./y2;
end
end

function r=r_cal(x)
                n=size(x,1);
                for jj=1:n
                    if x(jj,1)>5
                        r(jj,:)=25.*(x(jj,1)-6).^2;
                    elseif x(jj,1)<-5
                        r(jj,:)=25./4.*(x(jj,1)+7).^2;
                    else
                        r(jj,:)=(24/25).*x(jj,1).^2+1;
                    end
                end
            end