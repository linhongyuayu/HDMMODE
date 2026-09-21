classdef NMMF11 < PROBLEM
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
        q3
        Q1   
        Q2   
        Q3
        Q
    end
    methods

        function Setting(obj)
            obj.M = 2;
            obj.K=3;
            obj.q1=5; 
            obj.q2=5;
            obj.q3=5;
            if isempty(obj.D); obj.D = 18; end

            obj.lower(:,1:obj.K)    = zeros(1,obj.K);
            obj.lower(:,obj.K+1 : obj.K+obj.q1) = -10.*ones(1,obj.q1);
            obj.lower(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2) = -10.*ones(1,obj.q2);
            
            obj.upper    = ones(1,obj.K);
            obj.upper(:,obj.K+1 : obj.K+obj.q1)= 10.*ones(1,obj.q1);
            obj.upper(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2)= pi.*ones(1,obj.q2);
            obj.lower(:,obj.K+obj.q1+obj.q2+1 : obj.K+obj.q1+obj.q2+obj.q3) = -10.*ones(1,obj.q3);
            obj.upper(:,obj.K+obj.q1+obj.q2+1 : obj.K+obj.q1+obj.q2+obj.q3) = 10.*ones(1,obj.q3);

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
    PopDec(:,i)=(1+cos(pi./2.*(i./(obj.K+obj.q1)))).*(2.*Pop(:,i)-obj.lower(i))...
        -Pop(:,1).*(obj.upper(i)-obj.lower(i));
end

for i=obj.K+2 : 2 : obj.K+obj.q1
    PopDec(:,i)=(1+i./(obj.K+obj.q1)).*(2.*Pop(:,i)-obj.lower(i))...
        -obj.THETA_.*(obj.upper(i)-obj.lower(i));
end


for i=obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2
    PopDec(:,i)=(1+cos(pi./2.*(i./(obj.K+obj.q1+obj.q2)))).*(2.*Pop(:,i)-obj.lower(i))...
        -( (Pop(:,obj.K+1)-obj.lower(obj.K+1))./(obj.upper(obj.K+1)-obj.lower(obj.K+1))).*(obj.upper(i)-obj.lower(i));
end


for i=obj.K+obj.q1+obj.q2+1 : 2 : obj.K+obj.q1+obj.q2+obj.q3
    PopDec(:,i)=(1+cos(pi./2.*(i./(obj.K+obj.q1+obj.q2+obj.q3)))).*(2.*Pop(:,i)-obj.lower(i))...
        -obj.THETA_.*(obj.upper(i)-obj.lower(i));
end

for i=obj.K+obj.q1+obj.q2+2 : 2 : obj.K+obj.q1+obj.q2+obj.q3
    PopDec(:,i)=(1+i./(obj.K+obj.q1+obj.q2+obj.q3)).*(2.*Pop(:,i)-obj.lower(i))...
        -Pop(:,1).*(obj.upper(i)-obj.lower(i));
end


X_Q1=PopDec(:,obj.K+1:obj.K+obj.q1);
obj.Q1=SineWave1(X_Q1,obj);


X_Q2=PopDec(:,obj.K+obj.q1+1 : obj.K+obj.q1+obj.q2);
obj.Q2=Trigonometric(X_Q2,obj);

X_Q3=PopDec(:,obj.K+obj.q1+obj.q2+1 : obj.K+obj.q1+obj.q2+obj.q3);
obj.Q3=SineWave2(X_Q3,obj);


obj.Q=obj.Q1 + obj.Q2 + obj.Q3; 

            
obj.h = sum((PopDec(:,obj.K + obj.q1 + obj.q2 + obj.q3 +1 : obj.D) - OptX).^2,2);            

                T_=zeros(N,1);
                G_=zeros(N,M);
                for i=1:N
                    T_(i) = (1 - (T(i,1)^2 + T(i,2)^2)).^2  + obj.Q(i) + obj.h(i);
                    A = 2; 
                    G_ = 1-[ones(N,1) cumprod(sin(pi/2*obj.THETA_),2)] .* [cos(pi/2*obj.THETA_) ones(N,1)];% P2
                    G_(:,1) = obj.THETA_(:,1) - cos(2*pi*A*obj.THETA_(:,1)+pi/2)/2/A/pi;
                end
                PopObj = G_ .* repmat((1+T_),1,M) ;          
            end


        function R = GetOptimum(~,~)
            R = load('NMMF11_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                R = load('NMMF11_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF11_Reference_PSPF_data');
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

function z=SineWave1(x,obj)
      [N1,q1]=size(x);
      y=zeros(N1,1);
      for i=1:q1-1
         j1=obj.K+i; j2=obj.K+q1;
          z=y + (( abs(x(:,i)) - j1./j2).^2 + 10.*(x(:,i+1)- ((j1+1)./j1).*x(:,i)).^2).^(1/2).*...
              (sin( (50.*((abs(x(:,i))-j1./j2).^2 + 10.*(x(:,i+1)- ((j1+1)./j1).*x(:,i)).^2).^(1/4))).^2+100) ;
      end
end

function z=SineWave2(x,obj)
      [N1,q1]=size(x);
      y=zeros(N1,1);
      for i=1:q1-1
          j1=obj.K+obj.q1+obj.q2+i; j2=obj.K+obj.q1+obj.q2+obj.q3;
          z=y + (( abs(x(:,i)) - j1./j2).^2 + 100.*(x(:,i+1)- ((j1+1)./j1).*x(:,i)).^2).^(1/2).*...
              (sin( (50.*((abs(x(:,i))-j1./j2).^2 + 10.*(x(:,i+1)- ((j1+1)./j1).*x(:,i)).^2).^(1/4))).^2+100) ;
      end
end

function y=Trigonometric(x,obj)
         N_q2=size(x,1);
         D_q2=size(x,2);
         y1=sum(cos(x(:,1:D_q2)),2);
         y=zeros(N_q2,1);
         for i=1:D_q2
             j=obj.K+obj.q1+i;
             y=y + ( D_q2.*ones(N_q2,1) - y1 + j.*(1-cos(x(:,i))-sin(x(:,i)))).^2;
         end
end 