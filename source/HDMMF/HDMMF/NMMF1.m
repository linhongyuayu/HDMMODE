classdef NMMF1 < PROBLEM
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
        K   
        q1;   
        q;
        Q1
    end
    
    methods



        %% 
        function Setting(obj)
            obj.M = 2;
            obj.K=1; 
            obj.q1=6; 
            obj.q=4;   

            if isempty(obj.D); obj.D = 7; end
            obj.lower(:,1:obj.K) = 0*ones(1,obj.K);
            obj.upper(:,1:obj.K) = 1*ones(1,obj.K);
            obj.lower(:,obj.K+1:obj.D) = -10.*ones(1,obj.D-obj.K);
            obj.upper(:,obj.K+1:obj.D) = 10.*ones(1,obj.D-obj.K);
            obj.encoding = 'real';
        end

        %% 
        function PopObj = CalObj(obj,PopDec)
            X_Q1=PopDec(:,obj.K+1:obj.K+obj.q1);
            obj.Q1=Pathological(X_Q1);
            G = 1 + obj.Q1;


            h = 1- (PopDec(:,1)./G).^2 - (PopDec(:,1)./G).*sin(2*pi.*obj.q.*PopDec(:,1)) ;
            PopObj(:,1) = PopDec(:,1);
            PopObj(:,2) = G.*h;
        end


        function R = GetOptimum(~,~)
            R = load('NMMF1_Reference_PSPF_data.mat','PF');
            R=R.PF;
        end

        function R = GetPF(obj)
            if obj.M == 2
                 R = load('NMMF1_Reference_PSPF_data.mat','draw_pf');
                 R=R.draw_pf;
            else
                R = [];
            end
        end

        function score = CalMetric(obj,metName,Population)
            load('NMMF1_Reference_PSPF_data');
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

function y=Pathological(x)
      [N1,q1]=size(x);
      y=zeros(N1,1);
      r=r_cal(x(:,1));
      for i=2:q1
          if i==q1
            x(:,i+1)=x(:,1);
          end
          y1= 50.*sin(sqrt(100.* (x(:,i+1)-x(:,i)).^2 + (x(:,i)-x(:,i-1)).^2)).^2 + 100*r -2;
          y2= 1 + 1.*( (x(:,i+1)-x(:,i)).^2 - 2.*(x(:,i+1)-x(:,i)).*(x(:,i)-x(:,i-1)) + (x(:,i)-x(:,i-1)).^2).^2;
          y = y  + y1./y2 +2 ; 
      end
end

%%     
function r=r_cal(x)
n=size(x,1);
for jj=1:n
    if x(jj,1)>5
        r(jj,:)=25.*(x(jj,1)-6).^2;
    elseif x(jj,1)<-5
        r(jj,:)=25./4.*(x(jj,1)+7).^2;
    else
        if x(jj,1)>=0
            r(jj,:)= abs((0.508333333333333).*x(jj,1).^3 - (1.641666666666664).*x(jj,1).^2 + 5/2);
        else
            r(jj,:)=(3/2).*(x(jj,1)+1).^2+1;
        end
    end
end
end