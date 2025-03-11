function Offspring = OperatorDE_Single(Problem,Parent1,Parent2,Parent3,Parent4)


[CR,F,proM,disM] = deal(1,0.5,1,20);

if isa(Parent1(1),'SOLUTION')
    evaluated = true;
    Parent1   = Parent1.decs;
    Parent2   = Parent2.decs;
    Parent3   = Parent3.decs;
    Parent4   = Parent4.decs;
else
    evaluated = false;
end
[N,D] = size(Parent1);

%% Differental evolution
Site = rand(N,D) < CR;
Offspring       = Parent1;
Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent1(Site))   +    F*(Parent4(Site)-Parent3(Site));

%% Polynomial mutation
Lower = repmat(Problem.lower,N,1);
Upper = repmat(Problem.upper,N,1);
lower = repmat(Lower(1,:),N,1);
lower_index = Offspring < lower;
Offspring(lower_index) = (Parent1(lower_index) + lower(lower_index)) / 2;
upper = repmat(Upper(1,:),N,1);
upper_index = Offspring  > upper;
Offspring(upper_index) = (Parent1(upper_index) + upper(upper_index)) / 2;


mu    = rand(N,D);
temp  = Site & mu<=0.5;
Offspring       = min(max(Offspring,Lower),Upper);



Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
if evaluated
    Offspring = Problem.Evaluation(Offspring);
end
end