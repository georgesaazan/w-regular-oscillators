%this function outputs bounds on the upper bounds of lambda(rho) using Theorem 2,
%it takes the following arguments respectively: a set of matrices (ceil), rho, the fixed error for the dichotomy.
%ex: bounds_lambda_regular(osc(3,0.2),1.02,10^-3).
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi (https://github.com/SQLP/SeDuMi).
% works if number of oscillators<12
function c=bounds_lambda_regular(A,rhoo,eps)
m=length(A);
if m<3
error("Please enter a number of oscillators greater or equal than 3")
end
n=size(A{1},1);
variableQ.('Q') = sdpvar(n,n);
q=fieldnames(variableQ);
qu=[];
for k=0:m-3
    q=fieldnames(variableQ);
    qu={};
    for i=1:length(q)
    if length(q{i}(2:end))==k
        qu=[qu;q{i}];
    end
    end
for j=1:length(qu)
for i=setdiff(1:m,num2str(qu{j}(2:end))- '0')
variableQ.(strcat(qu{j},num2str(i))) = sdpvar(n,n); %initial state, step 1
end
end
end
q=fieldnames(variableQ);
delta=containers.Map; 
mx=find(cellfun(@(c) length(c)==m-1 ,q));
for i=1:length(q)
    if ~isempty(intersect(1:m,num2str(q{i}(2:end))- '0'))
    for j=intersect(1:m,num2str(q{i}(2:end))- '0')
    delta=[delta;containers.Map([q{i},j],q{i},'UniformValues',false)];
    end
    end
    if i<=mx(1)-1
    for j=setdiff(1:m,num2str(q{i}(2:end))- '0')
    delta=[delta;containers.Map([q{i},j],strcat(q{i},num2str(j)),'UniformValues',false)];
    end
    else
            for j=setdiff(1:m,num2str(q{i}(2:end))- '0')
    delta=[delta;containers.Map([q{i},j],'Q','UniformValues',false)];
            end
    end
end
a=0;
b=1;
err=b-a;
while(err>eps)%dichotomy algorithm
lambda=(a+b)/2;
F=[];
for i=1:length(q)
F=[F variableQ.(q{i})>=eye(n)];  %step to construct the LMIs
end

for i=1:length(q)
    for j=1:m
        if ~strcmp(delta([q{i},j]),'Q')
         F=[F,A{j}'*variableQ.(delta([q{i},j]))*A{j} <=rhoo^2*variableQ.(q{i})];
        else
         F=[F,A{j}'*variableQ.(delta([q{i},j]))*A{j} <=rhoo^2*lambda^2*variableQ.(q{i})];
        end
    end
end
su=0;
for i=i:length(q)
su=su+trace(variableQ.(q{i}));
end
optimize(F,su,sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
if(min(primalfeas)<-1e-3)%1e-5
   a=lambda
   err=abs(lambda-b);
else %%feasable lmi
    b=lambda
    err=abs(lambda-a);
end 
end
c=b;%output bounds