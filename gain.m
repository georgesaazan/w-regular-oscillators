%this provides the results of the numerical example of the chapter,
%this requires the JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi(https://github.com/SQLP/SeDuMi),
%% In the first part we give the design figure introduced in the example
%remark 1: the simulation may pause multiple times with "Press any key to
%proceed", to override this, comment out a "pause" in the JSR toolbox:
%JSR_louvain/Methods/jsr.m line 357.
%remark 2: the simulation may take a long time, however the results are saved
%in result.m, by default this part will be skipped if "result.mat" is present in this directory.
clear all;
gamma=0.1;m=3;
if ~isfile('result.mat')
result={};
for m=3:5
resultx=[]; 
for gamma=0:0.01:0.9 %we change the oscillators gain \gamma from 0 to 0.9
    if(gamma<=0.8) %we verify that for \gamma<=0.8 the jsr is 1.02 for all n
        rho=1.02;
    else
    t=jsr(oscillators(m,gamma)); % for each \gamma and n, we generate a set of matrices using oscillators.m
rho=t.bounds(1); %if \gamma>0.8 we calculate the jsr using the jsr toolbox
    end
b=bounds_lambda_2(oscillators(m,gamma),rho,10^-4) %upper bounds on \lambda 
x=-log(rho)/log(b); %we verify numerically that the infimum is reached at the JSR when x<1/(m-1)
if(x<1/(m-1)&&b~=1)
else x=1/(m-1); %we verify that when x>=1/(m-1) the infimum is 1/(m-1).
end 
resultx=[resultx ,x];
end
result=[result,resultx]
end
else
result=load('result.mat');
end
result=result.result;
C = {'r','b','g'}; %define colors used in the figure
for i=1:3
a(i)=plot(0:0.01:0.9,result{i},'color',C{i});
yline(1/(i+1),'--','LineWidth',2,'color',C{i});
xlabel("k");
ylabel("Inf -log(\rho)/log(\lambda)");
title("stabilizability vs k");
ylim([0 0.6]);
hold on;
end
legend([a(1),a(2),a(3)],'m=3 oscillators','m=4 oscillators','m=5 oscillators');
%% Construction of the automaton
A=osc(m,gamma);
m=length(A);
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
%% Simulation fast 
N=100; %number of steps
A=osc(3,0.1); 
x=[-1.5;-0.5;2;-1]; %initial state
traj=x; 
sigma=[];
k=[0];
st={'Q'};
for i=1:N
sigma=[sigma randi([1,m])];    
end
for i=1:N
traj=[traj A{sigma(i)}*traj(:,end)];
st=[st delta([st{end},sigma(i)])];
k=[k k(end)+strcmp(st{end},'Q')];
end
figure();
subplot(3,1,1);
stairs([0:N-1],sigma)
axis([0 N 0.9 3.1])
xlabel('t')
ylabel('\theta(t)')
subplot(3,1,2);
plot(0:N-1,k(1:N)./[1:N],'--')
xlabel('t')
ylabel('\kappa^{\theta(t)}/t')
ylim([0,0.5]);
subplot(3,1,3);
plot(0:N-1,traj(1,1:N),'--')
hold on;
plot(0:N-1,traj(2,1:N),'--')
hold on;
plot(0:N-1,traj(3,1:N),'--')
hold on;
plot(0:N-1,traj(4,1:N),'--')
xlabel('t')
ylabel('traj')
%% Simulation slow
traj=x; 
k=[0];
sigma=[];
st={'Q'};
for i=1:N
adm0=[];
for j=1:m
if ~strcmp(delta([st{end},j]),'Q')
adm0=[adm0 j];
end
end
if mod(i,25)==1
    adm0=[randi(m)];
end
sigma=[sigma adm0(randi(length(adm0)))];
traj=[traj A{sigma(end)}*traj(:,end)];
st=[st delta([st{end},sigma(end)])];
k=[k k(end)+strcmp(st{end},'Q')];
end
figure();
subplot(3,1,1);
stairs([0:N-1],sigma)
axis([0 N 0.9 3.1])
xlabel('t')
ylabel('\theta(t)')
subplot(3,1,2);
plot(0:N-1,k(1:N)./[1:N],'--')
xlabel('t')
ylabel('\kappa^{\theta(t)}/t')
subplot(3,1,3);
plot(0:N-1,traj(1,1:N),'--')
hold on;
plot(0:N-1,traj(2,1:N),'--')
hold on;
plot(0:N-1,traj(3,1:N),'--')
hold on;
plot(0:N-1,traj(4,1:N),'--')
xlabel('t')
ylabel('traj')