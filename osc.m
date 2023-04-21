%this function takes 2 arguments, m the number of oscillators and their
%fixed gain k.
%this function outputs the m state matrices for m oscillators.
%ex: oscillators(3,0.2).
function AA=osc(m,gamma)
phi=pi/6;
nu=1.02;
R=nu*[cos(phi) -sin(phi);
   sin(phi)   cos(phi)];
A={};
B=[];
K=gamma*eye(2);
for i=1:m-1
A=kron(eye(m-1),R);
A(2*i-1:2*i,2*i-1:2*i)=A(2*i-1:2*i,2*i-1:2*i)-2*K;
if i==1
A(2*i+1:2*i+2,2*i-1:2*i)=K;
elseif i==m-1
A(2*i-3:2*i-2,2*i-1:2*i)=K;
else
A(2*i+1:2*i+2,2*i-1:2*i)=K;  
A(2*i-3:2*i-2,2*i-1:2*i)=K;
end
B=[B {A}];
end
A=kron(eye(m-1),R);
KK=[];
for i=1:m-1
    KK=[KK,K];
end
A(1:2,1:end)=A(1:2,1:end)-KK;
A(end-1:end,1:end)=A(size(A,1)-1:size(A,1),1:end)-KK;
B=[B {A}];
AA=B;