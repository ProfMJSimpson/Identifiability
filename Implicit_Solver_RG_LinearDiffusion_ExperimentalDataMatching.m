function [e] =Implicit_Solver_RG_LinearDiffusion_ExperimentalDataMatching(Dr, Dg, kr, kg);
L=1242; %length of domain
dx=1.0; %spatial discretisation
dt=1.0; %temporal discretisation
tf=48.0; %duration of the experiment
tol=1e-10; 
maxsteps=tf/dt;
N=L/dx+1;


u=zeros(N,2);   %this is the solution that is marched through time
ut1=zeros(N,3); %this is the solution stored at t=0
ut2=zeros(N,3); %this is the solution stored at t=16
ut3=zeros(N,3); %this is the solution stored at t=32
ut4=zeros(N,3); %This is the solution stored at t=48

dut1=zeros(24,3); %this is the solution stored at t=0 at the same spatial location as the experiments
dut2=zeros(24,3); %this is the solution stored at t=16 at the same spatial location as the experiments
dut3=zeros(24,3); %this is the solution stored at t=32 at the same spatial location as the experiments
dut4=zeros(24,3); %This is the solution stored at t=48 at the same spatial location as the experiments

pu=zeros(N,2); %previous solution that is stored
u_update=zeros(N,2);
x=zeros(N,1);

idatat1=load('interpolateddatat0.txt');
idatat2=load('interpolateddatat16.txt');
idatat3=load('interpolateddatat32.txt');
idatat4=load('interpolateddatat48.txt');

edatat1=load('datat0.txt');
edatat2=load('datat16.txt');
edatat3=load('datat32.txt');
edatat4=load('datat48.txt');

u(:,1)=idatat1(:,2); %this is the initial synthetic red data interpolated
u(:,2)=idatat1(:,3); %this is the initial synthetic green data inerpolated

pu=u;

for i=1:N
    x(i,1)=0+(i-1)*dx;
end

a=zeros(N,1);
b=zeros(N,1);
c=zeros(N,1);
d=zeros(N,1);

t=0.0;

ut1(:,1)=x(:,1); % this is the solution at t= 0
ut1(:,2)=u(:,1); % this is the solution at t= 0
ut1(:,3)=u(:,2); % this is the solution at t= 0


for i=1:maxsteps
    
    t=t+dt;
    u_update=pu;
    diff=1.e10*ones(1,2);
    pic=0;
    
    while (norm(diff) > tol)
    pic = pic + 1;
    
 
%Solve r(x,t) equation
a(1,1)=0.0;
b(1,1)=1.0;
c(1,1)=-1.0;
d(1,1)=0.0;


a(N,1)=-1.0;
b(N,1)=1.0;
c(N,1)=0.0;
d(N,1)=0.0;

for j=2:N-1
a(j,1)=Dr/dx^2;
b(j,1)=-1.0/dt-kr-2*Dr/dx^2;
c(j,1)=Dr/dx^2;
d(j,1)=-pu(j,1)/dt-2.0*kg*u(j,2)*(1.0-u(j,1)-u(j,2));
end

u_update(:,1) = thomas(N,a,b,c,d);

%Solve g(x,t) equation
a(1,1)=0.0;
b(1,1)=1.0;
c(1,1)=-1.0;
d(1,1)=0.0;


a(N,1)=-1.0;
b(N,1)=1.0;
c(N,1)=0.0;
d(N,1)=0.0;

for j=2:N-1
a(j,1)=Dg/dx^2;
b(j,1)=-1.0/dt-kg*(1-u_update(j,1)-u(j,2))-2*Dg/dx^2;
c(j,1)=Dg/dx^2;
d(j,1)=-pu(j,2)/dt-kr*u_update(j,1);
end

u_update(:,2) = thomas(N,a,b,c,d);

diff=abs(max(u-u_update));
u = u_update;

    end
    pu=u_update;
    
   
    if abs(t-16)<tol
        ut2(:,1)=x(:,1); % this is the solution at t= 16
        ut2(:,2)=u(:,1); % this is the solution at t= 16
        ut2(:,3)=u(:,2); % this is the solution at t= 16
    elseif abs(t-32)<tol
        ut3(:,1)=x(:,1); % this is the solution at t= 32
        ut3(:,2)=u(:,1); % this is the solution at t= 32
        ut3(:,3)=u(:,2); % this is the solution at t= 32
    elseif abs(t-48)<tol
        ut4(:,1)=x(:,1); % this is the solution at t= 48
        ut4(:,2)=u(:,1); % this is the solution at t= 48
        ut4(:,3)=u(:,2); % this is the solution at t= 48
    end 

end
    
for i=1:24
dut1(i,1)=(i-1)*54;   
dut1(i,2)=ut1(1+(i-1)*54,2);
dut1(i,3)=ut1(1+(i-1)*54,3);

dut2(i,1)=(i-1)*54;   
dut2(i,2)=ut2(1+(i-1)*54,2);
dut2(i,3)=ut2(1+(i-1)*54,3);

dut3(i,1)=(i-1)*54;   
dut3(i,2)=ut3(1+(i-1)*54,2);
dut3(i,3)=ut3(1+(i-1)*54,3);

dut4(i,1)=(i-1)*54;   
dut4(i,2)=ut4(1+(i-1)*54,2);
dut4(i,3)=ut4(1+(i-1)*54,3);
end
    
    


   %Now compute the error
   e=0;
   %e=e+norm(ut1-idatat1)+norm(ut2-idatat2)+norm(ut3-idatat3)+norm(ut4-idatat4);
   e=e+norm(dut1-edatat1)+norm(dut2-edatat2)+norm(dut3-edatat3)+norm(dut4-edatat4);
    
   
   

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
x=zeros(N,1);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end



