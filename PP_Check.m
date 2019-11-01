NSamples=30;
opac=0.10;
Burned_kg=zeros(50000,1);
Burned_kr=zeros(50000,1);
Burned_Dr=zeros(50000,1);



Burned_kg=kg(41000:1:51000,1);
Burned_kr=kr(41000:1:51000,1);
Burned_Dr=Dr(41000:1:51000,1);
for i = 1:NSamples
   ind = randperm(numel(Burned_kg), 1); % select one element out of numel(x) elements, with the probability of occurrence of the element in x
   kg_Sample(i,1) = Burned_kg(ind);
   kr_Sample(i,1) = Burned_kr(ind);
   Dr_Sample(i,1) = Burned_Dr(ind);
   Burned_kg(Burned_kg==kg_Sample(i,1)) = []; % delete this element from the sample, such that the picked elements are unique
   Burned_kr(Burned_kr==kr_Sample(i,1)) = []; % delete this element from the sample, such that the picked elements are unique
   Burned_Dr(Burned_Dr==Dr_Sample(i,1)) = []; % delete this element from the sample, such that the picked elements are unique
end
 
for iii=1:NSamples

L=1242; %length of domain
dx=1.0; %spatial discretisation
dt=1.0; %temporal discretisation
tf=48.0; %duration of the experiment
Dr=Dr_Sample(iii,1); %Diffusivity of the red population
Dg=Dr_Sample(iii,1); %Diffusivity of the green population
kr=kr_Sample(iii,1); %rate of red-to-green
kg=kg_Sample(iii,1); %rate of green-to-red
tol=1e-10; 
maxsteps=tf/dt;
N=L/dx+1;

u=zeros(N,2);   %this is the solution that is marched through time
ut1=zeros(N,2); %this is the solution stored at t=0
ut2=zeros(N,2); %this is the solution stored at t=16
ut3=zeros(N,2); %this is the solution stored at t=32
ut4=zeros(N,2); %This is the solution stored at t=48
pu=zeros(N,2); %previous solution that is stored
u_update=zeros(N,2);
x=zeros(N,1);


data_interp0=zeros(N,2);  %this is the experimental data interpolated at t=0
data_interp16=zeros(N,2); %this is the experimental data interpolated at t=16
data_interp32=zeros(N,2); %this is the experimental data interpolated at t=32
data_interp48=zeros(N,2); %This is the experimental data interpolated at t=48


data0=load('datat0.txt');
data16=load('datat16.txt');
data32=load('datat32.txt');
data48=load('datat48.txt');

xi=(0:dx:1242)';
data_interp0(:,1)=interp1q(data0(:,1),data0(:,2),xi); %this is the experimental red   data interpolated at t=0
data_interp0(:,2)=interp1q(data0(:,1),data0(:,3),xi); %this is the experimental green data interpolater at t=0

data_interp16(:,1)=interp1q(data16(:,1),data16(:,2),xi); %this is the experimental red   data interpolated at t=16
data_interp16(:,2)=interp1q(data16(:,1),data16(:,3),xi); %this is the experimental green data interpolater at t=16

data_interp32(:,1)=interp1q(data32(:,1),data32(:,2),xi); %this is the experimental red   data interpolated at t=32
data_interp32(:,2)=interp1q(data32(:,1),data32(:,3),xi); %this is the experimental green data interpolater at t=32

data_interp48(:,1)=interp1q(data48(:,1),data48(:,2),xi); %this is the experimental red   data interpolated at t=48
data_interp48(:,2)=interp1q(data48(:,1),data48(:,3),xi); %this is the experimental green data interpolater at t=48


u(:,1)=data_interp0(:,1); %this is the initial experimental red data interpolated
u(:,2)=data_interp0(:,2); %this is the initial experimental green data inerpolated

pu=u;

for i=1:N
    x(i,1)=0+(i-1)*dx;
end

a=zeros(N,1);
b=zeros(N,1);
c=zeros(N,1);
d=zeros(N,1);

t=0.0;

ut1=u; % this is the solution at t= 0

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
        ut2=u; %this is the solution at t = 16
    elseif abs(t-32)<tol
        ut3=u; %this is the solution at t = 32
    elseif abs(t-48)<tol
    ut4=u; % this is the solution at t=48
    end 

end
    
    

    
    hold on
    subplot(2,2,1)
    p1= plot(x,ut1(:,1),'k-','LineWidth',1); 
    hold on
    p2= plot(x,ut1(:,2),'k-','LineWidth',1);
    p1.Color(4) = opac;
    p2.Color(4) = opac;
    xlabel('position')
    ylabel('density')
    title('Time 0 h')
    hold on
    axis([0 1242 0 1.2*max(ut4(:,1)+ut4(:,2))])
    
    
    subplot(2,2,2)
    p1= plot(x,ut2(:,1),'k-','LineWidth',1); 
    hold on
    p2= plot(x,ut2(:,2),'k-','LineWidth',1);
    p1.Color(4) = opac;
    p2.Color(4) = opac;
    xlabel('position')
    ylabel('density')
    title('Time 16 h')
    hold on
    axis([0 1242 0 1.2*max(ut4(:,1)+ut4(:,2))])
    
    
    subplot(2,2,3)
    p1= plot(x,ut3(:,1),'k-','LineWidth',1); 
    hold on
    p2= plot(x,ut3(:,2),'k-','LineWidth',1);
    p1.Color(4) = opac;
    p2.Color(4) = opac;
    xlabel('position')
    ylabel('density')
    title('Time 32 h')
    hold on
    axis([0 1242 0 1.2*max(ut4(:,1)+ut4(:,2))])
    
    
    subplot(2,2,4)
    p1= plot(x,ut4(:,1),'k-','LineWidth',1); 
    hold on
    p2= plot(x,ut4(:,2),'k-','LineWidth',1);
    p1.Color(4) = opac;
    p2.Color(4) = opac;
    xlabel('position')
    ylabel('density')
    title('Time 48 h')
    hold on
    axis([0 1242 0 1.2*max(ut4(:,1)+ut4(:,2))])

 

end
hold on
subplot(2,2,1)
plot(data0(:,1),data0(:,2),'ro','MarkerFaceColor','r')
plot(data0(:,1),data0(:,3),'go','MarkerFaceColor','g')
 subplot(2,2,2)
 plot(data16(:,1),data16(:,2),'ro','MarkerFaceColor','r')
 plot(data16(:,1),data16(:,3),'go','MarkerFaceColor','g')
 subplot(2,2,3)
plot(data32(:,1),data32(:,2),'ro','MarkerFaceColor','r')
plot(data32(:,1),data32(:,3),'go','MarkerFaceColor','g')

 subplot(2,2,4)
 plot(data48(:,1),data48(:,2),'ro','MarkerFaceColor','r')
 plot(data48(:,1),data48(:,3),'go','MarkerFaceColor','g')
 box on

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




