function [ut1,ut2,ut3,ut4,data_exp0,data_exp16,data_exp32,data_exp48] = forward_solve(Dr,Dg,kr,kg,do_synthetic)
%This MATLAB code reads in experimental data from Vittadello et al. (2018),
%solves the PDE model and plots the evolution of the solution of the PDE
%onto the data.  The code also calculates the total number of red and green
%cells from the numerical solution and plots that data
% Trial and error gives us Dg=Dr=700, kr=0.030 and kg=0.075 as excellent
% matches (visual, at least to me, August 6)
%
% OJM. 2 October 2019 - just wrapped Implicit_Solver_RG_LinearDiffusion in function.
% Returns solutions and data.
%
tic

L=1242; %length of domain
dx=1.0; %spatial discretisation
dt=1.0; %temporal discretisation
tf=48.0; %duration of the experiment
%Dr=695.0; %Diffusivity of the red population
%Dg=695.0; %Diffusivity of the green population
%kr=0.030; %rate of red-to-green
%kg=0.071; %rate of green-to-red
tol=1e-10;
maxsteps=tf/dt;
N=L/dx+1;
%n_r=zeros(4,1);
%n_g=zeros(4,1);
%n_tot=zeros(4,1);



u=zeros(N,2);   %this is the solution that is marched through time
ut1=zeros(N,2); %this is the solution stored at t=0
ut2=zeros(N,2); %this is the solution stored at t=16
ut3=zeros(N,2); %this is the solution stored at t=32
ut4=zeros(N,2); %This is the solution stored at t=48
pu=zeros(N,2); %previous solution that is stored
u_update=zeros(N,2);
x=zeros(N,1);


data_interp0=zeros(N,2);  %this is the experimental data interpolated at t=0
%data_interp16=zeros(N,2); %this is the experimental data interpolated at t=16
%data_interp32=zeros(N,2); %this is the experimental data interpolated at t=32
%data_interp48=zeros(N,2); %This is the experimental data interpolated at t=48


if do_synthetic == true
    %synthetic
    %'using synthetic'
    data0=load('datat0_syn.txt');
    data16=load('datat16_syn.txt');
    data32=load('datat32_syn.txt');
    data48=load('datat48_syn.txt');
    
else
    %actual
    data0=load('datat0.txt');
    data16=load('datat16.txt');
    data32=load('datat32.txt');
    data48=load('datat48.txt');
end

data_exp0 = data0(:,2:3);
data_exp16 = data16(:,2:3);
data_exp32 = data32(:,2:3);
data_exp48 = data48(:,2:3);


%NRed_data=[286 533 663 932];
%NGreen_data=[234 243 286 317];
%NTotal_data=NRed_data+NGreen_data;

xi=(0:dx:1242)';
data_interp0(:,1)=interp1q(data0(:,1),data0(:,2),xi); %this is the experimental red   data interpolated at t=0
data_interp0(:,2)=interp1q(data0(:,1),data0(:,3),xi); %this is the experimental green data interpolater at t=0

%data_interp16(:,1)=interp1q(data16(:,1),data16(:,2),xi); %this is the experimental red   data interpolated at t=16
%data_interp16(:,2)=interp1q(data16(:,1),data16(:,3),xi); %this is the experimental green data interpolater at t=16

%data_interp32(:,1)=interp1q(data32(:,1),data32(:,2),xi); %this is the experimental red   data interpolated at t=32
%data_interp32(:,2)=interp1q(data32(:,1),data32(:,3),xi); %this is the experimental green data interpolater at t=32

%data_interp48(:,1)=interp1q(data48(:,1),data48(:,2),xi); %this is the experimental red   data interpolated at t=48
%data_interp48(:,2)=interp1q(data48(:,1),data48(:,3),xi); %this is the experimental green data interpolater at t=48


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
%for i=1:N-1
%n_r(1,1)=n_r(1,1)+0.004*1745*dx*(ut1(i,1)+ut1(i+1,1))/2;
%n_g(1,1)=n_g(1,1)+0.004*1745*dx*(ut1(i,2)+ut1(i+1,2))/2;
%end
%n_tot(1,1)=n_g(1,1)+n_r(1,1);

for i=1:maxsteps
    
    t=t+dt;
    u_update=pu;
    diff=1.e10*ones(1,2);
    pic=0;
    
    while (norm(diff) > tol)
        pic = pic + 1;
        
        
        %Solve r(x,t) equation
        a(1,1)=0.0;
        b(1,1)=1.0; %no flux
        c(1,1)=-1.0; %no flux
        d(1,1)=0.0;
        
        
        a(N,1)=-1.0; % no flux
        b(N,1)=1.0; % no flux
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
    % if mod(i,100)==0
    % fprintf('Time %d\n',t);
    % fprintf('Iterations %d\n',pic);
    % end
    
    if abs(t-16)<tol
        ut2=u; %this is the solution at t = 16
    elseif abs(t-32)<tol
        ut3=u; %this is the solution at t = 32
    elseif abs(t-48)<tol
        ut4=u; % this is the solution at t=48
    end
    
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

%sol = [ut1,ut2,ut3,ut4];
end