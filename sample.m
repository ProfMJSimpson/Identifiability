%In this test case we use the MC to explore the identifiability with the
%constraint that Dr=Dg=D.  
%I have run the code several times with very differnt initial conditions
%for the MC and we appear to get very nice results
%kr=0.0306; %seems to be a good value of red-to-green
%kg=0.07065; %seems to be a good value of green-to-red
%Dr=688; %seems to be a good value of red diffusivity 
%Dg=688; %seems to be a good value of green diffusivity 
sigma = 0.05  ; %this is the value of the variance in the likelihood evaluation
Maxit=51000;
kg=zeros(Maxit,1);
kr=zeros(Maxit,1);
Dg=zeros(Maxit,1);
Dr=zeros(Maxit,1);
lhd=zeros(Maxit,1);

Dg(1,1)=0; %initial guess of Dg
Dr(1,1)=866; %initial guess of Dr
kg(1,1)=0.0552; %initial guess of kg
kr(1,1)=0.0096; %initial guess of kr
[e] =Implicit_Solver_RG_LinearDiffusion_ExperimentalDataMatching(Dr(1,1), Dr(1,1), kr(1,1), kg(1,1));
lhd(1,1)=exp(-0.5*((e)/sigma)^2)/(sigma*sqrt(2*pi));

mu=[0,0,0,0];
k=2;
sigma_prop=[100.0,0,0,0;  0,0.0,0,0;  0,0,0.000001,0;  0,0,0,0.000001];

while k < Maxit+1
R = mvnrnd(mu,sigma_prop,1);
Dr(k,1)=max(Dr(k-1,1)+R(1,1),0);
Dg(k,1)=max(Dg(k-1,1)+R(1,2),0);
kr(k,1)=max(kr(k-1,1)+R(1,3),0);
kg(k,1)=max(kg(k-1,1)+R(1,4),0);
[e] =Implicit_Solver_RG_LinearDiffusion_ExperimentalDataMatching(Dr(k,1), Dr(k,1), kr(k,1), kg(k,1));
lhd(k,1)=exp(-0.5*((e)/sigma)^2)/(sigma*sqrt(2*pi));    
alpha = min(1,lhd(k,1)/lhd(k-1,1));
RR=rand;    

    if RR <= alpha %accept this move
        k=k+1;
        
     if mod(k,100)==0
     fprintf('Iteration %d\n',k);
     e
     end   
        continue
    else
        %reject this move and try again
        Dr(k,1)=Dr(k-1,1);
        Dg(k,1)=Dg(k-1,1);
        kr(k,1)=kr(k-1,1);
        kg(k,1)=kg(k-1,1);
    end
   
end

%save('Markov_Chain_datax.mat','Dr','kg', 'kr', 'lhd')

x=1:1:Maxit;
figure
plot(x,Dr(:,1),'k')
title('Diffusivity Estimate')
xlabel('Markov Chain Iteration')
ylabel('D')
axis([0 Maxit 0 2000])

figure
plot(x,kr(:,1),'k')
title('Cell Cycle Rate Estimate')
xlabel('Markov Chain Iteration')
ylabel('kr')
axis([0 Maxit 0 0.1])



figure
plot(x,kg(:,1),'k')
title('Cell Cycle Rate Estimate')
xlabel('Markov Chain Iteration')
ylabel('kg')
axis([0 Maxit 0 0.1])

figure
[f,xi] = ksdensity(kr(1000:1:Maxit),'Bandwidth',0.001); 
[L,M,U]=CredibleInterval(xi,f)
plot(xi,f);
xx=[L,L];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy, 'r');
xx=[U,U];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'r');
xx=[M,M];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'g');
axis([0 0.1 0 1.2*max(f)])
ylabel('P(kr)')
xlabel('kr')


figure
[f,xi] = ksdensity(kg(1000:1:Maxit),'Bandwidth',0.001); 
[L,M,U]=CredibleInterval(xi,f)
plot(xi,f);
xx=[L,L];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy, 'r');
xx=[U,U];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'r');
xx=[M,M];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'g');
axis([0 0.1 0 1.2*max(f)])
ylabel('P(kg)')
xlabel('kg')

figure
[f,xi] = ksdensity(Dr(1000:1:Maxit),'Bandwidth',50.0); 
[L,M,U]=CredibleInterval(xi,f)
plot(xi,f);
xx=[L,L];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy, 'r');
xx=[U,U];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'r');
xx=[M,M];
yy=[0, 2.0*max(f)];
hold on
plot(xx,yy ,'g');
axis([0 2000 0 1.2*max(f)])
ylabel('Pr(D)')
xlabel('D')



figure
XX=[Dr(1000:1:Maxit),kg(1000:1:Maxit)];
[n,c] = hist3(XX, [10 10]);
contourf(c{1}, c{2}, n,'LineColor','none');
%yticks([ 0.025 0.030  ])
%xticks([ 0.070 0.075 ])
colormap(gray)
colormap(flipud(gray))
shading('interp')
ylabel('kg')
xlabel('D')
axis([0 2000 0 0.1])

box on

 
figure
XX=[Dr(1000:1:Maxit),kr(1000:1:Maxit)];
[n,c] = hist3(XX, [10 10]);
contourf(c{1}, c{2}, n,'LineColor','none');
%yticks([ 0.025 0.030  ])
%xticks([ 0.070 0.075 ])
colormap(gray)
colormap(flipud(gray))
shading('interp')
ylabel('kr')
xlabel('D')
axis([0 2000 0 0.1])

box on

figure
XX=[kr(1000:1:Maxit),kg(1000:1:Maxit)];
[n,c] = hist3(XX, [10 10]);
contourf(c{1}, c{2}, n,'LineColor','none');
%yticks([ 0.025 0.030  ])
%xticks([ 0.070 0.075 ])
colormap(gray)
colormap(flipud(gray))
shading('interp')
ylabel('kg')
xlabel('kr')
axis([0 0.1 0 0.1])

box on
 