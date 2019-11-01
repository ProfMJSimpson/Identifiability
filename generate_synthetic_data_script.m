%% Generate synthetic data

%% D same
D = 500
Dr = D
Dg = D
kr = 0.05
kg = 0.1

%% D different
Dr = 700
Dg = 300
kr = 0.05
kg = 0.1

%% Generate
% make sure use synthetic option is off for actually generating...
[ut1,ut2,ut3,ut4,data_exp0,data_exp16,...
    data_exp32,data_exp48] = forward_solve(Dr,Dg,kr,kg,false);

%coarse
data0=load('datat0.txt');
times = data0(:,1)
grid = linspace(1,size(ut1,1),size(data0,1));

%fine
%data0=load('interpolateddatat0.txt');
%times = data0(:,1)
%grid = linspace(1,size(ut1,1),size(data0,1));

datat0_syn = [times,ut1(grid,:)];
datat16_syn = [times,ut2(grid,:)];
datat32_syn = [times,ut3(grid,:)];
datat48_syn = [times,ut4(grid,:)];

%write to txt
writematrix(datat0_syn)
writematrix(datat16_syn)
writematrix(datat32_syn)
writematrix(datat48_syn)


