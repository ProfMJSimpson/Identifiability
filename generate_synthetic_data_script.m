%% Generate synthetic data

%Editing 17 Jan to make output at all times of interest.

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

%% Desired data

output_times_list = linspace(0,48,49);

%% Initial data
dx=1;
[data_xc0,data_c0,data_c0_interpolated] = load_exp_data_file('datat0.txt',dx);
initial_condition=data_c0_interpolated;

%% Generate
% make sure use synthetic option is off for actually generating...
[u_xtc_matrix] = forward_solve(Dr,Dg,kr,kg,initial_condition,output_times_list); %,false);

%% Processing

coarse_grid = false
if coarse_grid
    %coarse (real data)
    data0=load('datat0.txt');
    x = data0(:,1)
    x_grid = linspace(1,size(u_xtc_matrix,1),size(data0,1));
else
    data0=initial_condition;
    x = zeros(length(data0),1);
    x_grid = linspace(1,size(u_xtc_matrix,1),size(data0,1)); 
end
%fine
%data0=load('interpolateddatat0.txt');
%times = data0(:,1)
%grid = linspace(1,size(ut1,1),size(data0,1));

%% Extract at data

for i=1:length(output_times_list)

    data_xc_t_syn = [x,u_xtc_matrix(x_grid,i,1),u_xtc_matrix(x_grid,i,2)];
    %datat16_syn = [times,ut2(grid,:)];
    %datat32_syn = [times,ut3(grid,:)];
    %datat48_syn = [times,ut4(grid,:)];

    %write to txt
    writematrix(data_xc_t_syn,'data_xc_'+ string(output_times_list(i)) +'_syn.txt')
end
%writematrix(datat16_syn)
%writematrix(datat32_syn)
%writematrix(datat48_syn)


