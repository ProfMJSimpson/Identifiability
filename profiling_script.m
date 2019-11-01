

%% General implementation of profiling for target parameter
% Here: theta = (psi, lambda) i.e. (target,nuisance).
% Currently uses least squares formulation only.
% Each parameter is implemented separately below, calling same functions.

%% Overall max likelihood

options = optimoptions('lsqnonlin')
options.Display = 'iter';

k_ln_scale = true
if k_ln_scale
    theta0 = [600,600,log(0.03),log(0.07)]
    theta_lower_bounds = [0,0,log(0),log(0)]
    theta_upper_bounds = [5000,5000,5,5]
else
    theta0 = [600,600,0.03,0.07]
    theta_lower_bounds = [0,0,0,0]
    theta_upper_bounds = [5000,5000,1,1]
end

do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 1e-3 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end


f_theta = @(theta) vector_objective(theta(1),theta(2),theta(3),theta(4),...
    sd,do_synthetic,k_ln_scale)

[theta, fval] = lsqnonlin(f_theta,theta0,theta_lower_bounds,...
    theta_upper_bounds,options)

if k_ln_scale
    theta(3:4) = exp(theta(3:4));
end
theta

%% Dr
% target parameter details
psi_name = 'Dr'
n_psi = 100
%n_psi = 50
psi_min = 0
psi_max = 2000
psi_range = linspace(psi_min,psi_max,n_psi)

% nuisance parameter details
k_ln_scale = true
if k_ln_scale
    lambda_guess = [600,log(0.03),log(0.07)]
    lambda_lower_bounds = [0,log(0),log(0)]
    lambda_upper_bounds = [5000,5,5]
else
    lambda_guess = [600,0.03,0.07]
    lambda_lower_bounds = [0,0,0]
    lambda_upper_bounds = [5000,1,1]
end

% synthetic or not
do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 0.001 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end

% mapping (psi, lambda) -> (theta). Key user input
f_theta = @(psi, lambda) vector_objective(psi,lambda(1),lambda(2),lambda(3),sd,do_synthetic,k_ln_scale);

% display and plot options
do_display_it = true

% call main profiling function
[lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

% compute CIs and plot profile likelihood
n_fine = 10000
[fig,lower,upper] = plot_profile_with_CI(psi_range,obj,psi_name,n_fine)
if do_synthetic
    filename = ['./figures/',psi_name,'_synthetic']%,string(do_synthetic)]
else
    filename = ['./figures/',psi_name]%,string(do_synthetic)]
end
savefig(fig,filename)

% plot profile likelihood value vs lambda estimate of Dg (should be similar
% to below if ridge)
%n_fine = 10000
%plot_profile_with_CI(lambda_estimates(:,1)',obj,'Dg (nuisance) estimate',n_fine)

%% Dg
% target parameter details
psi_name = 'Dg'
%n_psi = 50
n_psi = 100
psi_min = 0
psi_max = 2000
psi_range = linspace(psi_min,psi_max,n_psi)

% nuisance parameter details
k_ln_scale = true
if k_ln_scale
    lambda_guess = [600,log(0.03),log(0.07)]
    lambda_lower_bounds = [0,log(0),log(0)]
    lambda_upper_bounds = [5000,5,5]
else
    lambda_guess = [600,0.03,0.07]
    lambda_lower_bounds = [0,0,0]
    lambda_upper_bounds = [5000,1,1]
end

% synthetic or not
do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 0.001 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end

% mapping (psi, lambda) -> (theta). Key user input
f_theta = @(psi, lambda) vector_objective(lambda(1),psi,lambda(2),lambda(3),sd,do_synthetic,k_ln_scale);

% display and plot options
do_display_it = true

% call main profiling function
[lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

% compute CIs and plot profile likelihood
n_fine = 10000
[fig,lower,upper] = plot_profile_with_CI(psi_range,obj,psi_name,n_fine)
if do_synthetic
    filename = ['./figures/',psi_name,'_synthetic'] ;%,string(do_synthetic)]
else
    filename = ['./figures/',psi_name];%,string(do_synthetic)]
end
savefig(fig,filename)

% plot profile likelihood value vs lambda estimate of Dr (should be similar
% to above if ridge)
%n_fine = 10000
%plot_profile_with_CI(lambda_estimates(:,1)',obj,'Dr (nuisance) estimate',n_fine)

%% kr
% target parameter details
psi_name = 'kr'

%k_ln_scale = false
k_ln_scale = true

%n_psi = 50
n_psi = 100

if k_ln_scale
    psi_min = log(0.01)
    psi_max = log(0.05)
    psi_range = linspace(psi_min,psi_max,n_psi)
else
    psi_min = 0.001
    psi_max = 0.05
    psi_range = linspace(psi_min,psi_max,n_psi) 
end
    
% nuisance parameter details
if k_ln_scale
    lambda_guess = [600,600,log(0.07)]
    lambda_lower_bounds = [0,0,log(0)]
    lambda_upper_bounds = [5000,5000,5]
else
    lambda_guess = [600,600,0.07]
    lambda_lower_bounds = [0,0,0]
    lambda_upper_bounds = [5000,5000,1]
end

% synthetic or not
do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 0.001 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end

% mapping (psi, lambda) -> (theta). Key user input
f_theta = @(psi, lambda) vector_objective(lambda(1),lambda(2),psi,lambda(3),sd,do_synthetic,k_ln_scale);

% display and plot options
do_display_it = true

% call main profiling function
[lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

% compute CIs and plot profile likelihood
n_fine = 10000
if k_ln_scale
    [fig,lower,upper] = plot_profile_with_CI(exp(psi_range),obj,psi_name,n_fine);
else
    [fig,lower,upper] = plot_profile_with_CI(psi_range,obj,psi_name,n_fine)
end
if do_synthetic
    filename = ['./figures/',psi_name,'_synthetic'] ;%,string(do_synthetic)]
else
    filename = ['./figures/',psi_name];%,string(do_synthetic)]
end
savefig(fig,filename)

%% kg
% target parameter details
psi_name = 'kg'

%k_ln_scale = false
k_ln_scale = true

n_psi = 100 %50

if k_ln_scale
    psi_min = log(0.04)
    psi_max = log(0.1)
    psi_range = linspace(psi_min,psi_max,n_psi)
else
    psi_min = 0.04
    psi_max = 0.1
    psi_range = linspace(psi_min,psi_max,n_psi) 
end
    
% nuisance parameter details
if k_ln_scale
    lambda_guess = [600,600,log(0.03)]
    lambda_lower_bounds = [0,0,log(0)]
    lambda_upper_bounds = [5000,5000,5]
else
    lambda_guess = [600,600,0.03]
    lambda_lower_bounds = [0,0,0]
    lambda_upper_bounds = [5000,5000,1]
end

% synthetic or not
do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 0.001 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end

% mapping (psi, lambda) -> (theta). Key user input
f_theta = @(psi, lambda) vector_objective(lambda(1),lambda(2),lambda(3),psi,sd,do_synthetic,k_ln_scale);

% display and plot options
do_display_it = true

% call main profiling function
[lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

% compute CIs and plot profile likelihood
n_fine = 10000
if k_ln_scale
    [fig,lower,upper] = plot_profile_with_CI(exp(psi_range),obj,psi_name,n_fine);
else
    [fig,lower,upper] = plot_profile_with_CI(psi_range,obj,psi_name,n_fine)
end
if do_synthetic
    filename = ['./figures/',psi_name,'_synthetic'] ;%,string(do_synthetic)]
else
    filename = ['./figures/',psi_name];%,string(do_synthetic)]
end
savefig(fig,filename)

%% Dr-Dg
% target parameter details
psi_name = 'Dr-Dg'
n_psi = 100
psi_min = -1000
psi_max = 1000
psi_range = linspace(psi_min,psi_max,n_psi)

% nuisance parameter details
k_ln_scale = false
%k_ln_scale = true

if k_ln_scale
    lambda_guess = [2000,log(0.03),log(0.07)]
    lambda_lower_bounds = [0,log(0),log(0)]
    lambda_upper_bounds = [5000,5,5]
else
    lambda_guess = [2000,0.03,0.07]
    lambda_lower_bounds = [0,0,0]
    lambda_upper_bounds = [5000,1,1]
end

% synthetic or not
do_synthetic = true
%do_synthetic = false
if do_synthetic
    sd = 0.001 %synthetic approx
else
    sd = 0.05
    %sd = 0.02
end

% mapping (psi, lambda) -> (theta). Key user input
f_theta = @(psi, lambda) vector_objective(lambda(1)+psi,lambda(1),lambda(2),lambda(3),sd,do_synthetic,k_ln_scale);

% display and plot options
do_display_it = true;

% call main profiling function
[lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

% compute CIs and plot profile likelihood
n_fine = 10000
[fig,lower,upper] = plot_profile_with_CI(psi_range,obj,psi_name,n_fine)
if do_synthetic
    filename = ['./figures/',psi_name,'_synthetic'] ;%,string(do_synthetic)]
else
    filename = ['./figures/',psi_name];%,string(do_synthetic)]
end
savefig(fig,filename)

