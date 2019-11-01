function [lambda_estimates,obj] = compute_profile_likelihood(f_theta,psi_range,...
    lambda_guess,lambda_lower_bounds,lambda_upper_bounds,do_display_it)

%theta = (psi,lambda), psi is target; lambda is nuisance and optimised out.

n_psi = length(psi_range);
obj = zeros(n_psi,1);
lambda_estimates = zeros(n_psi,size(lambda_guess,2));

%try least squares
options = optimoptions('lsqnonlin');%,'Algorithm','levenberg-marquardt');
if do_display_it
    options.Display = 'iter';
end

for i = 1:n_psi
    if do_display_it
       i 
    end
    psi_i = psi_range(i);
    f_lambda_given_psi = @(lambda) f_theta(psi_i,lambda);
    if i > 1
        lambda_guess = lambda_estimates(i-1,:);
    end
    [lambda_estimates(i,:),obj(i)] = lsqnonlin(f_lambda_given_psi,...
        lambda_guess,lambda_lower_bounds,lambda_upper_bounds,options);
end

end