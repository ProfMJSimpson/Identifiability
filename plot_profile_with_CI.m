function [fig,CI_lower,CI_upper] = plot_profile_with_CI(psi,obj,psi_name,n_fine)


%get confidence intervals
%d = psi;
like = exp(-obj(:,1))/max(exp(-obj));
    
psi_fine = linspace(psi(1),psi(length(psi)),n_fine);

like_fine = interp1(psi,like,psi_fine,'pchip');

CI_lower = psi_fine(find(like_fine >= 0.15, 1, 'first'));
CI_upper = psi_fine(find(like_fine >= 0.15, 1, 'last' ));

fig = figure;
box on;
hold on;
plot(psi_fine,like_fine,'LineWidth',1)

if find(like_fine <= 0.15)
    xline(CI_lower,'r')
    xline(CI_upper,'r')
end

ylim([0,1.1])
ylabel(['Profile Relative Likelihood L_p(',psi_name,'; y^o)'])

xlabel(psi_name)
xlim([min(psi),max(psi)])

end

