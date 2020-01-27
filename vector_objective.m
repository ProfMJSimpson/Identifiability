function obj = vector_objective(Dr,Dg,kr,kg,sd,data_c_all,initial_condition,output_times_list,k_ln_scale)

    if k_ln_scale
        u_xtc_array = forward_solve(Dr,Dg,exp(kr),exp(kg),initial_condition,output_times_list);
    else
        u_xtc_array = forward_solve(Dr,Dg,kr,kg,initial_condition,output_times_list);
    end
    %split out in case want to weight differently etc.
    %grid = linspace(1,
    %assume uniform and extract sim at measured
    x_grid = linspace(1,size(u_xtc_array,1),size(data_c_all,1));
    %grid = linspace(1,size(ut1,1),size(data_exp0,1));
    
    residual = data_c_all - u_xtc_array(x_grid,:,:);
    d = reshape(residual,[numel(residual),1]);
    
    obj = d/sd;
    %d2 = ut2(x_grid,:)-data_exp16;
    %d3 = ut3(x_grid,:)-data_exp32;
    %d4 = ut4(x_grid,:)-data_exp48;
    
    %approx scale factor
    %data
    %sd = 0.02;
    %synthetic exact
    %sd = 1e-5;
    
    %obj = [d2(:);d3(:);d4(:)]/sd;
end