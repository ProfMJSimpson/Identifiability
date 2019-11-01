function obj = vector_objective(Dr,Dg,kr,kg,sd,do_synthetic,k_ln_scale)
    if k_ln_scale
        [ut1,ut2,ut3,ut4,data_exp0,data_exp16,data_exp32,data_exp48] = forward_solve(Dr,Dg,exp(kr),exp(kg),do_synthetic);
    else
        [ut1,ut2,ut3,ut4,data_exp0,data_exp16,data_exp32,data_exp48] = forward_solve(Dr,Dg,kr,kg,do_synthetic);
    end
    %split out in case want to weight differently etc.
    %grid = linspace(1,
    %assume uniform and extract sim at measured
    grid = linspace(1,size(ut1,1),size(data_exp0,1));
    
    d2 = ut2(grid,:)-data_exp16;
    d3 = ut3(grid,:)-data_exp32;
    d4 = ut4(grid,:)-data_exp48;
    
    %approx scale factor
    %data
    %sd = 0.02;
    %synthetic exact
    %sd = 1e-5;
    
    obj = [d2(:);d3(:);d4(:)]/sd;
end