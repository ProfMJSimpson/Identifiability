function [data_xc,data_c,data_c_interp] = load_exp_data_file(filename,dx)

data_xc=load(filename);
data_c = data_xc(:,2:3);

xi=(0:dx:1242)';
data_c_interp = zeros(length(xi),2);
data_c_interp(:,1)=interp1q(data_xc(:,1),data_xc(:,2),xi); %this is the experimental red data interpolated
data_c_interp(:,2)=interp1q(data_xc(:,1),data_xc(:,3),xi); %this is the experimental green data interpolated

end

