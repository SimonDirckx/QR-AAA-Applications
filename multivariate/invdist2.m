%% sep. rat. approx to a non-rational function
clear
run init.m
clc

%odd #pts
nx=151;
ny=151;
nz=151;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);

powvec = 1:10;
err_grid_vec = zeros(numel(powvec),1);
err_off_grid_vec = zeros(numel(powvec),1);
info_vec = cell(numel(powvec),1);
for powind=1:numel(powvec)
    pow = powvec(powind);
    del = 2^(-pow);
    f = @(x,y,z) 1./sqrt(x.^2+2*y.^2+3*z.^2+del);
    
    [X,Y,Z] = ndgrid(x,y,z);
    F = tensor(f(X,Y,Z));
    options.interp=true;
    options.sep=true;
    options.tolAAA = 1e-11;
    options.tol_qr = 1e-15;
    [rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);
    info_vec{powind} = info_interp;
    R = rxy_interp.eval({x,y,z});

    err_grid_vec(powind) = norm(F(:)-R(:))/norm(F(:));
    
    xtest = linspace(-1,1.,300);
    ytest = linspace(-1,1.,300);
    ztest = linspace(-1,1.,300);
    [Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);

    Fexact = f(Xtest,Ytest,Ztest);
    
    Finterp =  rxy_interp.eval({xtest,ytest,ztest});
    
    err_off_grid_vec(powind) = norm(Finterp-Fexact,'fro')/norm(Fexact,"fro");

end

%%
figure(1);clf;
semilogy(err_grid_vec)
hold on
semilogy(err_off_grid_vec)

degr_x_vec = cellfun(@(x) numel(x.Supp{1}),info_vec);
degr_y_vec = cellfun(@(x) numel(x.Supp{2}),info_vec);
degr_z_vec = cellfun(@(x) numel(x.Supp{3}),info_vec);

figure(2);clf;
plot(degr_x_vec)
hold on
plot(degr_y_vec)
plot(degr_z_vec)
