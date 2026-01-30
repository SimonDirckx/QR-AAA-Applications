%% sep. rat. approx to a non-rational function
clear
run init.m
clc

% odd #pts
nx=401;
ny=200;
nz=200;

x = linspace(-1,1,nx);
y = linspace(0,1,ny);
z = linspace(0,1,nz);


mmax_vec = 5:5:50;
info_vec = cell(1,numel(mmax_vec));
err_off_grid_vec=zeros(numel(mmax_vec),1);
err_grid_vec=zeros(numel(mmax_vec),1);

f = @(x,y,z) sqrt(x.^2+2*y+3*z);
[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));
%%


for mm_ind=1:numel(mmax_vec)
    mmax = mmax_vec(mm_ind);
    disp("mmax = ")
    mmax
    
    options.mmax=mmax;
    options.sep=true;
    options.tolAAA = 1e-8;
    options.tol_qr = 1e-9; %avoid overfitting
    [rxy_interp,info_interp] = construct_multi_AAA_ID({x,y,z},F,f,options);
    info_vec{mm_ind} = info_interp;
    R = rxy_interp.eval({x,y,z});
    err_grid_vec(mm_ind) = norm(F(:)-R(:))/norm(F(:));
    xtest = linspace(-1,1.,351);
    ytest = linspace(0,1.,351);
    ztest = linspace(0,1.,351);
    

    [Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
    Fexact = f(Xtest,Ytest,Ztest);
    Finterp =  rxy_interp.eval({xtest,ytest,ztest});
    err_off_grid_vec(mm_ind) = norm(Finterp-Fexact,'fro')/norm(Fexact,"fro");
    
end
%%
degr_vec = cellfun(@(x) max([ numel(x.Supp{1}),numel(x.Supp{2}),numel(x.Supp{3}) ]),info_vec);

[~,ii] = unique(degr_vec,'stable');


figure(1);clf;
semilogy(degr_vec(ii),err_off_grid_vec(ii))
hold on
semilogy(degr_vec(ii),err_grid_vec(ii))
legend('off\_grid','grid')

%%
clear
clc

nx=401;
ny=200;
nz=200;

x = linspace(-1,1,nx);
y = linspace(0,1,ny);
z = linspace(0,1,nz);

f = @(x,y,z) sqrt(x.^2+2*y+3*z);
[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));

options.tolAAA = 1e-8;
options.tol_qr = 1e-9;
[rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);
R = rxy_interp.eval({x,y,z});
err_grid = norm(F(:)-R(:))/norm(F(:))
xtest = linspace(-1,1.,351);
ytest = linspace(0,1.,351);
ztest = linspace(0,1.,351);


[Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
Fexact = f(Xtest,Ytest,Ztest);
Finterp =  rxy_interp.eval({xtest,ytest,ztest});
err_off_grid = norm(Finterp-Fexact,inf)/norm(Fexact,inf)
%% isolate tol = 10^-8
clear
nx=301;
ny=150;
nz=150;

x = linspace(-1,1,nx);
y = linspace(0,1,ny);
z = linspace(0,1,nz);


f = @(x,y,z) sqrt(x.^2+2*y+3*z);
[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));

options.tolAAA = 1e-8;
options.tol_qr = 1e-9; 
options.twostep = true;
[rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);

xtest = linspace(-1,1.,351);
ytest = linspace(0,1.,351);
ztest = linspace(0,1.,351);


[Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
Fexact = f(Xtest,Ytest,Ztest);
Finterp =  rxy_interp.eval({xtest,ytest,ztest});
norm(Finterp-Fexact,'fro')/norm(Fexact,"fro")
