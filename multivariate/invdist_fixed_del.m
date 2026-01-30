%% sep. rat. approx to a non-rational function
clear
run init.m
clc

% odd #pts
nx=151;
ny=151;
nz=151;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);


mmax_vec = 5:5:75;
err_off_grid_vec=zeros(numel(mmax_vec),1);
err_grid_vec=zeros(numel(mmax_vec),1);
info_vec = cell(1,numel(mmax_vec));

del = 1/64;
f = @(x,y,z) 1./sqrt(x.^2+2*y.^2+3*z.^2+del);

[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));

for mm_ind=1:numel(mmax_vec)
    mmax = mmax_vec(mm_ind);
    disp("mmax = ")
    mmax
    options.mmax=mmax;
    options.sep=true;
    options.tolAAA = 1e-11;
    options.tol_qr = 1e-15;
    [rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);
    info_vec{mm_ind} = info_interp;
    R = rxy_interp.eval({x,y,z});
    err_grid_vec(mm_ind) = norm(F(:)-R(:))/norm(F(:));
    xtest = linspace(-1,1.,300);
    ytest = linspace(-1,1.,300);
    ztest = linspace(-1,1.,300);
    

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


%% isolate tol = 10^-8, compare to p-AAA
clear
run init.m
nx=151;
ny=151;
nz=151;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);

del = 1/64;
f = @(x,y,z) 1./sqrt(x.^2+2*y.^2+3*z.^2+del);

[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));


options.tolAAA = 1e-9;
options.tol_qr = 1e-10;
[rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);

xtest = linspace(-1,1.,351);
ytest = linspace(-1,1.,351);
ztest = linspace(-1,1.,351);


[Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
Fexact = f(Xtest,Ytest,Ztest);
Finterp =  rxy_interp.eval({xtest,ytest,ztest});
norm(Finterp(:)-Fexact(:),inf)/norm(Fexact(:),inf)

%% Compare to two-step
clear
run init.m
nx=151;
ny=151;
nz=151;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);

del = 1/64;
f = @(x,y,z) 1./sqrt(x.^2+2*y.^2+3*z.^2+del);

[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));


options.tolAAA = 1e-8;
options.mmax={15,15,15};
options.tol_qr = 1e-9;
options.twostep = true;
options.paaapts = {x(1:10:end),y(1:10:end),z(1:10:end)};
options.tolpaaa = 1e-8;
[rxy_interp,info_interp] = construct_multi_AAA({x,y,z},F,f,options);

xtest = linspace(-1,1.,351);
ytest = linspace(-1,1.,351);
ztest = linspace(-1,1.,351);


[Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
Fexact = f(Xtest,Ytest,Ztest);
Finterp =  rxy_interp.eval({xtest,ytest,ztest});
norm(Finterp(:)-Fexact(:),inf)/norm(Fexact(:),inf)
%% compare to lr_pAAA (pAAA crashes)
clear
nx=51;
ny=51;
nz=51;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);

del = 1/64;
f = @(x,y,z) 1./sqrt(x.^2+2*y.^2+3*z.^2+del);

[X,Y,Z] = ndgrid(x,y,z);
F = tensor(f(X,Y,Z));


options.tolpAAA = 1e-8;

[rxy_pAAA,info_pAAA] = lr_paaa(F.data,{x,y,z},options.tolpAAA,7);

xtest = linspace(-1,1.,351);
ytest = linspace(-1,1.,351);
ztest = linspace(-1,1.,351);


[Xtest,Ytest,Ztest] = ndgrid(xtest,ytest,ztest);
Fexact = f(Xtest,Ytest,Ztest);
Fpaaa =  rxy_pAAA.eval({xtest,ytest,ztest});
norm(Fpaaa-Fexact,'fro')/norm(Fexact,"fro")


