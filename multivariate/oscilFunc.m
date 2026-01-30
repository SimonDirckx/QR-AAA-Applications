%% sep. rat. approx to a non-rational function
clear
run init.m
clc

nx=51;
ny=51;
nz=51;
nt=51;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);
t = linspace(-1,1,nt);

mmax_vec = 5:5:50;
err_off_grid_vec=zeros(numel(mmax_vec),1);
err_grid_vec=zeros(numel(mmax_vec),1);
info_vec = cell(1,numel(mmax_vec));

f = @(x,y,z,t) (cos(3*pi*x).*sin(4*pi*y)+sin(2*pi*x).*sin(5*pi*y))./(x.*(y.^2)+(z.^3).*t+3);

[X,Y,Z,T] = ndgrid(x,y,z,t);
F = tensor(f(X,Y,Z,T));

for mm_ind=1:numel(mmax_vec)
    mmax = mmax_vec(mm_ind);
    disp("mmax = ")
    mmax
    options.mmax=mmax;
    options.sep=true;
    options.tolAAA = 1e-12;
    options.tol_qr = 1e-15;
    [rxy_interp,info_interp] = construct_multi_AAA({x,y,z,t},F,f,options);
    info_vec{mm_ind} = info_interp;
    R = rxy_interp.eval({x,y,z,t});
    err_grid_vec(mm_ind) = norm(F(:)-R(:))/norm(F(:));
    xtest = linspace(-1,1.,81);
    ytest = linspace(-1,1.,81);
    ztest = linspace(-1,1.,81);
    ttest = linspace(-1,1.,81);

    [Xtest,Ytest,Ztest,Ttest] = ndgrid(xtest,ytest,ztest,ttest);
    Fexact = f(Xtest,Ytest,Ztest,Ttest);
    Finterp =  rxy_interp.eval({xtest,ytest,ztest,ttest});
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

%% isolate tol = 10^-8
clear
run init.m
clc

nx=51;
ny=51;
nz=51;
nt=51;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);
t = linspace(-1,1,nt);


f = @(x,y,z,t) (cos(3*pi*x).*sin(4*pi*y)+sin(2*pi*z).*sin(5*pi*t))./(x.*(y.^2)+(z.^3).*t+4);
[X,Y,Z,T] = ndgrid(x,y,z,t);
F = tensor(f(X,Y,Z,T));


options.tolAAA = 1e-8;
options.tol_qr = 1e-9; 
[rxy_interp,info_interp] = construct_multi_AAA({x,y,z,t},F,f,options);

xtest = linspace(-1,1.,81);
ytest = linspace(-1,1.,81);
ztest = linspace(-1,1.,81);
ttest = linspace(-1,1.,81);

[Xtest,Ytest,Ztest,Ttest] = ndgrid(xtest,ytest,ztest,ttest);
Fexact = f(Xtest,Ytest,Ztest,Ttest);
Finterp =  rxy_interp.eval({xtest,ytest,ztest,ttest});
norm(Finterp-Fexact,'fro')/norm(Fexact,"fro")


%% compare to paaa

clear
run init.m
% lr_paaa crashes if grid is too fine
nx=21;
ny=21;
nz=21;
nt=21;

x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);
t = linspace(-1,1,nt);
f = @(x,y,z,t) (cos(3*pi*x).*sin(4*pi*y)+sin(2*pi*z).*sin(5*pi*t))./(x.*(y.^2)+(z.^3).*t+4);
[X,Y,Z,T] = ndgrid(x,y,z,t);
F = tensor(f(X,Y,Z,T));

optionspaaa.max_iter = 50;
tolpaaa = 1e-8;
rkpaaa = 7;
tic
[rxy_pAAA,info_pAAA] = lr_paaa(F.data,{x,y,z,t},tolpaaa,rkpaaa,optionspaaa);
toc
xtest = linspace(-1,1.,81);
ytest = linspace(-1,1.,81);
ztest = linspace(-1,1.,81);
ttest = linspace(-1,1.,81);


[Xtest,Ytest,Ztest,Ttest] = ndgrid(xtest,ytest,ztest,ttest);
Fexact = f(Xtest,Ytest,Ztest,Ttest);
Fpaaa =  rxy_pAAA.eval({xtest,ytest,ztest,ttest});
norm(Fpaaa(:)-Fexact(:),2)/norm(Fexact(:),2)