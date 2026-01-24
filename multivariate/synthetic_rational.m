%% two-step rational approximation for 'adversarial example': order 50 bivariate system
clear
run init.m
clc

nx=500;
ny=50;

x = logspace(0,4,nx)*1i;
y = logspace(-1.5,0,ny);

% from pAAA

% warning: set N = 50 in 'evaluate_synthetic_parametric.m', otherwise grids 
% and initial degrees must be adapted and experiments take VERY LONG

f = @(x,y) evaluate_synthetic_parametric(x,y);

F = tensor(f(x,y));



%% separable qr-AAA approach, interpolative -> fails!
clear options
options.interp=true;
options.twostep=false;
options.tolAAA = {1e-13,1e-13};
options.tol_qr = {1e-15,1e-15};
options.valpts = {x,y};
tic
[rxy_interp,info_interp] = construct_multi_AAA({x,y},F,f,options);
toc

%% test separable.
xtest = logspace(0,4,500)*1i;
ytest = logspace(-1.5,0,100);
Fexact = f(xtest,ytest);
Finterp =  rxy_interp.eval({xtest,ytest});
norm(Finterp-Fexact,'fro')/norm(Fexact,"fro")
norm(Finterp(:)-Fexact(:),inf)/norm(Fexact(:),inf)


%% two-step approach

clear options

options.interp=true;
options.twostep=true;
options.tolAAA = {1e-8,1e-5};
options.mmax = {30,100};
options.tol_qr = {1e-15,1e-15};
% toggle for convergence information of paaa
%options.valpts = {logspace(0,4,2*nx)*1i,logspace(-1.5,0,2*ny)};
%options.Fvalpts = tensor(f(options.valpts{:}));
tic
[rxy_2step,info_2step] = construct_multi_AAA({x,y},F,f,options);
toc

%% test two-step
xtest = logspace(0,4,1000);
ytest = logspace(-1.5,0,200);
Fexact = f(xtest,ytest);
F2step =  rxy_2step.eval({xtest,ytest});
norm(F2step-Fexact,'fro')/norm(Fexact,"fro")
norm(F2step(:)-Fexact(:),inf)/norm(Fexact(:),inf)



%% pAAA (for comparison)
clear options
options.max_iter = 100;
valpts = {logspace(0,4,2*nx)*1i,logspace(-1.5,0,2*ny)};
Fvalpts = f(valpts{:});
options.validation.sampling_values = valpts;
options.validation.samples = Fvalpts;

tic
[rxy_pAAA,info_pAAA] = paaa(F.data,{x,y},1e-10,options);
toc

%% test pAAA


xtest = logspace(0,4,1000);
ytest = logspace(-1.5,0,200);


Fexact = f(xtest,ytest);

Fpaaa =  rxy_pAAA.eval({xtest,ytest});

norm(Fpaaa-Fexact,'fro')/norm(Fexact,"fro")
norm(Fpaaa(:)-Fexact(:),inf)/norm(Fexact(:),inf)
