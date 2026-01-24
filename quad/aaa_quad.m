%% AAA_QUAD_sing_log
% This script tests QR-AAA based quadrature on functions of the form
% x^alpha+x^beta log(x)
% for alpha \in (0,30], beta \in [0,1]
run init.m
clear
clc
fs = @(x,alph) x.^alph;
int_class = @(alph) 1./(1+alph);
use_exact_int = true;
%use_exact_int = false;

tolvec = [1e-5,1e-6,1e-7,1e-8,1e-9,1e-10];
err1 = zeros(numel(tolvec),1);
err2 = zeros(numel(tolvec),1);
err3 = zeros(numel(tolvec),1);

W=zeros(100,numel(tolvec));
wghts_sz_vec = zeros(numel(tolvec),1);

for indtol=1:numel(tolvec)
    tol = tolvec(indtol);
p=5000;
xpts = chebpts(p).';
xpts = (1+xpts)/2.;
xref = logspace(-8,-1,50);%avoid quadGK instabilities poles, xref(1) not TOO small
xpts  = [xpts,xref];

xpts = unique(xpts);
xpts = sort(xpts).';

alph = linspace(0,50,3000);

Fmat = [fs(xpts,alph)];

IF = int_class(alph);

for j= 1:size(Fmat,2)
    if max(abs(Fmat(:,j)))>0 && max(abs(Fmat(:,j)))<Inf
        mx = max(abs(Fmat(:,j)));
        Fmat(:,j)=Fmat(:,j)/mx;
        IF(j) = IF(j)/mx;
    end
end
N = numel(xpts);
tic
[Q,R,P] = qr(Fmat,"econ");
toc
if use_exact_int
    tol_qr = 1e-12;
else
    tol_qr = 1e-12;
end
k = sum(abs(diag(R))>tol_qr*abs(R(1,1)));
Qk = Q(:,1:k);
Rk = R(1:k,:);
Gamm=diag(diag(R(1:k,1:k)));
Fk_qr = Q(:,1:k)*diag(diag(R(1:k,1:k)));
R = diag(diag(R(1:k,1:k)))\R(1:k,:);



tic
[r_qr, pol_qr, res_qr, zer_qr, z_qr, ff_qr, w_qr, errvec_qr] = aaa_sv(Fk_qr,xpts,'tol',tol);

%remove spurious poles
ii = find(real(pol_qr)>-1e-14 & real(pol_qr)<1 & abs(imag(pol_qr))<1e-14);
[r_qr,pol_qr,res_qr,zer_qr]=cleanup(z_qr,ff_qr,xpts,Fk_qr,pol_qr,ii);
toc





dlmwrite('quad_poles.dat',[real(pol_qr),imag(pol_qr)]);
I_qr = zeros(numel(w_qr),1);
ind_qr = [];
for i=1:numel(z_qr)
    ind_qr=[ind_qr,find(xpts==z_qr(i))];
end

if use_exact_int
    Fm = Fmat(ind_qr,:);
    [U,S,V] = svd(Fm);
    s = diag(S);
    del = max(1e-15,tol/N);
    k = sum(s>del*s(1));

    int_wghts_test = IF*V(:,1:k)*(S(1:k,1:k)\U(:,1:k)');
    int_wghts_test=int_wghts_test.';
    int_wghts_qr = int_wghts_test;
    ii = find(abs(int_wghts_qr)<1e-15);
    z_qr(ii) = [];
    int_wghts_qr(ii) = [];
else
    for j=1:numel(w_qr)
        ff00 = zeros(numel(w_qr),1);
        ff00(j) = 1.;
        r00 = @(x) reval(x, z_qr, ff00, w_qr);
        I_1 = quadgk(r00,0,1e-7,'RelTol',1e-13);
        I_2 = quadgk(r00,1e-7,1,'RelTol',1e-13);
        I_qr(j) = I_1+I_2;
    end
    int_wghts_qr = I_qr;
    ii = find(abs(int_wghts_qr)<1e-15);
    z_qr(ii) = [];
    int_wghts_qr(ii) = [];
end
disp('order is')
numel(int_wghts_qr)
W(1:numel(int_wghts_qr),indtol) = int_wghts_qr;
wghts_sz_vec(indtol) = numel(int_wghts_qr);

%-----------------------------------------------------


err1(indtol) = test_accuracy('rootsin',z_qr,int_wghts_qr,false);
err2(indtol) = test_accuracy('rootcos',z_qr,int_wghts_qr,false);
err3(indtol) = test_accuracy('besselj',z_qr,int_wghts_qr,false);
end
%%
figure(1);clf;
loglog(tolvec.',err1)
hold on
loglog(tolvec.',err2)
loglog(tolvec.',err3)
loglog(tolvec.',tolvec,'--')
legend('rootsin','rootcos','besselj')

sum_wghts = ones(1,size(W,1))*abs(W);
if use_exact_int
    dlmwrite('quad_errors_use_exact_int.dat',[tolvec.',err1,err2,err3]);
    dlmwrite('quad_sum_wghts_use_exact_int.dat',[tolvec.',sum_wghts.']);
    dlmwrite('quad_supp_exact_int.dat',z_qr);
else
    dlmwrite('quad_errors.dat',[tolvec.',err1,err2,err3]);
    dlmwrite('quad_sum_wghts.dat',[tolvec.',sum_wghts.']);
    dlmwrite('quad_sing_supp.dat',z_qr);
end


figure(2);clf;
loglog(tolvec,sum_wghts)


function err_max=test_accuracy(type,z_qr,int_wghts_qr,bool_plot)
aa_vec = linspace(0,10,1000);  
err = zeros(numel(aa_vec),1);

switch type
    case 'rootsin'
        f = @(x,alpha) (x.^(alpha)).*sin(x);
        If = @(x,a) x.^(a+1).*(log(x)/(a+1)-(1/(a+1)^2));
        Iexact = @(a) quadgk(@(x)f(x,a),0,1,'RelTol',1e-13);
    case 'besselj'
        f = @(x,a) besselj(0,a*x);
        PHIj = @(x)(pi*x/2).*(besselj(1,x).*struve(0,x)-besselj(0,x).*struve(1,x));
        intj = @(x) x.*besselj(0,x)+PHIj(x);
        Iexact = @(a)intj(a)/a;
    case 'rootcos'
        f = @(x,alpha) (x.^(alpha)).*cos(x);
        If = @(x,a) x.^(a+1).*(log(x)/(a+1)-(1/(a+1)^2));
        Iexact = @(a) quadgk(@(x)f(x,a),0,1,'RelTol',1e-13);
end

for i =1:numel(aa_vec)
I = sum(int_wghts_qr.*f(z_qr,aa_vec(i)));
err(i) = abs(I-Iexact(aa_vec(i)) );
end


if bool_plot
    figure(101);clf;
    semilogy(err)
end

err_max = max(err);


end

%% test accuracy
%aa_vec = linspace(0,alph(end),1000);
aa_vec = linspace(0,10,1000);
err_qr = zeros(numel(aa_vec),1);
bnd = zeros(numel(aa_vec),1);
f = @(x,alpha) (x.^(alpha)).*log(x)+2*x.^alpha+cos(2*x);
If1 = @(x,a) x.^(a+1).*(log(x)/(a+1)-(1/(a+1)^2));
If2 = @(x,a) 2*(x.^(a+1))/(a+1);
If3 = @(x) sin(2*x)/2;



for i =1:numel(aa_vec)
bnd(i)  = 1e-10*(sum(abs(Fk_qr(ind_qr,:)\f(xpts(ind_qr),aa_vec(i)))));
I = If1(1,aa_vec(i))+If2(1,aa_vec(i))+If3(1);
Ihat_qr = sum(int_wghts_qr.*f(z_qr,aa_vec(i)));
err_qr(i) = abs(I-Ihat_qr)/abs(I);
end
figure(1)
clf
semilogy(aa_vec,err_qr)
dlmwrite('quad_sing_err_logfunc.dat',[aa_vec.',err_qr]);

%% test on Bessel
clc
khvec = linspace(.1,10,5000);

PHIy = @(x)(pi*x/2).*(bessely(1,x).*struve(0,x)-bessely(0,x).*struve(1,x));
inty = @(x) x.*bessely(0,x)+PHIy(x);

PHIj = @(x)(pi*x/2).*(besselj(1,x).*struve(0,x)-besselj(0,x).*struve(1,x));
intj = @(x) x.*besselj(0,x)+PHIj(x);



erry = zeros(numel(khvec),1);
errj = zeros(numel(khvec),1);
for i=1:numel(khvec)
y = @(x) bessely(0,khvec(i)*x);%cos_lim(x,-.31138,5.16552);
j = @(x) besselj(0,khvec(i)*x);


Iexacty = inty(khvec(i))/khvec(i);
Iexactj = intj(khvec(i))/khvec(i);
%tic
%Ireg = integral(f,0,1);
%toc
%tic
Igky = Iexacty;%quadgk(y,0,1,'RelTol',1e-13);
Igkj = Iexactj;%quadgk(j,0,1,'RelTol',1e-13);
%toc

%tic
Ihat_qry = sum(int_wghts_qr.*y(z_qr));
Ihat_qrj = sum(int_wghts_qr.*j(z_qr));
%toc

erry(i) = abs(Igky-Ihat_qry);
errj(i) = abs(Igkj-Ihat_qrj);
%abs(Iexact-Ireg)/abs(Iexact)
%abs(Iexact-Igk)/abs(Iexact)
end
figure(2)
clf
semilogy(khvec,erry)
hold on
semilogy(khvec,errj)
legend('erry','errj')
dlmwrite('quad_sing_err_besselfunc.dat',[khvec.',errj,erry]);


%% test on strong singularity
clc
alphavec = linspace(0,10,5000);

err = zeros(numel(alphavec),1);
for i=1:numel(alphavec)
f = @(x) cos(alphavec(i)*x).*log(x)+x.^(alphavec(i));

Iexact = quadgk(f,0,1,'RelTol',1e-13);%(-sinint(alphavec(i))+sinint(0))/alphavec(i);
Ihat = sum(int_wghts_qr.*f(z_qr));
err(i) = abs(Iexact-Ihat);
end
figure(2)
clf
semilogy(alphavec,err)
%dlmwrite('quad_sing_err_besselfunc.dat',[khvec.',errj,erry]);



%% Evaluate rational function in barycentric form.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REVAL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = reval(zz, zj, fj, wj)
%   Construct function handle to evaluate rational function in barycentric form.

zv = zz(:);                         % vectorize zz if necessary
CC = 1./(zv-zj.');                  % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);         % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL.

function [pol, res, zer] = prz(r, zj, fj, wj)
% Compute poles, residues, and zeros of rational function in barycentric form.
m = length(wj);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

% Compute residues via discretized Cauchy integral:
dz = 1e-5*exp(2i*pi*(1:4)/4);

pp = pol+dz;
rvals = r(pp(:));
res = zeros(length(pol),size(rvals,2));
for it = 1:size(rvals,2)
    res(:,it) = reshape(rvals(:,it),[],4)*dz.'/4;
end

% Compute zeros via generalized eigenvalue problem:
for it = 1:size(fj,2)
    E = [0 (wj.*fj(:,it)).'; ones(m, 1) diag(zj)];
    zer{it} = eig(E, B);
    % Remove zeros of numerator at infinity:
    zer{it} = zer{it}(~isinf(zer{it}));
end
end % End of PRZ().


function I = int_from_pzr(zj, fj, wj,zz,func)

[pol, ~, ~] = prz(zj, fj, wj);

C = pol.'./((pol.')-zz);
C = [ones(size(C,1),1),C];
y = func(zz);
resvec = C\y;

I = resvec(1);
for i = 2:numel(resvec)
    I = I + resvec(i)*pol(i-1)*( log( abs( 1-pol(i-1) ) )-log( abs( pol(i-1) ) ) );
end

end

function b=spurious_pole(pol)
    b = any(real(pol)>0 & real(pol)<1 & abs(imag(pol))<1e-14);
end

function [r,pol,res,zer]=cleanup(z,f,Z,F,pol,ii)
ni=numel(ii);
% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj,:) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    F(Z == z(jj),:) = [];
    Z(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
% A: slow build form now
A=[];
for i=1:size(F,2)
    Fi = F(:,i);
    fi = f(:,i);
    Asub = Fi.*C-C.*fi.';
    A=[A;Asub];
end


[~, S, V] = svd(A, 0);
s=diag(S);
ind0 = find(s==s(end));
w = V(:,ind0)*ones(numel(ind0),1)/sqrt(numel(ind0));

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval_sv(zz, z, f, w);
[pol, res, zer] = prz(r,z, f, w);
end
function r = reval_sv(zz, zj, fj, wj)
% Evaluate rational function in barycentric form.
l = length(zz);
zv = zz(:);                             % vectorize zz if necessary
CC = 1./(zv-zj.');   % Cauchy matrix
r = CC*(wj.*fj)./(CC*wj);

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv),:) = kron(ones(sum(isinf(zv)),1),sum(wj.*fj,1)./sum(wj));

% Deal with NaN:
ii = find(isnan(r));
[row,col] = ind2sub(size(r),ii);
ii = [row,col];
ii(ii(:,1) == 0) = l;
for jj = 1:size(ii,1)
    if ( isnan(zv(ii(jj,1))) || ~any(zv(ii(jj,1)) == zj) )
        disp('REACHED')
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj,1),ii(jj,2)) = fj(zv(ii(jj,1)) == zj,ii(jj,2));
    end
end

% Reshape to input format:
% r = reshape(r, length(zz),size(fj,2));

end % End of REVAL().
