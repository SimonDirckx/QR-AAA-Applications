%% AAA_QUAD_sing_compare
% This script tests QR-AAA based quadrature on functions of the form
% x^alpha+x^beta log(x)
% for alpha \in (0,30], beta \in [0,1]
% and compares to the same approximation to the integral, computed using
% sketchAAA
run init.m
clear
clc
fs = @(x,alph) x.^alph;
fl = @(x,alph) log(x).*x.^alph;
fl2 = @(x,alph) log(alph.*x);


p=1000;
xpts = legpts(p).';
xpts = (1+xpts)/2.;
xref = logspace(-8,-1,50);
xpts  = [xpts,xref];

xpts = unique(xpts);
xpts = sort(xpts).';

alph = linspace(0,30,500);
beta = logspace(-14,-2,100);
gamm = linspace(.5,10,500);

Fmat = [fl2(xpts,gamm),fs(xpts,alph),fl(xpts,beta)];%,fc(xpts.',alph_beta)];
for j= 1:size(Fmat,2)
    if max(abs(Fmat(:,j)))>0 && max(abs(Fmat(:,j)))<Inf
        Fmat(:,j)=Fmat(:,j)/max(abs(Fmat(:,j)));
    end
end
tic
[Q,R,P] = qr(Fmat,"econ");
toc
tol = 1e-8;
k = sum(abs(diag(R))>tol*abs(R(1,1)));
Fk_qr = Q(:,1:k)*diag(diag(R(1:k,1:k)));
R = diag(diag(R(1:k,1:k)))\R(1:k,:);



tic
[r_qr, pol_qr, res_qr, zer_qr, z_qr, ff_qr, w_qr, errvec_qr] = aaa_sv(Fk_qr,xpts,'tol',1e-10);
toc
%norm(r_qr(xpts)-Fk)/norm(Fk)
%norm(r_qr(xpts)*R/P-Fmat)/norm(Fmat)


Fk = Fmat*randn(size(Fmat,2),20);
opts.tol = 1e-10; 
[z_sketch,w_sketch,ind_sketch,stats_sketch,stats_T_sketch] = sketchAAA(Fk.',xpts,100,opts);


I_qr = zeros(numel(w_qr),1);
I_sketch = zeros(numel(w_sketch),1);

for j=1:numel(w_qr)
    ff00 = zeros(numel(w_qr),1);
    ff00(j) = 1.;
    r00 = @(x) reval(x, z_qr, ff00, w_qr);
    I_1 = quadgk(r00,0,1e-8,'RelTol',1e-13);
    I_2 = quadgk(r00,1e-8,1,'RelTol',1e-13);
    I_qr(j) = I_1+I_2;%quadgk(r00,0,1,'RelTol',1e-15);

end


for j=1:numel(w_sketch)
    ff00 = zeros(numel(w_sketch),1);
    ff00(j) = 1.;
    r00 = @(x) reval(x, z_sketch, ff00, w_sketch);
    I_1 = quadgk(r00,0,1e-8,'RelTol',1e-13);
    I_2 = quadgk(r00,1e-8,1,'RelTol',1e-13);
    I_sketch(j) = I_1+I_2;%integral(r00,0,1,'RelTol',1e-13);
end

%-----------------------------------------------------
int_wghts_qr = I_qr;
int_wghts_sketch = I_sketch;



%% test log
f = @(x) log(x)+1;
I = integral(f,0,1);
Ihat_qr = sum(int_wghts_qr.*f(z_qr));
Ihat_sketch = sum(int_wghts_sketch.*f(z_sketch));

%% test accuracy
%aa_vec = linspace(0,alph(end),1000);
aa_vec = linspace(0,1,10000);
err_qr = zeros(numel(aa_vec),1);
err_sketch = zeros(numel(aa_vec),1);
f = @(x,alpha) (x.^(alpha)).*log(x)+2*x.^alpha+cos(2*x);
If1 = @(x,a) x.^(a+1).*(log(x)/(a+1)-(1/(a+1)^2));
If2 = @(x,a) 2*(x.^(a+1))/(a+1);
If3 = @(x) sin(2*x)/2;


for i =1:numel(aa_vec)
I = If1(1,aa_vec(i))+If2(1,aa_vec(i))+If3(1);
Ihat_qr = sum(int_wghts_qr.*f(z_qr,aa_vec(i)));
Ihat_sketch = sum(int_wghts_sketch.*f(z_sketch,aa_vec(i)));
err_qr(i) = abs(I-Ihat_qr)/abs(I);
err_sketch(i) = abs(I-Ihat_sketch)/abs(I);
end
figure(1)
clf
semilogy(aa_vec,err_qr)
hold on
semilogy(aa_vec,err_sketch)
legend('qr','sketch')




%% test on Bessel
f = @(x) bessely(0,x);%cos_lim(x,-.31138,5.16552);

PHIy = @(x)(pi*x/2).*(bessely(1,x).*struve(0,x)-bessely(0,x).*struve(1,x));
inty = @(x) x.*bessely(0,x)+PHIy(x);
Iexact = inty(1);

I = integral(f,0,1);
Ihat_qr = sum(int_wghts_qr.*f(z_qr));
Ihat_sketch = sum(int_wghts_sketch.*f(z_sketch));
abs(I-Ihat_qr)/abs(I)
abs(I-Ihat_sketch)/abs(I)






%% compute error bands
N=100;
perc = 0;
err = zeros(numel(aa_vec),N);

f = @(x,alpha) (x.^(alpha)).*log(x)+2*x.^alpha;
If1 = @(x,a) x.^(a+1).*(log(x)./(a+1)-(1./(a+1).^2));
If2 = @(x,a) 2*(x.^(a+1))./(a+1);


for ind=1:N
    disp(ind)
Fk = Fmat*randn(size(Fmat,2),10);
opts.tol = 1e-10;
[z_sketch,w_sketch,ind_sketch,stats_sketch,stats_T_sketch] = sketchAAA(Fk.',xpts,100,opts);


I_sketch = zeros(numel(w_sketch),1);

for j=1:numel(w_sketch)
    ff00 = zeros(numel(w_sketch),1);
    ff00(j) = 1.;
    r00 = @(x) reval(x, z_sketch, ff00, w_sketch);
    I_1 = quadgk(r00,0,1e-7,'RelTol',1e-13);
    I_2 = quadgk(r00,1e-7,1,'RelTol',1e-13);
    I_sketch(j) = I_1+I_2;%integral(r00,0,1,'RelTol',1e-10);
end


int_wghts_sketch = I_sketch;



I = If1(1,aa_vec)+If2(1,aa_vec);
Ihat_sketch = sum(int_wghts_sketch.*f(z_sketch,aa_vec));
err(:,ind) = abs(I-Ihat_sketch)./abs(I);

end
freq=sum(perc)/N

%%
figure(1)
clf
semilogy(aa_vec,err,'Color', [.6, .6, .6])
hold on
semilogy(aa_vec,err_qr,'Color',[1.,0.,0.],'LineWidth',2)



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

function [pol, res, zer] = prz(zj, fj, wj)
%   Compute poles, residues, and zeros of rational fun in barycentric form.

% Compute poles via generalized eigenvalue problem:
m = length(wj);
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
pol = pol(~isinf(pol));

% Compute residues via formula for res of quotient of analytic functions:
% (This is not always accurate and starting in 2025 is overwritten in
% the calling function by a least-squares computation of residues.)
N = @(t) (1./(t-zj.')) * (fj.*wj);
Ddiff = @(t) -((1./(t-zj.')).^2) * wj;
res = N(pol)./Ddiff(pol);

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
zer = zer(~isinf(zer));

end % End of PRZ.


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


