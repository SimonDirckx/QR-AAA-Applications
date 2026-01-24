%% new L-shape script
% uses qr-aaa and trig interp

clear
run init.m
f = @(x,y)sin(2.*pi*x).*sin(2.*pi*y);
f_polar = @(r,th) f(real(gamm(r,th)),imag(gamm(r,th)));
rmin = 0;
th = chebpts(150);
th = (1+th)/2;
th1 = 2*th;
th2 = 2+th;
th3 = 3+th;
th4 = 4+th;
th5 = 5+th;
th6 = 6+2*th;
th = [th1;th2;th3;4+th;5+th;6+2*th];
%%
rvec=linspace(rmin,1,1000).';
nt = numel(th);
nr = numel(rvec);


Fmat = f_polar(rvec,th.');
F=Fmat.';
F=F(:);
GAMM = gamm(rvec,th.').';
GAMM = GAMM(:);
bdry = gamm(1,th.');

g = gamm(.5,th);

%%


rho = zeros(numel(rvec)*numel(th),2);
rho(:,1) = kron(rvec,ones(size(th)));
rho(:,2) = kron(ones(size(rvec)),th);
XY = rho(:,1).*[cos(rho(:,2)),sin(rho(:,2))];

mn = min(real(GAMM(:)));
mx = max(real(GAMM(:)));
xy = linspace(mn-.5,mx+.5,500);
[X,Y] = meshgrid(xy,xy);
FF = f(X,Y);
figure(1)
clf
hold on
imagesc(xy,xy,FF)
plot(real(bdry),imag(bdry),'k','LineWidth',4)
axis('equal')
axis('off')
%%
degr = 100;

PHI1 = build_PHI(degr,th1);
PHI2 = build_PHI(degr,th2);
PHI3 = build_PHI(degr,th3);
PHI4 = build_PHI(degr,th4);
PHI5 = build_PHI(degr,th5);
PHI6 = build_PHI(degr,th6);


PHI = blkdiag(PHI1,PHI2,PHI3,PHI4,PHI5,PHI6);



Ftild = (PHI\Fmat')';

%%
tol_qr = 1e-16;
tic
[Q,R,P] = qr(Ftild,'econ','vector');

k_qr = sum(abs(diag(R))>tol_qr*abs(R(1,1)));
Fk = Q(:,1:k_qr)*diag(diag(R(1:k_qr,1:k_qr)));
[r, pol, res, zer, z, ff, w, errvec] = aaa_sv(Fk,rvec,'tol',1e-13);

ind_supp = [];
for i=1:numel(z)
    ind_supp=[ind_supp,find(rvec==z(i))];
end
r = @(zz) reval_sv(zz, z, Ftild(ind_supp,:), w);




%%




nr_ext = 2*nr;
rext = linspace(0,3,nr_ext).';
thext = chebpts(200);
thext = (1+thext)/2;
th1ext = 2*thext;
th2ext = 2+thext;
th3ext = 3+thext;
th4ext = 4+thext;
th5ext = 5+thext;
th6ext = 6+2*thext;
thext = [th1ext;th2ext;th3ext;th4ext;th5ext;th6ext];
%thext = th;%(1+thext)/2;

PHI1ext = build_PHI(degr,th1ext);
PHI2ext = build_PHI(degr,th2ext);
PHI3ext = build_PHI(degr,th3ext);
PHI4ext = build_PHI(degr,th4ext);
PHI5ext = build_PHI(degr,th5ext);
PHI6ext = build_PHI(degr,th6ext);

PHI = blkdiag(PHI1ext,PHI2ext,PHI3ext,PHI4ext,PHI5ext,PHI6ext);



Rtest = (r(rext));
Rtest = PHI*(Rtest.');
Fext=Rtest(:);
Fexact = f_polar(rext,thext.').';




GAMMext = gamm(rext,thext.');
GAMMext = GAMMext.';
GAMMext = GAMMext(:);

figure(3)
clf
scatter(real(GAMMext),imag(GAMMext),[],log10(abs(Fext-Fexact(:))),'filled')
hold on
plot(real(bdry),imag(bdry),'k')
colorbar
axis('equal')

figure(4)
clf
scatter(real(GAMMext),imag(GAMMext),[],Fext,'filled')
hold on
plot(real(bdry),imag(bdry),'k')
colorbar
axis('equal')

%%
addpath ../chebfun-master
th0 = thext(100); %select th not in approx domain
rmin_aaa = 0;
rext = linspace(rmin,3,1000).';
rvec_aaa = linspace(rmin_aaa,1,nr).';
gamm_vec_aaa = gamm(rvec_aaa,th0);
fvec_aaa = f_polar(rvec_aaa,th0);
[r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(fvec_aaa,rvec_aaa,'tol',1e-13);




PHIvec = PHI(100,:);

fext_qr = r(rext);
fext_qr = PHIvec*(fext_qr');


gamm_ext = gamm(rext,th0);
fext_aaa = r_aaa(rext);


%%
figure(12);clf;
plot(abs(gamm_ext),fext_qr)
hold on
plot(abs(gamm_ext),fext_aaa)
plot(abs(gamm_ext),f(real(gamm_ext),imag(gamm_ext)))
legend('qr','aaa','ref')


fext_qr = r(rext);
fext_qr = (PHIvec*(fext_qr'))';%what is happening here???


gamm_ext = gamm(rext,th0);
fext_aaa = r_aaa(rext);
f_exact = f(real(gamm_ext),imag(gamm_ext));


loss_ind_aaa=0;
loss_ind_qr=0;
err = 0;
while err<1e-6 && loss_ind_aaa<numel(fext_aaa)
    loss_ind_aaa=loss_ind_aaa+1;
    err = abs(f_exact(loss_ind_aaa)-fext_aaa(loss_ind_aaa));
end
err=0;
while err<1e-6 && loss_ind_qr<numel(fext_qr)
    loss_ind_qr=loss_ind_qr+1;
    err = abs(f_exact(loss_ind_qr)-fext_qr(loss_ind_qr));
end


loss_pt_aaa = abs(gamm_ext(loss_ind_aaa));
loss_pt_qr = abs(gamm_ext(loss_ind_qr));
bdry = abs(gamm(1.,th0)); ycut = [1e-20 10*max(abs(f_exact-fext_aaa))];

figure(13);clf
colors = get(gca, 'colororder');
semilogy(abs(gamm_ext),abs(f_exact-fext_qr),'Color',colors(1,:))
hold on
semilogy(abs(gamm_ext),abs(f_exact-fext_aaa),'Color',colors(2,:))

semilogy([bdry bdry],ycut,'k--',LineWidth=2)
semilogy([loss_pt_aaa loss_pt_aaa],ycut,'Color',colors(2,:),'LineStyle','--',LineWidth=2)
semilogy([loss_pt_qr loss_pt_qr],ycut,'Color',colors(1,:),'LineStyle','--',LineWidth=2)
legend('qr','aaa','a','b','c')
%%
loss_prec =1e-8;
loss_pts_qr = zeros(numel(thext),1);
loss_pts_aaa = zeros(numel(thext),1);
rext = linspace(rmin,4,1000).';
mx_n_pol = 0;
avg_n_pol = 0;

for ind=1:numel(thext)
    th0 = thext(ind);
    rmin_aaa = 0;
    rvec_aaa = linspace(rmin_aaa,1,nr).';
    gamm_vec_aaa = gamm(rvec_aaa,th0);
    fvec_aaa = f(real(gamm_vec_aaa),imag(gamm_vec_aaa));
    [r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(fvec_aaa,rvec_aaa,'tol',1e-13);
    mx_n_pol = max(numel(pol_aaa),mx_n_pol);
    avg_n_pol = avg_n_pol+numel(pol_aaa);
    PHIvec = PHI(ind,:);
    
    fext_qr = r(rext);
    fext_qr = (PHIvec*(fext_qr.')).';
    
    
    gamm_ext = gamm(rext,th0);
    fext_aaa = r_aaa(rext);
    f_exact = f(real(gamm_ext),imag(gamm_ext));
    loss_ind_aaa=0;
    loss_ind_qr=0;
    err = 0;
    while err<loss_prec && loss_ind_aaa<numel(f_exact)
        loss_ind_aaa=loss_ind_aaa+1;
        err = abs(f_exact(loss_ind_aaa)-fext_aaa(loss_ind_aaa));
    end
    err=0;
    while err<loss_prec && loss_ind_qr<numel(f_exact)
        loss_ind_qr=loss_ind_qr+1;
        err = abs(f_exact(loss_ind_qr)-fext_qr(loss_ind_qr));
    end
    
    
    loss_pts_aaa(ind) = min(abs(rext(loss_ind_aaa)),4);
    loss_pts_qr(ind) = min(abs(rext(loss_ind_qr)),4);

end
avg_n_pol = avg_n_pol/numel(thext);
%%
figure(14);clf
colors = get(gca,"ColorOrder");
plot(real(gamm(loss_pts_aaa,thext)),imag(gamm(loss_pts_aaa,thext)),'LineStyle','-','LineWidth',2)
hold on
plot(real(gamm(loss_pts_qr,thext)),imag(gamm(loss_pts_qr,thext)),'LineStyle','-','LineWidth',2)


mn_aaa = min(abs(gamm(loss_pts_aaa,thext)));
mn_qr = min(abs(gamm(loss_pts_qr,thext)));
bdry = gamm(1.,thext);
bdry_qr = mn_qr*exp(1i*thext);
bdry_aaa = mn_aaa*exp(1i*thext);

plot(real(bdry),imag(bdry),'k','LineWidth',4)
%plot(real(bdry_aaa),imag(bdry_aaa),'Color',colors(1,:),'LineWidth',2,'LineStyle','--')
%plot(real(bdry_qr),imag(bdry_qr),'Color',colors(2,:),'LineWidth',2,'LineStyle','--')
%legend('aaa','qr')
axis('equal')
axis('off')
%%

rr= 1.75;
[~,indr] = min(abs(rext-rr));

rr_ext = rext(indr);
fext_qr = r(rext);
fext_qr = fext_qr(indr,:).';
fext_qr = PHI*fext_qr;
f_exact = f_polar(rr_ext,thext.');

figure(15);clf
plot(thext,f_exact)
hold on
plot(thext,fext_qr)
legend('exact','qr')



%%

function g=gamm(r,t0)

t= t0(:);
%t = (1+(t/pi))*4;
t2 = t(t<2);
t3 = t(t>=2 & t<3);
t4 = t(t>=3 & t<4);
t5 = t(t>=4 & t<5);
t6 = t(t>=5 & t<6);
t8 = t(t>=6 & t<=8);


g2 = [t2,zeros(numel(t2),1)];
g3 = [2*ones(numel(t3),1),t3-2];
g4 = [-(t4-3)+2,ones(numel(t4),1)];
g5 = [ones(numel(t5),1),t5-4+1];
g6 = [-(t6-5)+1,2*ones(numel(t6),1)];
g8 = [zeros(numel(t8),1),-(t8-7)+1];
g=[g2;g3;g4;g5;g6;g8];
avg = [.5,.5];%sum(g,1)/size(g,1);
g = g-avg;
g = g(:,1)+(1i)*g(:,2);
if size(t0,2)>1
    g=g.';
end
g = r.*g;
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


function PHI = build_PHI(degr,th)
PHI = zeros(numel(th),degr);
th0 = (th-min(th))/(max(th)-min(th));
th0 = 2*th0-1;
for i=1:degr
    cc = zeros(degr,1);
    cc(i) = 1;
    c = chebfun(cc,'coeffs');
    PHI(:,i) = c(th0);
end
end