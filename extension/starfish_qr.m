%% new starfish script
% uses qr-aaa and trig interp

clear
run init.m
nr = 500;
nt = 500;
rmin = 1e-1;%.5; %starting point for strip along the boundary
rvec = linspace(rmin,1,nr).';
nr = numel(rvec);
th = linspace(-pi,pi,nt).';
nt = numel(th);


gamm = @(r,th) r.*(1+.3*sin(5*th)).*exp(1i*th);
f = @(x,y)-sin(2*pi*x).*sin(2*pi*y);
f_polar = @(r,th) f(real(gamm(r,th)),imag(gamm(r,th)));


c1=(1+(.3)*sin(5*th)).*exp(1i*th);

rho = zeros(numel(rvec)*numel(th),2);
rho(:,1) = kron(rvec,ones(size(th)));
rho(:,2) = kron(ones(size(rvec)),th);
XY = rho(:,1).*[cos(rho(:,2)),sin(rho(:,2))];


gXY = zeros(size(rho));
gXY(:,1) = real(gamm(rho(:,1),rho(:,2)));
gXY(:,2) = imag(gamm(rho(:,1),rho(:,2)));
%F = f(gXY(:,1),gXY(:,2));
tic
Fmat = f_polar(rvec,th.');
toc
F=Fmat.';
F=F(:);
%Fmat=reshape(F,[nt,nr]).';
figure(1)
clf
scatter(XY(:,1),XY(:,2))
hold on
plot(real(c1),imag(c1),'r','linewidth',4)
axis('equal')

figure(2)
clf
scatter(gXY(:,1),gXY(:,2))
hold on
plot(real(c1),imag(c1),'r','linewidth',4)
axis('equal')

figure(3)
clf
scatter(gXY(:,1),gXY(:,2),[],F)
axis('equal')

figure(4)
clf
scatter(XY(:,1),XY(:,2),[],F)
axis('equal')

tic
om_max = min(200,(nt-1)/2);
PHI = exp(1i*th.*(-om_max:om_max))/sqrt(nt);
Ftild = (PHI\Fmat.').';
toc

%% Apply QR-AAA with trig interp


use_qr=true;

if use_qr
    tol_qr = 1e-13;
    tic
    [Q,R,P] = qr(Ftild,'econ','vector');
    
    k_qr = sum(abs(diag(R))>tol_qr*abs(R(1,1)));
    Fk = Q(:,1:k_qr)*diag(diag(R(1:k_qr,1:k_qr)));
    toc
    
else
    tol = 1e-14;
    dg = zeros(size(Ftild,2),1);
    for i=1:size(Ftild,2)
    dg(i) = norm(Ftild(:,i),inf);
    end
    ii = find(abs(dg)<tol);
    dg(ii)=[];
    Ftild0=Ftild;
    Ftild0(:,ii)=[];
    Fk = Ftild0;
    R=eye(size(Ftild,2));
end
tic
[r, pol, res, zer, z, ff, w, errvec] = aaa_sv(Fk,rvec,'tol',1e-13);
toc
ind_supp = [];
for i=1:numel(z)
    ind_supp=[ind_supp,find(rvec==z(i))];
end
r = @(zz) reval_sv(zz, z, Ftild(ind_supp,:), w);


%%
norm(r(rvec)-Ftild)/norm(Ftild)


nr_ext = 2*nr;
rext = linspace(rmin,2.5,nr_ext).';
thext = linspace(0,2*pi,3000).';

PHIext = zeros(numel(thext),2*om_max+1);
for i=-om_max:om_max
    PHIext(:,om_max+i+1) = exp(1i*i*thext)/sqrt(nt);
end

Rext = r(rext);
Rext = PHIext*(Rext.');
Fapprox = Rext(:);

rho = zeros(numel(rext)*numel(thext),2);
rho(:,1) = kron(rext,ones(size(thext)));
rho(:,2) = kron(ones(size(rext)),thext);
XYext = rho(:,1).*[cos(rho(:,2)),sin(rho(:,2))];


gXYext = zeros(size(rho));
gXYext(:,1) = real(gamm(rho(:,1),rho(:,2)));
gXYext(:,2) = imag(gamm(rho(:,1),rho(:,2)));
Fext = f(gXYext(:,1),gXYext(:,2));

norm(Fext-Fapprox,inf)


figure(5)
clf
scatter(gXYext(:,1),gXYext(:,2),[],Fext)
hold on
plot(real(c1),imag(c1),'k','linewidth',4)
axis('equal')
axis('off')

figure(6)
clf
scatter(XYext(:,1),XYext(:,2),[],Fext)
axis('equal')

figure(7)
clf
scatter(gXYext(:,1),gXYext(:,2),[],Fapprox)
axis('equal')

figure(8)
clf
scatter(gXYext(:,1),gXYext(:,2),[],Fapprox-Fext)
hold on
plot(real(c1),imag(c1),'r','linewidth',4)
axis('equal')
colorbar

figure(9)
clf
scatter(gXYext(:,1),gXYext(:,2),[],log10(abs(Fapprox-Fext)))
hold on
plot(real(c1),imag(c1),'r','linewidth',4)
axis('equal')
colorbar

figure(10)
clf
scatter(real(pol),imag(pol))
axis('equal')

figure(11)
clf
plot(real(c1),imag(c1),'r','linewidth',3)
hold on
for i=1:numel(z)
    g = gamm(z(i),thext);
    plot(real(g),imag(g),'k--','linewidth',1)
end
axis('equal')
%% compare along rays
addpath ../chebfun-master
th0 = thext(10); %select th not in approx domain
rmin_aaa = 0;
rext = linspace(rmin,4,1000).';
rvec_aaa = linspace(rmin_aaa,1,nr).';
gamm_vec_aaa = gamm(rvec_aaa,th0);
fvec_aaa = f(real(gamm_vec_aaa),imag(gamm_vec_aaa));
[r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(fvec_aaa,rvec_aaa,'tol',1e-13);



PHIvec = zeros(1,2*om_max+1);
for i=-om_max:om_max
    PHIvec(:,om_max+i+1) = exp(1i*i*th0)/sqrt(nt);
end

fext_qr = r(rext);
fext_qr = PHIvec*(fext_qr.');


gamm_ext = gamm(rext,th0);
fext_aaa = r_aaa(rext);
figure(12);clf;
plot(abs(gamm_ext),fext_qr)
hold on
plot(abs(gamm_ext),fext_aaa)
plot(abs(gamm_ext),f(real(gamm_ext),imag(gamm_ext)))
legend('qr','aaa','ref')


fext_qr = r(rext);
fext_qr = (PHIvec*(fext_qr.')).';


gamm_ext = gamm(rext,th0);
fext_aaa = r_aaa(rext);
f_exact = f(real(gamm_ext),imag(gamm_ext));


loss_ind_aaa=0;
loss_ind_qr=0;
err = 0;
while err<1
    loss_ind_aaa=loss_ind_aaa+1;
    err = abs(f_exact(loss_ind_aaa)-fext_aaa(loss_ind_aaa));
end
err=0;
while err<1 && loss_ind_qr<numel(fext_qr)
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

%% compare loss contours at specified precision

loss_prec =1e-1;
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
    PHIvec = zeros(1,2*om_max+1);
    for i=-om_max:om_max
        PHIvec(:,om_max+i+1) = exp(1i*i*thext(ind))/sqrt(nt);
    end
    
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
