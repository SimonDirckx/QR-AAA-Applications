%% extension for the norm function
% uses qr-aaa and trig interp

clear
run init.m
nx = 200;
ny = 401;

x=linspace(-4,-2,nx);
y=linspace(-4,4,ny);
f = @(x,y)sqrt(x.^2+y.^2);

[X,Y] = ndgrid(x,y);
F = tensor(f(X,Y));




%% Apply multivariate QR-AAA
clear options
options.interp=true;
options.tolID = 1e-14;
options.tolAAA = {1e-11,1e-11};
[rxy,infoxy] = construct_multi_AAA_ID({x,y},F,f,options);
R = rxy.eval({x,y});
norm(F(:)-R(:))/norm(F(:))


%% Apply 2step QR-AAA
clear options
options.interp=true;
options.twostep = true;
options.tol_qr = 1e-15;
options.mmax = {4,12}; %(half of support in qr)
options.tolpaaa = 1e-11;
options.sep=true;
[rxy_paaa,info_paaa] = construct_multi_AAA({x,y},F,f,options);
R = rxy_paaa.eval({x,y});
norm(F(:)-R(:))/norm(F(:))



%%



nr_ext = 2*nx;
xext = linspace(-4,4,nr_ext);
yext = linspace(-4,4,2*nr_ext);

Fapprox = rxy.eval({xext,yext});

[Xext,Yext] = ndgrid(xext,yext);
Fext = f(Xext,Yext);

norm(Fext-Fapprox,inf)


figure(1);clf;
contour(Xext,Yext,log10(abs(Fapprox-Fext)),30,'LineWidth',2)
hold on
plot([-2 -2], [-4,4],'k',Linewidth=4)
colorbar
%% compare at xpt
xpt = .1;

fx = zeros(numel(yext),1);

for i = 1:size(fx,1)

[r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(f(x,yext(i)),x,'tol',1e-13);
fx(i) = r_aaa(xpt);

end

f_exact = f(xpt,yext);
f_xy = rxy.eval({xpt,yext});
f_xy_paaa = rxy_paaa.eval({xpt,yext});

figure(2);clf;
semilogy(yext,abs(f_exact-fx.'))
hold  on
semilogy(yext,abs(f_exact-f_xy))
semilogy(yext,abs(f_exact-f_xy_paaa))
legend('aaa','qr-aaa','2step')


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
