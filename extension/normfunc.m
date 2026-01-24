%% qr-aaa extension for sqrt(x^2+y^2)

f = @(x,y)sqrt(x.^2+y.^2);

xpts = linspace(-4,-2,500).';



ypts = 2*chebpts(100);
ypts = [(ypts-2);(2+ypts)];
ypts = sort(unique(ypts));
XY = zeros(numel(xpts)*numel(ypts),2); 

XY(:,1) = kron(xpts,ones(size(ypts)));
XY(:,2) = kron(ones(size(xpts)),ypts);

Fmat = f(xpts,ypts.');

[Q,R,P]=qr(Fmat,"econ","vector");
Gamm = diag(R);
k = sum(abs(Gamm)>1e-15*Gamm(1));
Qk = Q(:,1:k)*diag(Gamm(1:k));
Rk = diag(Gamm(1:k))\R(1:k,:);
Rk(:,P)=Rk;


%Qk = Fmat*randn(size(Fmat,2),5);

[r, pol, res, zer, z, ff, w, errvec] = aaa_sv(Qk,xpts,'tol',1e-11);

ind_qr = [];
for i=1:numel(z)
    ind_qr = [ind_qr,find(xpts==z(i))];
end

r=@(zz) reval_sv(zz, z, Fmat(ind_qr,:), w);
xext = linspace(-4,4,1000).';
yext = ypts;%linspace(-1,1,500).';
Fapprox = r(xext).';



[X,Y] = meshgrid(xext,yext);
%figure(1)
%contour(X,Y,Fapprox.',.25:.25:3,'blue')
figure(2);clf;
contour(X,Y,log10(abs(Fapprox-f(xext,yext.').')),30)
colorbar
F_aaa=zeros(numel(yext),numel(xext));
for i=1:numel(yext)
    fvec = f(xpts,yext(i));
    [r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(fvec,xpts,'tol',1e-13);
    F_aaa(i,:) = r_aaa(xext);
end
figure(3);clf;

contour(X,Y,log10(abs(F_aaa-f(xext,yext.').')),30)
colorbar
%%
[~,yind] = min(abs(ypts-1.)); 
y0 = ypts(yind);
figure(4);clf;
plot(xext,F_aaa(yind,:))
hold on
plot(xext,Fapprox(yind,:))
plot(xext,f(xext,y0))
legend('aaa','qr','f')

figure(5);clf;
semilogy(xext,abs(F_aaa(yind,:).'-f(xext,y0))./abs(f(xext,y0)))
hold on
semilogy(xext,abs(Fapprox(yind,:).'-f(xext,y0))./abs(f(xext,y0)))
legend('aaa','qr')





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




