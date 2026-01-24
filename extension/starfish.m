%% ray-AAA approximation for the starfish
clear
run init.m
nr = 100;
nt = 501;
nt_coarse = 80;
rmin = 0;
rmax_ext = 2.5;
rvec = linspace(rmin,1,nr).';
nr = numel(rvec);
th = linspace(-pi,pi,nt).';
th_coarse = linspace(-pi,pi,nt_coarse).';
gamm = @(r,th) r.*(1+.3*sin(5*th)).*exp(1i*th);
f = @(x,y)-sin(2*pi*x).*sin(2*pi*y);
c1=gamm(1,th);

rho = zeros(numel(rvec)*numel(th),2);
rho(:,1) = kron(rvec,ones(size(th)));
rho(:,2) = kron(ones(size(rvec)),th);
XY = rho(:,1).*[cos(rho(:,2)),sin(rho(:,2))];

gXY = zeros(size(rho));
gXY(:,1) = real(gamm(rho(:,1),rho(:,2)));
gXY(:,2) = imag(gamm(rho(:,1),rho(:,2)));

F = f(gXY(:,1),gXY(:,2));
Fmat=reshape(F,[nt,nr]).';
z = rmax_ext*exp(1i*th_coarse);

figure(1)
clf
plot(real(c1),imag(c1),'k','linewidth',4)
hold on
for i =1:numel(z)
plot([0,real(z(i))],[0,imag(z(i))],'k','LineWidth',2)
end
axis('equal')
axis('off')

figure(2)
clf
hold on
plot(real(c1),imag(c1),'k','linewidth',4)
hold on
for i =1:numel(z)
plot([0,real(z(i))],[0,imag(z(i))],'k','LineWidth',2)
end
scatter(gXY(:,1),gXY(:,2),[],F,'filled')
axis('equal')
axis('off')

%%
Next = 500;

Fext = zeros(numel(th),Next);
Fexact = zeros(numel(th),Next);

rext = linspace(0,rmax_ext,Next).';
rhoext = zeros(numel(rext)*numel(th),2);
rhoext(:,1) = kron(rext,ones(size(th)));
rhoext(:,2) = kron(ones(size(rext)),th);
XYext = rhoext(:,1).*[cos(rhoext(:,2)),sin(rhoext(:,2))];


for ind=1:numel(th)
    th0 = th(ind);
    rmax = (1+.3*sin(5*th0));
    rvec = linspace(0,rmax,nr);
    
    fvec = f(rvec*cos(th0),rvec*sin(th0));
    [r_aaa, pol_aaa, res_aaa, zer_aaa, z_aaa, ff_aaa, w_aaa, errvec_aaa] = aaa(fvec,rvec);
    
    rr = r_aaa(rext);
    rr(rr>1.)=1.;
    rr(rr<-1.)=-1.;
    Fext(ind,:) = rr;
    Fexact(ind,:) = f(rext*cos(th0),rext*sin(th0));

end

%%
figure(503)
clf
scatter(XYext(:,1),XYext(:,2),[],Fext(:) ,'filled')
hold on
plot(real(c1),imag(c1),'k','linewidth',4)
axis('equal')
axis('off')

%%
figure(504)
clf
scatter(gXYext(:,1),gXYext(:,2),[],Fext(:),'filled' )
hold on
plot(real(c1),imag(c1),'k','linewidth',4)
axis('equal')
axis('off')

figure(505)
clf
plot(th,Fext(:,900))