function [r,w,smin] = compute_rat(Zx,Zy,Ix,Iy,F,f)
    Zsuppx = Zx(Ix);
    Zsuppy = Zy(Iy);
    Zrestx = Zx;
    Zresty = Zy;
    Zrestx(Ix) = [];
    Zresty(Iy) = [];
    Fsupp = F(Ix,Iy).';
    Frest = F;
    Frest(Ix,:)=[];
    Frest(:,Iy)=[];
    Frest=Frest.';
    frest = Frest(:);
    fsupp=Fsupp(:);
    Cx = 1./(Zrestx.'-Zsuppx);
    Cy = 1./(Zresty.'-Zsuppy);
    Cxy = kron(Cx,Cy);
    Cx1 = kron(Cx,eye(numel(Zsuppy)));
    C1y = kron(eye(numel(Zsuppx)),Cy);
    
    fx = f(Zrestx,Zsuppy.');
    fy = f(Zsuppx,Zresty.');
    Lx1 = fx(:).*Cx1 - Cx1.*fsupp.';
    L1y = fy(:).*C1y - C1y.*fsupp.';
    L = (frest.*Cxy-Cxy.*fsupp.');
    L=[L;Lx1;L1y];
    disp('rank L')
    rank(L)
    size(L)
    [~,S,V] = svd(L,0);
    s = diag(S);
    smin = s(end);
    ii = find(s<=max(1e-15*s(1),s(end)));
    disp('nii')
    numel(ii)
    smin = smin*s(1);
    w = V(:,ii);
    w = w*ones(size(w,2),1)/sqrt(size(w,2));
    %wmat = zeros(numel(Ix),numel(Iy));
    %wmat(:) = w;
    %[U,S,V] = svd(wmat);
    %s = diag(S);
    %k = sum(s>1e-12*s(1));
    %wmat = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    %w = wmat(:);
    r=@(x,y) rat_eval(x,y,w,Zx,Zy,Ix,Iy,F);
end

function r=rat_eval(x,y,w,Zx,Zy,Ix,Iy,F)
    x=x(:);
    y=y(:);
    Zx=Zx(:);
    Zy=Zy(:);
    Zsuppx = Zx(Ix);
    Zsuppy = Zy(Iy);
    
    r = zeros(numel(x),numel(y));
    iix = [];
    
    for i = 1:numel(Zsuppx)
        iix = [iix,find(x==Zsuppx(i))];
    end
    iiy = [];
    for i = 1:numel(Zsuppy)
        iiy = [iiy,find(y==Zsuppy(i))];
    end
    jjx = [];
    for i=1:numel(iix)
        jjx=[jjx,find(x(iix(i))==Zsuppx)];
    end
    jjy = [];
    for i=1:numel(iiy)
        jjy=[jjy,find(y(iiy(i))==Zsuppy)];
    end
    
    xrest = x;
    yrest = y;
    xrest(iix) = [];
    yrest(iiy) = [];
    IIx = 1:numel(x);
    IIy = 1:numel(y);
    IIx(iix) = [];
    IIy(iiy) = [];
    Fsupp = F(Ix,Iy);
    ff = Fsupp.';
    ff = ff(:);
    Cx = 1./( xrest- Zsuppx.');
    Cy = 1./( yrest- Zsuppy.');
    Cxy = kron(Cx,Cy);
    Cx1 = kron(Cx,eye(numel(Zsuppy)));
    C1y = kron(eye(numel(Zsuppx)),Cy);


    N=Cxy*(w.*ff);
    D=Cxy*w;
    rvec = N./D;


    Nx1=Cx1*(w.*ff);
    Dx1=Cx1*w;
    rvecx1 = Nx1./Dx1;

    N1y=C1y*(w.*ff);
    D1y=C1y*w;
    rvec1y = N1y./D1y;


    r(IIx,IIy) = reshape(rvec,[numel(IIy),numel(IIx)]).';
    r(IIx,iiy) = reshape(rvecx1,[numel(iiy),numel(IIx)]).';
    r(iix,IIy) = reshape(rvec1y,[numel(IIy),numel(iix)]).';
    r(iix,iiy) = f(Zsuppx,Zsuppy.');
end