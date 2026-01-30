function [rxy,info] = construct_multi_AAA_ID(pts,F,f,options)

% Function that computes multivariate rational approximation
%
%
%
%
%       INPUTS      
%                   pts :       1xd cell array
%                   F   :       N1x...xNd tensor of function values
%                   f   :       function that can be evaluated at arbitrary
%                               points (currently needed for twostep approximation)
%                   options:    options for type of rational approx
%                               options.tolID   :   tol for the initial ID decompositions, can be cell array or single number for all dimensions
%                               options.tolAAA  :   tol for component-wise AAA approximations, can be cell array or single number for all dimensions 
%                               options.mmax    :   max degree, can be cell array or single number for all directions
%                               options.twostep :   boolean, whether to use
%                                                   twostep approx
%                               options.paaapts :   1xd cell array of grid for construction of twostep approx
%                               options.validation: whether to pass validation points to paaa, used for convergence stats
%                               options.valpts  :   grid of validation pts
%                               options.Fvalpts :   corresponding functions values


dim = numel(size(F));
[twostep,tolID,tolAAA,mmax] = parse_opts(options,dim);

dim = size(pts,2);

% init cells
wcell = cell(1,dim);
Icell = cell(1,dim);
nodes = cell(1,dim);
Supp = cell(1,dim);
info.time_aaa = 0;
info.rk = cell(1,dim);
tic
INDS = ID(F,tolID);
info.time_ID = toc;


toc
INDSTOT = cell(1,dim);
for d=1:dim
    INDSTOT{d} = 1:size(F,d);
    info.rk{d} = numel(INDS{d});
end
for d=1:dim
    tic
    INDSsub = INDS;
    INDSsub{d} = INDSTOT{d};
    Fsub = F(INDSsub{:});
    Fd = tenmat(Fsub,d).data;
    %Fd = Fd./vecnorm(Fd,inf,1);%rescale if functions in F are badly scaled
    tic
    [~, ~, ~, ~, zd, ~, wd, ~] = aaa_sv(Fd,pts{d},'tol',tolAAA{d}*sqrt(size(Fd,2)),'mmax',mmax{d});
    info.time_aaa = info.time_aaa+toc;
    Supp{d} = zd;
    Id = [];
    for i=1:numel(zd)
        Id = [Id,find(pts{d}==zd(i))];
    end
    Icell{d} = Id;
    wcell{d} = wd;
    nodes{d} = pts{d}(Id);
end
info.Supp = Supp;

ff = F(Icell{:});
w = (full(ktensor(1,wcell{:})));
rxy = BarycentricForm(nodes,w.*ff,w);
if twostep
    paaapts = cell(1,dim);
    itpl_part = cell(1,dim);

    if isfield(options,'paaapts')
        ppts = options.paaapts;
    else
        ppts = pts;
    end
    for d=1:dim
        paaapts{d} = unique([nodes{d},ppts{d}],'stable');
        itpl_part{d} = 1:numel(nodes{d});
    end
    % toggle depending on if function takes cell input or ndgrid input
    PAAAPTS = cell(dim,1);
    [PAAAPTS{:}] = ndgrid(paaapts{:});
    Fpaaa = f(PAAAPTS{:});
    %Fpaaa = f(paaapts{:});
    optionspaaa.itpl_part = itpl_part;
    optionspaaa.max_iter = 50;
    if validation
        optionspaaa.validation.sampling_values = valpts;
        optionspaaa.validation.samples = Fvalpts.data;
    end
    [rxy,infopaaa] = paaa(Fpaaa,paaapts,tolpaaa,optionspaaa);
    info.Supp = rxy.itpl_nodes;
    info.infopaaa = infopaaa;
end



end


function I = ID(T,tol)
dim = numel(size(T));
    I = cell(1,dim);
    Isub = cell(1,dim);
    for d=1:dim
        Isub{d} = 1:size(T,d);
    end
    for d=1:dim
        Tsub = T(Isub{:});
        [~,R,piv] = qr((tenmat(Tsub,d).data).','econ','vector');
        dr = abs(diag(R));
        k = sum(dr>tol{d}*dr(1));
        I{d} = piv(1:k);
        Isub{d} = piv(1:k);
    end
    
end



function [twostep,tol_qr,tolAAA,mmax] = parse_opts(options,dim)

if ~isfield(options,'twostep')
    twostep=false;
else
    twostep=options.twostep;
end

if ~isfield(options,'tol_qr')
    tol_qr = cell(1,dim);
    tol_qr(:)={1e-15};
else
    if ~iscell(options.tol_qr)
        tol_qr = cell(1,dim);
        tol_qr(:) = {options.tol_qr};
    else
        tol_qr = options.tol_qr;
    end
end

if ~isfield(options,'tolAAA')
    tolAAA = cell(1,dim);
    tolAAA(:)={1e-13};
else
    if ~iscell(options.tolAAA)
        tolAAA = cell(1,dim);
        tolAAA(:) = {options.tolAAA};
    else
        tolAAA = options.tolAAA;
    end
end

if ~isfield(options,'mmax')
    mmax = cell(1,dim);
    mmax(:)={100};
else
    if ~iscell(options.mmax)
        mmax = cell(1,dim);
        mmax(:) = {options.mmax};
    else
        mmax = options.mmax;
    end
end



end