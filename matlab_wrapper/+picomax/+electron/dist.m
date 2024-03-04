function [] = dist(o,qvec,lambda,varargin)

%-- input parser
p = inputParser;
addParameter(p,'bloch',true);
addParameter(p,'resolution',51); % resolution of spatial grids
addParameter(p,'grid',0);
parse(p,varargin{:});
p = p.Results;


%-- real-space grids
N = p.resolution; % resolution of spatial grids
if p.grid==0 % primitive cell
    df = 1/N;
    [f1,f2,f3] = ndgrid(0:df:1-df,0:df:1-df,0:df:1-df);
    a1 = o.sys.lattice_vector(:,1).';
    a2 = o.sys.lattice_vector(:,2).';
    a3 = o.sys.lattice_vector(:,3).';
    r = a1.*f1(:)+a2.*f2(:)+a3.*f3(:);
    Nr = size(r,1);
elseif p.grid==1
    dx = 1/N;
    x = 0:dx:1-dx;
    y = x;
    z = x;
    [x,y,z] = ndgrid(x,y,z);
    r = [x(:),y(:),z(:)];
    Nr = size(r,1);
end


%-- plane-wave grids
gvec = o.dat.electron.gvec;
NPW = size(gvec,2);

%-- q vector
if ~isnumeric(qvec)
    qvec = picomax.highsym.fcc(qvec);
end
dQ = vecnorm(o.num.QVEC-qvec(:),1,1);
[~,q] = min(dQ);
qvec = o.num.QVEC(:,q);

% eigenvector component
Cql = o.dat.electron.V(:,lambda,q);

psi = zeros(Nr,1);
for i=1:NPW
    qg = qvec+gvec(:,i);
    psi = psi+Cql(i)*exp(1i*2*pi*r*qg);
end

%-- save
o.dat.electron.psi = reshape(psi,[N N N]);
o.dat.electron.x = reshape(r(:,1),[N N N]);
o.dat.electron.y = reshape(r(:,2),[N N N]);
o.dat.electron.z = reshape(r(:,3),[N N N]);
o.dat.electron.Q = qvec;

