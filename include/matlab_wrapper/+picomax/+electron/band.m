% compute electron band structure

function [] = band(o,varargin)

% input parser
p = inputParser;
addParameter(p,'sort',true);
addParameter(p,'offset',true);
parse(p,varargin{:})
p = p.Results;




% generate plane-wave-grids
KCON = 3.809982111485962;
GCUT = sqrt(o.num.encut/KCON)/(2*pi/o.sys.a);
gvec = picomax.kpoint.pwgrid_fcc(GCUT);
gvec = gvec.'; % g [3 x NPW]
NPW = size(gvec,2);

% generate q-vectors
picomax.qvector.generate(o);
qvec = o.num.qvec; % q [3 x NQ]
NQ = size(qvec,2);
[~,qg] = min(sum(abs(qvec).^2,1)); % find gamma point index

% construct hamiltonian and solve EVP
H = zeros([NPW NPW NQ]);
for q=1:NQ
    H(:,:,q) = picomax.electron.HamiltonianEPM(gvec,qvec(:,q),o.sys.epm,o.sys.a);
end
% V1 [NPW x NPW x NQ] eigenvectors
% D1 [NPW x NQ] eigenvalues

% sort bands and truncate to NBAND
NBAND = min([NPW o.num.nband]);
if p.sort
    [V1,D1] = picomax.util.eigenshuffle(H);
    [~,idx] = sort(D1(:,qg));
    idx = idx(1:NBAND);
    V = V1(:,idx,:);
    D = D1(idx,:);
else
    V = zeros([NPW NBAND NQ]);
    D = zeros([NBAND NQ]);
    for q=1:NQ
        [V1,D1] = eig(H(:,:,q));
        V(:,:,q) = V1(:,1:NBAND);
        D1 = diag(D1);
        D(:,q) = D1(1:NBAND);
    end
end


% energy offset
if p.offset
    Eoffset = D(4,qg);
    D = D-Eoffset;
    o.dat.electron.Eoffset = Eoffset;
end

% save
o.dat.electron.V = V;
o.dat.electron.D = D;
o.dat.electron.gvec = gvec;
