function [] = compute(o)


% generate plane-wave-grids
gvec = picomax.kpoint.pwgrid_fcc(o.num.GCUT);
gvec = gvec.'; % g [3 x NPW]

% generate q-vectors
picomax.qvector.generate(o);
qvec = o.num.QVEC; % q [3 x NQ]
NQ = size(qvec,2);

NPW = size(gvec,2);
V = zeros(NPW,NPW,NQ);
D = zeros(NPW,NQ);

for q=1:NQ
    H = picomax.electron.HamiltonianEPM(gvec,qvec(:,q),o.sys.epm,o.sys.a);
    [V(:,:,q),D1] = eig(H);
    D(:,q) = diag(D1);
end

o.dat.electron.V = V;
o.dat.electron.D = D;
o.dat.electron.gvec = gvec;







