








function [] = polarizationdist(o,Q,omega,lambda,varargin)

%-- input parser
p = inputParser;
addParameter(p,'bloch',true);
addParameter(p,'resolution',51); % resolution of spatial grids
addParameter(p,'longitudinal',true);
addParameter(p,'grid',0);
parse(p,varargin{:});
p = p.Results;


% preallocation
NEPS = o.num.NEPS;
N = p.resolution; % resolution of spatial grids


%-- grids
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


g = [0     0     0
    -1    -1    -1
    -1    -1     1
    -1     1    -1
    -1     1     1
     1    -1    -1
     1    -1     1
     1     1    -1
     1     1     1];







dQ = vecnorm(o.num.QVEC-Q(:),1,1);
[~,q] = min(dQ);
Q = o.num.QVEC(:,q);
o.dat.Q = Q;

f = find(o.dat.omega>=omega,1);



% % V0 (f = 0)
% chimn = zeros(NEPS);
% for m=1:NEPS
%     for n=1:m
%         idx = (m-1)*m/2+n;
%         chimn(m,n) = o.dat.realepsL(1,q,idx)+1i*o.dat.imagepsL(1,q,idx);
%         if (m~=n)
%             chimn(n,m) = conj(chimn(m,n));
%         end
%         if (m==n)
%             chimn(n,m) = chimn(n,m)-1;
%         end
%     end
% end
% [V0,~] = eig(chimn);

if p.longitudinal
    V = o.tmp.VL(:,:,f,q);
    alpha = o.tmp.DL(:,f,q)-1;
else
    V = o.tmp.VT(:,:,f,q);
    alpha = o.tmp.DT(:,f,q)-1;
end





% chimn = zeros(NEPS);
% for m=1:NEPS
%     for n=1:m
%         idx = (m-1)*m/2+n;
%         chimn(m,n) = o.dat.realepsL(f,q,idx)+1i*o.dat.imagepsL(f,q,idx);
%         if (m~=n)
%             chimn(n,m) = conj(chimn(m,n));
%         end
%         if (m==n)
%             chimn(n,m) = chimn(n,m)-1;
%         end
%     end
% end
% [V,~] = eig(chimn);
% alphaL = diag(V0'*chimn*V0);
% [V,alphaL] = eig(chimn);
% alphaL = diag(alphaL);
% alphaL = diag(V0'*chimn*V0);

%
px = zeros(Nr,1);
py = zeros(Nr,1);
pz = zeros(Nr,1);
for m=1:NEPS
    n = Q(:).'+g(m,:);
    if p.longitudinal % longitudinal polarization
        if norm(n)==0
            n = [1 0 0];
        end
        normp = n/norm(n);
    else % transverse polarization
        if norm(n)==0
            n = [0 1 0];
        end
        normp = picomax.util.tangentialvector(n);
    end
    if ~p.bloch
        n = n-Q(:).';
    end
    tmp = normp*V(m,lambda).*exp(1i*2*pi*n.*r);
    px = px+tmp(:,1);
    py = py+tmp(:,2);
    pz = pz+tmp(:,3);
end
px = (alpha(lambda))/(4*pi)*px;
py = (alpha(lambda))/(4*pi)*py;
pz = (alpha(lambda))/(4*pi)*pz;


o.dat.px = reshape(px,[N N N]);
o.dat.py = reshape(py,[N N N]);
o.dat.pz = reshape(pz,[N N N]);
o.dat.x = reshape(r(:,1),[N N N]);
o.dat.y = reshape(r(:,2),[N N N]);
o.dat.z = reshape(r(:,3),[N N N]);





end

