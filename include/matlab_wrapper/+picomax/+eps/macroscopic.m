








function [] = macroscopic(o)

neps = o.num.neps;
nfreq = o.num.nfreq;
nq = sum(o.num.qnum)+1;

epsmnL = o.dat.realepsL+1i*o.dat.imagepsL;
epsmnT = o.dat.realepsT+1i*o.dat.imagepsT;
if nq==1 && norm(o.num.qvec)==0
    epsmnL(:,:,1) = epsmnT(:,:,1);
end

%-- longitudinal
invepsmnL = zeros(nfreq,nq,neps,neps);
for f=1:nfreq
    for q=1:nq
        epsmn_L = zeros(neps);
        for m=1:neps
            for n=1:m
                idx = (m-1)*m/2+n;
                epsmn_L(m,n) = epsmnL(f,q,idx);
                if (m~=n)
                    epsmn_L(n,m) = conj(epsmn_L(m,n));
                end
            end
        end
        invepsmnL(f,q,:,:) = inv(epsmn_L);
    end
end

%-- transverse
g = [0     0     0
    -1    -1    -1
    -1    -1     1
    -1     1    -1
    -1     1     1
     1    -1    -1
     1    -1     1
     1     1    -1
     1     1     1];
a = o.sys.a*1e-10; % lattice constant [m]
c0 = 3e8; % speed of light [m/s]
hbar_eV = 4.135667696e-15; % planck constant [eV*s]
w = o.dat.omega(:)/hbar_eV; % frequency [rad/s]
qvec = o.num.qvec;
invepsmnT = zeros(nfreq,nq,neps,neps);
for f=1:nfreq
    for q=1:nq
        epsmn_T = zeros(neps);
        for m=1:neps
            for n=1:m
                idx = (m-1)*m/2+n;
                epsmn_T(m,n) = epsmnT(f,q,idx);
                if (m~=n)
                    epsmn_T(n,m) = conj(epsmn_T(m,n));
                end
            end
        end
        % diagonal terms
        for m=1:neps
            qg = qvec(:,q)+g(m,:).';
            qg = qg*2*pi/a;
            norm2qg = sum(abs(qg).^2,1);
            if norm2qg==0
                norm2qg = 1e-15;
            end
            epsmn_T(m,m) = epsmn_T(m,m) ...
                -c0^2/w(f)^2*norm2qg;
        end

        invepsmnT(f,q,:,:) = inv(epsmn_T);
    end
end



%-- save
o.dat.eps.eps00 = epsmnL(:,:,1,1);
o.dat.eps.epsM = 1./invepsmnL(:,:,1,1);
o.dat.eps.epsM_T = 1./invepsmnT(:,:,1,1);
end