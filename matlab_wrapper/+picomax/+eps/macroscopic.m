








function [] = macroscopic(o)

NEPS = o.num.NEPS;
NFREQ = o.num.NFREQ;
NQ = sum(o.num.QNUM)+1;

epsmnL = o.dat.realepsL+1i*o.dat.imagepsL;
epsmnT = o.dat.realepsT+1i*o.dat.imagepsT;
if NQ==1 && norm(o.num.QVEC)==0
    epsmnL(:,:,1) = epsmnT(:,:,1);
end

%-- longitudinal
invepsmnL = zeros(NFREQ,NQ,NEPS,NEPS);
for f=1:NFREQ
    for q=1:NQ
        epsmn_L = zeros(NEPS);
        for m=1:NEPS
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
qvec = o.num.QVEC;
invepsmnT = zeros(NFREQ,NQ,NEPS,NEPS);
for f=1:NFREQ
    for q=1:NQ
        epsmn_T = zeros(NEPS);
        for m=1:NEPS
            for n=1:m
                idx = (m-1)*m/2+n;
                epsmn_T(m,n) = epsmnT(f,q,idx);
                if (m~=n)
                    epsmn_T(n,m) = conj(epsmn_T(m,n));
                end
            end
        end
        % diagonal terms
        for m=1:NEPS
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