








function [] = band(o)
NEPS = o.num.NEPS;
NFREQ = o.num.NFREQ;
NQ = sum(o.num.QNUM)+1;
% o.initialize

o.sys.initialize
picomax.qvector.generate(o)

%-- construct [[eps]] matrix
f = 1;
epsmn = zeros([NEPS NEPS NQ]); % f = 0
for q=1:NQ
    % fill lower triangle
    for m=1:NEPS
        for n=1:m
            idx = (m-1)*m/2+n;
            epsmn(m,n,q) = o.dat.realepsL(f,q,idx)+1i*o.dat.imagepsL(f,q,idx);
        end
    end

    % enforce symmetry
    % epsmn(:,:,q) = picomax.eps.enforcesym(epsmn(:,:,q),o.num.QVEC(:,q));

    % fill uppder triangle
    for m=1:NEPS
        for n=1:m-1
            epsmn(n,m,q) = conj(epsmn(m,n,q));
        end
    end
end
[V0,D0] = picomax.util.eigenshuffle(epsmn);


epsmn = zeros([NEPS NEPS NFREQ NQ]);
for q=1:NQ
    for f=1:NFREQ
        % fill lower triangle
        for m=1:NEPS
            for n=1:m
                idx = (m-1)*m/2+n;
                epsmn(m,n,f,q) = o.dat.realepsL(f,q,idx)+1i*o.dat.imagepsL(f,q,idx);
            end
        end

        % enforce symmetry
        % epsmn(:,:,f,q) = picomax.eps.enforcesym(epsmn(:,:,f,q),o.num.QVEC(:,q));

        % fill uppder triangle
        for m=1:NEPS
            for n=1:m-1
                epsmn(n,m,f,q) = conj(epsmn(m,n,f,q));
            end
        end
    end
end
epsmn = reshape(epsmn,[NEPS NEPS NFREQ*NQ]);
[V,D] = picomax.util.eigenshuffle(epsmn);
V = reshape(V,[NEPS NEPS NFREQ NQ]);
D = reshape(D,[NEPS NFREQ NQ]);

tol = 1e-10;
for q=1:NQ
    idx = zeros(1,NEPS);
    for m=1:NEPS
        idx(m) = find(abs(D(:,1,q)-D0(m,q))<tol,1);
    end
    D(:,:,q) = D(idx,:,q);
    V(:,:,:,q) = V(:,idx,:,q);
end
o.dat.alphaL = permute(D,[2 3 1]);
o.tmp.VL = reshape(V,[NEPS NEPS NFREQ NQ]);
o.tmp.DL = reshape(D,[NEPS NFREQ NQ]);



%-- construct [[eps]] matrix
f = 1;
epsmn = zeros([NEPS NEPS NQ]); % f = 0
for q=1:NQ
    % fill lower triangle
    for m=1:NEPS
        for n=1:m
            idx = (m-1)*m/2+n;
            epsmn(m,n,q) = o.dat.realepsT(f,q,idx)+1i*o.dat.imagepsT(f,q,idx);
        end
    end

    % enforce symmetry
    % epsmn(:,:,q) = picomax.eps.enforcesym(epsmn(:,:,q),o.num.QVEC(:,q));

    % fill uppder triangle
    for m=1:NEPS
        for n=1:m-1
            epsmn(n,m,q) = conj(epsmn(m,n,q));
        end
    end
end
[V0,D0] = picomax.util.eigenshuffle(epsmn);

epsmn = zeros([NEPS NEPS NFREQ NQ]);
for q=1:NQ
    for f=1:NFREQ
        % fill lower triangle
        for m=1:NEPS
            for n=1:m
                idx = (m-1)*m/2+n;
                epsmn(m,n,f,q) = o.dat.realepsT(f,q,idx)+1i*o.dat.imagepsT(f,q,idx);
            end
        end

        % enforce symmetry
        % epsmn(:,:,f,q) = picomax.eps.enforcesym(epsmn(:,:,f,q),o.num.QVEC(:,q));

        % fill uppder triangle
        for m=1:NEPS
            for n=1:m-1
                epsmn(n,m,f,q) = conj(epsmn(m,n,f,q));
            end
        end
    end
end

epsmn = reshape(epsmn,[NEPS NEPS NFREQ*NQ]);
[V,D] = picomax.util.eigenshuffle(epsmn);
V = reshape(V,[NEPS NEPS NFREQ NQ]);
D = reshape(D,[NEPS NFREQ NQ]);

tol = 1e-10;
for q=1:NQ
    idx = zeros(1,NEPS);
    for m=1:NEPS
        idx(m) = find(abs(D(:,1,q)-D0(m,q))<tol,1);
    end
    D(:,:,q) = D(idx,:,q);
    V(:,:,:,q) = V(:,idx,:,q);
end
o.dat.alphaT = permute(reshape(D,[NEPS NFREQ NQ]),[2 3 1]);
o.tmp.VT = reshape(V,[NEPS NEPS NFREQ NQ]);
o.tmp.DT = reshape(D,[NEPS NFREQ NQ]);
end

