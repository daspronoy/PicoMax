








function [] = band(o,varargin)

%** inputparser
p = inputParser;
addParameter(p,'forcesym',false);
addParameter(p,'tol',1e-10);
addParameter(p,'pol',{'L','T'});
addParameter(p,'originshift',0);
addParameter(p,'sort',true);
parse(p,varargin{:});
p = p.Results;
neps = o.num.neps;
nfreq = o.num.nfreq;
nq = sum(o.num.qnum)+1;
o.sys.initialize
picomax.qvector.generate(o)

%% Longitudinal
if any(strcmp(p.pol,'L'))
    %-- construct [[eps]] matrix
    if p.sort
        f = 1;
        epsmn = zeros([neps neps nq]); % f = 0
        for q=1:nq
            % fill lower triangle
            for m=1:neps
                for n=1:m
                    idx = (m-1)*m/2+n;
                    epsmn(m,n,q) = o.dat.realepsL(f,q,idx)+1i*o.dat.imagepsL(f,q,idx);
                end
            end

            % enforce symmetry
            if p.forcesym
                epsmn(:,:,q) = picomax.eps.enforcesym(epsmn(:,:,q),o.num.qvec(:,q));
            end

            % fill uppder triangle
            for m=1:neps
                for n=1:m-1
                    epsmn(n,m,q) = conj(epsmn(m,n,q));
                end
            end
        end
        [~,D0] = picomax.util.eigenshuffle(epsmn);
    end


    epsmn = zeros([neps neps nfreq nq]);
    for q=1:nq
        for f=1:nfreq
            % fill lower triangle
            for m=1:neps
                for n=1:m
                    idx = (m-1)*m/2+n;
                    epsmn(m,n,f,q) = o.dat.realepsL(f,q,idx)+1i*o.dat.imagepsL(f,q,idx);
                end
            end

            % enforce symmetry
            if p.forcesym
                epsmn(:,:,f,q) = picomax.eps.enforcesym(epsmn(:,:,f,q),o.num.qvec(:,q));
            end

            % fill uppder triangle
            for m=1:neps
                for n=1:m-1
                    epsmn(n,m,f,q) = conj(epsmn(m,n,f,q));
                end
            end
        end
    end
    if p.sort
        epsmn = reshape(epsmn,[neps neps nfreq*nq]);
        [V,D] = picomax.util.eigenshuffle(epsmn);
        V = reshape(V,[neps neps nfreq nq]);
        D = reshape(D,[neps nfreq nq]);
        for q=1:nq
            idx = zeros(1,neps);
            for m=1:neps
                idx(m) = find(abs(D(:,1,q)-D0(m,q))<p.tol,1);
            end
            D(:,:,q) = D(idx,:,q);
            V(:,:,:,q) = V(:,idx,:,q);
        end
    else
        V = zeros([neps neps nfreq nq]);
        D = zeros([neps nfreq nq]);
        for q=1:nq
            for f=1:nfreq
                [V(:,:,f,q),D1] = eig(epsmn(:,:,f,q));
                D(:,f,q) = diag(D1);
            end
        end
    end
    o.dat.alphaL = permute(D,[2 3 1]);
    o.tmp.VL = reshape(V,[neps neps nfreq nq]);
    o.tmp.DL = reshape(D,[neps nfreq nq]);
end



%% Transverse
if any(strcmp(p.pol,'T'))
    %-- construct [[eps]] matrix
    if p.sort
        f = 1;
        epsmn = zeros([neps neps nq]); % f = 0
        for q=1:nq
            % fill lower triangle
            for m=1:neps
                for n=1:m
                    idx = (m-1)*m/2+n;
                    epsmn(m,n,q) = o.dat.realepsT(f,q,idx)+1i*o.dat.imagepsT(f,q,idx);
                end
            end

            % enforce symmetry
            if p.forcesym
                epsmn(:,:,q) = picomax.eps.enforcesym(epsmn(:,:,q),o.num.qvec(:,q));
            end

            % fill uppder triangle
            for m=1:neps
                for n=1:m-1
                    epsmn(n,m,q) = conj(epsmn(m,n,q));
                end
            end
        end
        [~,D0] = picomax.util.eigenshuffle(epsmn);
    end

    epsmn = zeros([neps neps nfreq nq]);
    for q=1:nq
        for f=1:nfreq
            % fill lower triangle
            for m=1:neps
                for n=1:m
                    idx = (m-1)*m/2+n;
                    epsmn(m,n,f,q) = o.dat.realepsT(f,q,idx)+1i*o.dat.imagepsT(f,q,idx);
                end
            end

            % enforce symmetry
            if p.forcesym
                epsmn(:,:,f,q) = picomax.eps.enforcesym(epsmn(:,:,f,q),o.num.qvec(:,q));
            end

            % fill uppder triangle
            for m=1:neps
                for n=1:m-1
                    epsmn(n,m,f,q) = conj(epsmn(m,n,f,q));
                end
            end
        end
    end
    if p.sort
        epsmn = reshape(epsmn,[neps neps nfreq*nq]);
        [V,D] = picomax.util.eigenshuffle(epsmn);
        V = reshape(V,[neps neps nfreq nq]);
        D = reshape(D,[neps nfreq nq]);

        for q=1:nq
            idx = zeros(1,neps);
            for m=1:neps
                idx(m) = find(abs(D(:,1,q)-D0(m,q))<p.tol,1);
            end
            D(:,:,q) = D(idx,:,q);
            V(:,:,:,q) = V(:,idx,:,q);
        end
    else
        V = zeros([neps neps nfreq nq]);
        D = zeros([neps nfreq nq]);
        for q=1:nq
            for f=1:nfreq
                [V(:,:,f,q),D1] = eig(epsmn(:,:,f,q));
                D(:,f,q) = diag(D1);
            end
        end
    end
    o.dat.alphaT = permute(reshape(D,[neps nfreq nq]),[2 3 1]);
    o.tmp.VT = reshape(V,[neps neps nfreq nq]);
    o.tmp.DT = reshape(D,[neps nfreq nq]);
end
end

