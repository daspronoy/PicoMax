filename = ['chi_11_G.dat'];
fid = fopen(filename);
data= fread(fid,'double');
% fclose(fid);

%z11 -> left z12 -> right


nfreq=1001;
dfreq=0.02;
nq=1;
neps=9;
ne = neps*(neps+1)/2;
N = 3*3*ne*nq*nfreq;
realepsij = permute(reshape(data(1:N),[nfreq ne nq 3 3]),[1 3 2 5 4]);
imagepsij = permute(reshape(data(N+1:2*N),[nfreq ne nq 3 3]),[1 3 2 5 4]);



f=20;
chiij = zeros([3*neps 3*neps nq]);   
for q=1:nq
    %lower triangle
    for i=1:3
        for j=1:3
            for m=1:neps
                for n=1:m
                    idx = (m-1)*m/2+n;
                    chiij((i-1)*neps+m,(j-1)*neps+n,q) = realepsij(f,q,idx,i,j) ...
                        +1i*imagepsij(f,q,idx,i,j);
                end
            end
        end
    end
    %upper triangle
    for i=1:3
        for j=1:3
            for m=1:neps
                for n=1:m
                    chiij((j-1)*neps+n,(i-1)*neps+m,q) = conj(chiij((i-1)*neps+m,(j-1)*neps+n,q));
                end
            end
        end
    end
end
epsij = eye(3*neps)+chiij;
[V,D]=eigenshuffle(epsij);
D=sort(D,'descend','ComparisonMethod','real');
