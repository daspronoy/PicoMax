


function [k] = generate_fcc(N)

err = 1e-10;


df = 1/N;
[f1,f2,f3] = ndgrid(0:df:1-df,0:df:1-df,0:df:1-df);

% transform to cartesian coordinate
b1 = [-1 1 1];
b2 = [1 -1 1];
b3 = [1 1 -1];
k = b1.*f1(:)+b2.*f2(:)+b3.*f3(:);

% translate to lower |k|
[g1,g2,g3] = ndgrid(-2:2,-2:2,-2:2);
i000 = find(g1==0 & g2==0 & g3==0);
g1(i000) = [];
g2(i000) = [];
g3(i000) = [];
Ng = numel(g1);

Nk = size(k,1);
for i=1:Nk
    ki = k(i,:);
    normki = norm(ki);
    kj = ki+b1.*g1(:)+b2.*g2(:)+b3.*g3(:);
    normkj = zeros(Ng,1);
    for j=1:Ng
        normkj(j) = norm(kj(j,:));
    end
    [minnormkj,idx] = min(normkj);
    if normki>minnormkj
        k(i,:) = kj(idx,:);
    end
end

% inside 1st BZ
for i=Nk:-1:1
    if (abs(k(i,1))+abs(k(i,2))+abs(k(i,3))-1.5>err) ... [111] direction
            || (abs(k(i,1))-1.0>err) ... [100] direction
            || (abs(k(i,2))-1.0>err) ... [010]
            || (abs(k(i,3))-1.0>err)    %[001]
        disp('hmmf')
    end
end