


function [k] = generate_fcc(N)

tol = 1e-10;

% reciprocal lattice vectors (in units of 2*pi/a)
b1 = [-1 1 1].';
b2 = [1 -1 1].';
b3 = [1 1 -1].';


% generate fractional coordinate grids
df = 1/N;
fi = 0:df:1;
[f1,f2,f3] = ndgrid(fi,fi,fi);
f1 = f1(:).';
f2 = f2(:).';
f3 = f3(:).';

% transform to cartesian coordinate
k = b1.*f1+b2.*f2+b3.*f3; % [3 x N^3]

% translate by a g-vector, so that |k| is minimized
[g1,g2,g3] = ndgrid(-1:1,-1:1,-1:1);
i000 = find(g1==0 & g2==0 & g3==0);
g1(i000) = [];
g2(i000) = [];
g3(i000) = [];
Nk = size(k,2);
for i=1:Nk
    ki = k(:,i);
    kj = ki+b1.*g1+b2.*g2+b3.*g3;
    normkj = vecnorm(kj);
    [minnormkj,idx] = min(normkj);
    if norm(ki)>minnormkj
        k(:,i) = kj(:,idx);
    end
end

% check if all the points are inside the 1st BZ
for i=1:Nk
    if (abs(k(1,i))+abs(k(2,i))+abs(k(3,i))-1.5>tol) ... [111] direction
            || (abs(k(1,i))-1.0>tol) ... [100] direction
            || (abs(k(2,i))-1.0>tol) ... [010]
            || (abs(k(3,i))-1.0>tol)    %[001]
        disp('generate_fcc :: a point outside the 1st BZ was found');
    end
end

% % duplicate edge points
% for i=1:Nk
%     if (abs(abs(k(1,i))+abs(k(2,i))+abs(k(3,i))-1.5)<tol) ... [111] direction
%             || (abs(abs(k(1,i))-1.0)<tol) ... [100] direction
%             || (abs(abs(k(2,i))-1.0)<tol) ... [010]
%             || (abs(abs(k(3,i))-1.0)<tol)    %[001]
%         ki = k(:,i);
%         kj = ki+b1.*g1+b2.*g2+b3.*g3;
% 
%     end
% end









