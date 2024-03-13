











function [kx,ky,kz,w] = cohenchadi(n)

Nk = 2^(n-1)*(2^(n-1)+1)*(2^n+1)/3;


if n==1
    kx = [0.75 0.25];
    ky = [0.25 0.25];
    kz = [0.25 0.25];
    w = [0.75 0.25];
    return
end

if n==2
    kx = [0.875 0.875 0.625 0.625 0.625 0.625 0.375 0.375 0.375 0.125];
    ky = [0.375 0.125 0.625 0.375 0.375 0.125 0.375 0.375 0.125 0.125];
    kz = [0.125 0.125 0.125 0.375 0.125 0.125 0.375 0.125 0.125 0.125];
    w = [0.18750 0.09375 0.09375 0.09375 0.18750 0.09375 0.03125 0.09375 0.09375 0.03125];
    return;
end


kx = zeros(1,Nk);
ky = zeros(1,Nk);
kz = zeros(1,Nk);
w = zeros(1,Nk);
k = 1;
for i=1:1:2^(n-1)
    for j=1:i
        for l=1:j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=2^(n-1)+2:2:2^n
    for j=1:3*2^(n-2)-0.5*i
        for l=1:j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=2^(n-2)+2:2:3*2^(n-2)
    for j=3*2^(n-2)-0.5*i+1:i
        for l=1:3*2^(n-1)+1-i-j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=3*2^(n-2)+2:2:2^n
    for j=3*2^(n-2)-0.5*i+1:3*2^(n-1)-i
        for l=1:3*2^(n-1)+1-i-j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=2^(n-1)+1:2:2^n-1
    for j=1:3*2^(n-2)-0.5*(i-1)
        for l=1:j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=2^(n-1)+1:2:3*2^(n-2)-1
    for j=3*2^(n-2)-0.5*(i-1)+1:i
        for l=1:3*2^(n-1)+1-i-j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end
for i=3*2^(n-2)+1:2:2^n-1
    for j=3*2^(n-2)-0.5*(i-1)+1:3*2^(n-1)-i
        for l=1:3*2^(n-1)+1-i-j
            kx(k) = (2*i-1)/2^(n+1);
            ky(k) = (2*j-1)/2^(n+1);
            kz(k) = (2*l-1)/2^(n+1);
            w(k) = (3*(1-delta(i,j))+3*(1-delta(j,l))+delta(i,l))/(4*8^(n-1));
            k = k+1;
        end
    end
end

fprintf('Nk = %d\n',Nk);
fprintf('Nk = %d\n',numel(kx));
fprintf('weight sum = %f\n',sum(w));
function [d] = delta(i,j)
if i==j
    d = 1;
else
    d = 0;
end
end
end