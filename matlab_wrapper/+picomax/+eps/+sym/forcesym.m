function [e] = forcesym(e)

N = size(e,1);

for m=1:N
    for n=1:m-1
        e(m,n) = e(n,m);
    end
end
end

