function [e] = initialize(N)

e = sym([N N]);
for m=1:N
    for n=1:N
        if m<10 && n<10
            e(m,n) = sym(sprintf('e%d%d',m-1,n-1));
        else
            e(m,n) = sym(sprintf('e%d_%d',m-1,n-1));
        end
    end
end

end

