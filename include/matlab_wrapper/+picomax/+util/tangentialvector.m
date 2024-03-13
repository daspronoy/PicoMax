



function [t] = tangentialvector(n)

if numel(n)~=3
    error('tangentialvector::input should be a vector');
end
t = zeros(size(n));

nx = n(1);
ny = n(2);
nz = n(3);

if ny~=0
    t = [-ny nx 0];
    t = t/norm(t);
elseif ny==0 && nz~=0
    t = [nz 0 -nx];
    t = t/norm(t);
elseif ny==0 && nz==0
    t = [0 1 0];
else
    error('tangentialvector::check input vector');
end





