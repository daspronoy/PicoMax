%picomax.sym.reflection | returns the symmetry operation matrix for
%reflection w.r.t. \vec{u}


function [S] = reflection(u,gvec)

tol = 1e-10;
N = size(gvec,2);
u = u(:)/norm(u);
R = eye(3)-2*u*u.';
gvec_transform = R*gvec; % reflected gvec


S = zeros(N);
for i=1:N
    tmp = vecnorm(gvec-gvec_transform(:,i));
    S(i,tmp<tol) = 1;
end

if det(S)==0
    error('determinant of the symmetry operation matrix is zero');
end