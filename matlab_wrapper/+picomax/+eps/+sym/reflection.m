function [e,S] = reflection(e,u,gvec,skip)

if nargin<4
    skip = false;
end


N = size(gvec,2);
u = u(:)/norm(u);
R = eye(3)-2*u*u.';

gvec_transform = R*gvec; % reflected gvec
S = zeros(N);
for i=1:N
    tmp = sum(abs(gvec-gvec_transform(:,i)).^2,1);
    S(i,tmp<1e-10) = 1;
end
if det(S)==0
    error('determinant of the symmetry operation matrix is zero');
end

if skip
    e = S;
    return;
end

e_rot = S*e/S;
for m=1:N
    for n=1:N
        s1 = solve(e(m,n)==e_rot(m,n),e(m,n));
        s2 = solve(e(m,n)==e_rot(m,n),e_rot(m,n));
        if s1~=0 && s2~=0
            e = subs(e,s1,s2);
        end
    end
end