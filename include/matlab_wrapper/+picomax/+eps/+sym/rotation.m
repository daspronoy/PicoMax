function [e,S] = rotation(e,u,t,gvec,skip)
if nargin<5
    skip = false;
end
N = size(gvec,2);
u = u(:)/norm(u); % rotation axis
t = t*pi/180; % rotation angle [rad]

% rotation matrix
R = [cos(t)+u(1)^2*(1-cos(t)), u(1)*u(2)*(1-cos(t))-u(3)*sin(t), u(1)*u(3)*(1-cos(t))+u(2)*sin(t);
    u(2)*u(1)*(1-cos(t))+u(3)*sin(t), cos(t)+u(2)^2*(1-cos(t)), u(2)*u(3)*(1-cos(t))-u(1)*sin(t);
    u(3)*u(1)*(1-cos(t))-u(2)*sin(t), u(3)*u(2)*(1-cos(t))+u(1)*sin(t), cos(t)+u(3)^2*(1-cos(t))];
gvec_transform = R*gvec; % rotated gvec
S = zeros(N); % symmetry operation matrix
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
