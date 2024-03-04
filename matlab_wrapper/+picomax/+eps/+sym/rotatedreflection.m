function [e,S] = rotatedreflection(e,u_rot,t,u_ref,gvec,skip)

if nargin<7
    skip = false;
end
N = size(gvec,2);


% rotation matrix
u_rot = u_rot(:)/norm(u_rot); % rotation axis
t = t*pi/180; % rotation angle [rad]
R_rot = [cos(t)+u_rot(1)^2*(1-cos(t)), u_rot(1)*u_rot(2)*(1-cos(t))-u_rot(3)*sin(t), u_rot(1)*u_rot(3)*(1-cos(t))+u_rot(2)*sin(t);
    u_rot(2)*u_rot(1)*(1-cos(t))+u_rot(3)*sin(t), cos(t)+u_rot(2)^2*(1-cos(t)), u_rot(2)*u_rot(3)*(1-cos(t))-u_rot(1)*sin(t);
    u_rot(3)*u_rot(1)*(1-cos(t))-u_rot(2)*sin(t), u_rot(3)*u_rot(2)*(1-cos(t))+u_rot(1)*sin(t), cos(t)+u_rot(3)^2*(1-cos(t))];

% reflection matrix
u_ref = u_ref(:)/norm(u_ref);
R_ref = eye(3)-2*u_ref*u_ref.';


% transfromed gvec
gvec_transform = R_ref*R_rot*gvec;

% construct symmetry operation matrix
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

% compare [[eps]] elements
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
