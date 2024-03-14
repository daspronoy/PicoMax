%picomax.sym.rotation | returns the symmetry operation matrix foQNUMr
%rotation by t (degrees) w.r.t. \vec{u}


function [S] = rotation(u,t,gvec)

tol = 1e-10;
N = size(gvec,2);
u = u(:)/norm(u); % rotation axis
cos_t = cosd(t);
sin_t = sind(t);


% rotation matrix
R = [cos_t+u(1)^2*(1-cos_t), u(1)*u(2)*(1-cos_t)-u(3)*sin_t, u(1)*u(3)*(1-cos_t)+u(2)*sin_t;
    u(2)*u(1)*(1-cos_t)+u(3)*sin_t, cos_t+u(2)^2*(1-cos_t), u(2)*u(3)*(1-cos_t)-u(1)*sin_t;
    u(3)*u(1)*(1-cos_t)-u(2)*sin_t, u(3)*u(2)*(1-cos_t)+u(1)*sin_t, cos_t+u(3)^2*(1-cos_t)];
gvec_transform = R*gvec; % rotated gvec

% symmetry operation matrix
S = zeros(N); 
for i=1:N
    tmp = sum(abs(gvec-gvec_transform(:,i)).^2,1);
    S(i,tmp<tol) = 1;
end

% check
if det(S)==0
    error('determinant of the symmetry operation matrix is zero');
end