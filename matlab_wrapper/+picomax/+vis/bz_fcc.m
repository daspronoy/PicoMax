function [v,f,e] = bz_fcc()




% primitive lattice vectors
% a = 5.43;
a1 = [0;1;1];
a2 = [1;0;1];
a3 = [1;1;0];

% reciprocal lattice vectors
b1 = 2*cross(a2,a3)/dot(a1,cross(a2,a3));
b2 = 2*cross(a3,a1)/dot(a2,cross(a3,a1));
b3 = 2*cross(a1,a2)/dot(a3,cross(a1,a2));


hkl = [
    0 0 1
    0 0 -1
    0 1 0
    0 -1 0
    0 1 1
    0 -1 -1
    1 0 0
    -1 0 0
    1 0 1
    -1 0 -1
    1 1 0
    -1 -1 0
    1 1 1
    -1 -1 -1
    1 1 -1
    -1 -1 1
    1 -1 1
    -1 1 -1
    -1 1 1
    1 -1 -1
    -1 1 0
    1 -1 0
    1 0 -1
    -1 0 1
    0 1 -1
    0 -1 1];

G = zeros(26,4);
for i=1:26
    G(i,1:3) = b1*hkl(i,1)+b2*hkl(i,2)+b3*hkl(i,3);
    G(i,4) = 1;
    d_gamma = norm(G(i,1:3))/2; % distance to Gamma
    for j=1:14
        if j~=i
            GG = norm(G(i,1:3)/2-G(j,1:3));
            if GG<d_gamma
                G(i,4) = 0;
                break;
            end
        end
    end
end

%-- intersection corners
% c = [x,y,z,i,j,k]
Np = sum(G(:,4)); % number of planes
c = [];
for i=3:Np
    for j=2:i
        for k=1:j
            A = [G(i,1:3);G(j,1:3);G(k,1:3)];
            b = sum(A.^2,2)/2;

            x = det(A);
            if x~=0
                x = A\b;
                c = [c;[x.',i,j,k]];
            end
        end
    end
end

for i=size(c,1):-1:1
    d_gamma = norm(c(i,1:3)); % distance to Gamma
    for j=1:26
        if j~=i
            DGG = norm(c(i,1:3)-G(j,1:3));
            if DGG<d_gamma
                c(i,:) = [];
                break;
            end
        end
    end
end

Nc = sum(c(:,4)); % number of corners


v = c(:,1:3); % vertices [3 x Nv]

%-- find edge information
% e = [e1x e1y e1z e2x e2y e2z p q]
% {p,q} index of planes associated with the edge
% {e1x,e1y,e1z} start coordinate
% {e2x,e2y,e2z} end coordinate
e = zeros(8,0);
for p=1:Np
    % index of edges corresponding to the face $p$
    i1 = unique([find(c(:,4)==p);find(c(:,5)==p);find(c(:,6)==p)]);
    for q=p+1:Np
        i2 = unique([find(c(:,4)==q);find(c(:,5)==q);find(c(:,6)==q)]);
        ic = intersect(i1,i2);
        if ~isempty(ic)
            e = [e,[c(ic(1),1:3) c(ic(2),1:3) p q].'];
        end
    end

end
e = e.';

% face index
f = nan(6,Np);
for p=1:Np
    % index of edges associated with the face $p$
    % idx_face_p = unique([find(c(:,4)==p);find(c(:,5)==p);find(c(:,6)==p)]);
    
    idx_face_p = unique([find(e(:,7)==p);find(e(:,8)==p)]);


    % coordinates of edges associated with the face $p$
    Ne = numel(idx_face_p);
    edge_face_p = e(idx_face_p,1:6).';

    % edge_face_p = zeros(6,0);
    idx_edge = nan(1,Ne);
    idx_edge(1) = find(vecnorm(v.'-edge_face_p(1:3,1))<1e-6,1);
    idx_edge(2) = find(vecnorm(v.'-edge_face_p(4:6,1))<1e-6,1);
    edge_end = edge_face_p(4:6,1);
    edge_face_p(:,1) = [];
    
    for i=3:Ne
        j = find(vecnorm(edge_face_p(1:3,:)-edge_end)<1e-6,1);
        if ~isempty(j)
            edge_end = edge_face_p(4:6,j);
        else
            j = find(vecnorm(edge_face_p(4:6,:)-edge_end)<1e-6,1);
            edge_end = edge_face_p(1:3,j);
        end
        edge_face_p(:,j) = [];


        idx_edge(i) = find(vecnorm(v.'-edge_end)<1e-6,1);
        
    end
    f(1:Ne,p) = idx_edge(1:Ne);
end

% modify data for plotting
f = f.';
e(:,7:8) = [];
