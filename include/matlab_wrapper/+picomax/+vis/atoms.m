



function [] = atoms(varargin)

%-- input parser
p = inputParser;
addParameter(p,'UnitCell',[1 1 1]);
addParameter(p,'LatticeVector',[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0].');
addParameter(p,'Tau',[0 0 0; 0.25 0.25 0.25].');
addParameter(p,'FaceColor',[0.4 0.4 0.9; 0.4 0.9 0.4]);
addParameter(p,'FaceAlpha',1);
addParameter(p,'Radius',[0.1 0.1]);
addParameter(p,'LineWidth',1.5);
addParameter(p,'LineColor','k');
addParameter(p,'Origin',[]);
addParameter(p,'AmbientStrength',0.7);
addParameter(p,'DiffuseStrength',0.6);
addParameter(p,'SpecularStrength',0.9);
addParameter(p,'SpecularExponent',10);
addParameter(p,'SpecularColorReflectance',1);
addParameter(p,'res',20)
parse(p,varargin{:});
p = p.Results;

if isempty(p.Origin)
    p.Origin = 0.5*p.UnitCell.';
end


% %--
% R1 = p.Radius(1);
% if numel(p.Radius)>1
%     R2 = p.Radius(2);
% else
%     R2 = R1;
% end
% FaceColor1 = p.FaceColor(1,1:3);
% if size(p.FaceColor,1)>1
%     FaceColor2 = p.FaceColor(2,1:3);
% else
%     FaceColor2 = FaceColor1;
% end

%-- generate atom coordinates
NATOM = size(p.Tau,2);
r = cell(1,NATOM);
NPOS = zeros(1,NATOM);
for a=1:NATOM
    r{a} = zeros(3,prod(p.UnitCell*4+1));
    cnt = 1;
    for i=-2*p.UnitCell(1):2*p.UnitCell(1)
        for j=-2*p.UnitCell(2):2*p.UnitCell(2)
            for k=-2*p.UnitCell(3):2*p.UnitCell(3)
                pos = p.LatticeVector*[i;j;k] + p.Tau(:,a);
                if abs(pos(1)-p.Origin(1))<=p.UnitCell(1)/2 ...
                        && abs(pos(2)-p.Origin(2))<=p.UnitCell(2)/2 ...
                        && abs(pos(3)-p.Origin(3))<=p.UnitCell(3)/2
                    r{a}(:,cnt) = pos;
                    cnt = cnt+1;
                end
            end
        end
    end
    r{a}(:,cnt:end) = [];
    NPOS(a) = cnt-1;
end


% plot atoms
[xs,ys,zs] = sphere(p.res);
for a=1:NATOM
    for i=1:NPOS(a)
        surf(p.Radius(a)*xs+r{a}(1,i),p.Radius(a)*ys+r{a}(2,i),p.Radius(a)*zs+r{a}(3,i), ...
            'EdgeColor','none','FaceColor',p.FaceColor(a,:),'FaceAlpha',p.FaceAlpha, ...
            'AmbientStrength',p.AmbientStrength, ...
            'DiffuseStrength',p.DiffuseStrength, ...
            'SpecularStrength',p.SpecularStrength, ...
            'SpecularExponent',p.SpecularExponent, ...
            'SpecularColorReflectance',p.SpecularColorReflectance);
    end
end

% connection
for a=1:NATOM
    for b=a+1:NATOM
        for i=1:NPOS(a)
            for j=1:NPOS(b)
                vec = r{a}(:,i)-r{b}(:,j);
                dist = norm(vec);
                if dist<0.25*sqrt(3)+eps
                    line([r{a}(1,i) r{b}(1,j)], ...
                        [r{a}(2,i) r{b}(2,j)], ...
                        [r{a}(3,i) r{b}(3,j)], ...
                        'Color',p.LineColor,'LineWidth',p.LineWidth);
                end
            end
        end
    end
end






% 
% 
% 
% 
% 
% 
% 
% dx = p.Origin(1);
% dy = p.Origin(2);
% dz = p.Origin(3);
% 
% x1 = [0 0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1 1];
% y1 = [0 0 0.5 1 1 0 0.5 0.5 1 0 0 0.5 1 1];
% z1 = [0 1 0.5 0 1 0.5 0 1 0.5 0 1 0.5 0 1];
% tau = [1/4 1/4 1/4];
% x2 = [0 0.5 0.5 0]+tau(1);
% y2 = [0 0.5 0   0.5]+tau(2);
% z2 = [0 0   0.5 0.5]+tau(3);
% connection = [ ...
%     2 3 6 8
%     3 4 7 9
%     6 7 10 12
%     8 9 12 14];
% 
% 
% % atom 1
% for i=1:numel(x1)
%     surf(R1*xs+x1(i)+dx,R1*ys+y1(i)+dy,R1*zs+z1(i)+dz, ...
%         'EdgeColor','none','FaceColor',FaceColor1)
% end
% % atom 2
% for i=1:numel(x2)
%     surf(R2*xs+x2(i)+dx,R2*ys+y2(i)+dy,R2*zs+z2(i)+dz, ...
%         'EdgeColor','none','FaceColor',FaceColor2)
% end
% for i=1:numel(x2)
%     for j=1:numel(x1)
%         dist = norm([x1(j)-x2(i);y1(j)-y2(i);z1(j)-z2(i)]);
%         if dist<0.25*sqrt(3)+eps
%             line([x2(i) x1(j)]+dx, ...
%                 [y2(i) y1(j)]+dy, ...
%                 [z2(i) z1(j)]+dz, ...
%                 'Color',p.LineColor,'LineWidth',p.LineWidth);
%         end
%     end
% end





