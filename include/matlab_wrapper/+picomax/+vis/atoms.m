



function [] = atoms(varargin)

%-- input parser
p = inputParser;
addParameter(p,'unitcell',1);
addParameter(p,'FaceColor',[0.4 0.4 0.9; 0.4 0.9 0.4]);
addParameter(p,'Radius',[0.1 0.1]);
addParameter(p,'LineWidth',1.5);
addParameter(p,'LineColor','k');
addParameter(p,'Shift',[0 0 0]);
parse(p,varargin{:});
p = p.Results;

%--
R1 = p.Radius(1);
if numel(p.Radius)>1
    R2 = p.Radius(2);
else
    R2 = R1;
end
FaceColor1 = p.FaceColor(1,1:3);
if size(p.FaceColor,1)>1
    FaceColor2 = p.FaceColor(2,1:3);
else
    FaceColor2 = FaceColor1;
end


dx = p.Shift(1);
dy = p.Shift(2);
dz = p.Shift(3);
if p.unitcell==1
    x1 = [0 0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1 1];
    y1 = [0 0 0.5 1 1 0 0.5 0.5 1 0 0 0.5 1 1];
    z1 = [0 1 0.5 0 1 0.5 0 1 0.5 0 1 0.5 0 1];
    tau = [1/4 1/4 1/4];
    x2 = [0 0 0.5 0.5]+tau(1);
    y2 = [0 0.5 0 0.5]+tau(2);
    z2 = [0.5 0 0 0.5]+tau(3);
    connection = [ ...
        2 3 6 8
        3 4 7 9
        6 7 10 12
        8 9 12 14];
    [xs,ys,zs] = sphere(200);

    % atom 1
    for i=1:numel(x1)
        surf(R1*xs+x1(i)+dx,R1*ys+y1(i)+dy,R1*zs+z1(i)+dz, ...
            'EdgeColor','none','FaceColor',FaceColor1)
    end
    % atom 2
    for i=1:numel(x2)
        surf(R2*xs+x2(i)+dx,R2*ys+y2(i)+dy,R2*zs+z2(i)+dz, ...
            'EdgeColor','none','FaceColor',FaceColor2)
    end
    for i=1:numel(x2)
        for j=1:4
            line([x2(i) x1(connection(i,j))]+dx, ...
                [y2(i) y1(connection(i,j))]+dy, ...
                [z2(i) z1(connection(i,j))]+dz, ...
                'Color',p.LineColor,'LineWidth',p.LineWidth);
        end
    end
else

end





