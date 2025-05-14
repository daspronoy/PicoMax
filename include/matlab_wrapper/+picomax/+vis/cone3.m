function [h] = cone3(fx,fy,fz,x_,y_,z_,ratio,varargin)

%** input parser
p = inputParser;
addParameter(p,'FaceAlpha',1.0);
addParameter(p,'FaceColor',[1.0 0.2 0.2]);
addParameter(p,'AmbientStrength',0.6);
addParameter(p,'DiffuseStrength',0.8);
addParameter(p,'sym',false);
addParameter(p,'complex',0);
parse(p,varargin{:});
p = p.Results;

%** range
s = sqrt(abs(fx).^2+abs(fy).^2+abs(fz).^2);
if p.complex==1 && p.sym
    maxs = max(abs(s(:)));
    mins = -maxs;
    fx = ratio/maxs*real(fx);
    fy = ratio/maxs*real(fy);
    fz = ratio/maxs*real(fz);
elseif p.complex==2 && p.sym
    maxs = max(abs(s(:)));
    mins = -maxs;
    fx = ratio*maxs*imag(fx);
    fy = ratio*maxs*imag(fy);
    fz = ratio*maxs*imag(fz);
elseif p.complex==0 && ~p.sym
    mins = min(s(:));
    maxs = max(s(:));
elseif p.complex==0 && p.sym
    mins = min(s(:));
    maxs = max(s(:));
    maxs = max([-mins maxs]);
    mins = -maxs;
end

fx = permute(fx,[2 1 3]);
fy = permute(fy,[2 1 3]);
fz = permute(fz,[2 1 3]);


x_ = permute(x_,[2 1 3]);
y_ = permute(y_,[2 1 3]);
z_ = permute(z_,[2 1 3]);


h = coneplot(x_,y_,z_,fx,fy,fz,0,'nointerp');
set(h,'EdgeColor','none')
set(h,'FaceColor',p.FaceColor)
set(h,'FaceAlpha',p.FaceAlpha)
set(h,'AmbientStrength',p.AmbientStrength)
set(h,'DiffuseStrength',p.DiffuseStrength)
