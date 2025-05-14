



function [h] = isosurf(x,y,z,s,f,varargin)

p = inputParser;
addParameter(p,'FaceAlpha',0.5);
addParameter(p,'FaceColor',[1.0 0.7 0.7]);
addParameter(p,'sym',false);
addParameter(p,'complex',0);
addParameter(p,'arg',false);
addParameter(p,'phase',0);
addParameter(p,'abs',0);
parse(p,varargin{:});
p = p.Results;

if p.phase
    s = s*exp(1i*p.phase*pi/180);
end

if p.arg
    s = angle(s);
elseif p.complex==1 && p.sym
    maxs = max(abs(s(:)));
    mins = -maxs;
    s = real(s);
elseif p.complex==2 && p.sym
    maxs = max(abs(s(:)));
    mins = -maxs;
    s = imag(s);
elseif p.complex==0 && ~p.sym
    mins = min(s(:));
    maxs = max(s(:));
elseif p.complex==0 && p.sym
    mins = min(s(:));
    maxs = max(s(:));
    maxs = max([-mins maxs]);
    mins = -maxs;
end

if p.arg
    v = 0;
elseif p.abs
    v = f*maxs;
else
    v = mins+f*(maxs-mins);
end



[f_,v_] = isosurface(x,y,z,s,v);
h = patch('Faces',f_,'Vertices',v_, ...
    'FaceColor',p.FaceColor, ...
    'EdgeColor','none', ...
    'FaceAlpha',p.FaceAlpha);

