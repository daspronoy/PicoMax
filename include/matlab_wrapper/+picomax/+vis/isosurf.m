



function [h] = isosurf(x,y,z,s,f,varargin)

p = inputParser;
addParameter(p,'FaceAlpha',0.5);
addParameter(p,'FaceColor',[1.0 0.7 0.7]);
parse(p,varargin{:});
p = p.Results;


mins = min(s(:));
maxs = max(s(:));
v = mins+f*(maxs-mins);

[f_,v_] = isosurface(x,y,z,s,v);
h = patch('Faces',f_,'Vertices',v_, ...
    'FaceColor',p.FaceColor, ...
    'EdgeColor','none', ...
    'FaceAlpha',p.FaceAlpha);

