function [h] = cone3(fx,fy,fz,x_,y_,z_,varargin)

fx = permute(fx,[2 1 3]);
fy = permute(fy,[2 1 3]);
fz = permute(fz,[2 1 3]);


x_ = permute(x_,[2 1 3]);
y_ = permute(y_,[2 1 3]);
z_ = permute(z_,[2 1 3]);


h = coneplot(x_,y_,z_,fx,fy,fz,0.075,'nointerp');

set(h,'EdgeColor','none')
set(h,'FaceColor','red')
set(h,'AmbientStrength',0.6)
set(h,'DiffuseStrength',0.8)
set(h,'FaceAlpha',0.6)
