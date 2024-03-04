% returns tangential vectors in certain way...





function [t] = tangential(n)



nx = n(1,:);
ny = n(2,:);
nz = n(3,:);
N = numel(nx);




% G-X (100)?
t = zeros(size(n));
for i=1:N
    if ny(i)==0 && nz(i)==0
        t(:,i) = [0;1;0];
    else
        
        tx = 2*(nx(i)^2+1)/nx(i);
        ty = -(nx(i)^2+1);
        t(:,i) = [tx;ty;ty];
    end
end
t = t./vecnorm(t);






% G-X (100)?
% t = zeros(size(n));
% for i=1:N
%     if ny(i)==0 && nz(i)==0
%         t(:,i) = [0;1;0];
%     elseif abs(ny(i)-nz(i))<1e-10
%         t(:,i) = [1;-nx(i)/(ny(i)+nz(i));-nx(i)/(ny(i)+nz(i))];
%     elseif abs(ny(i)+nz(i))<1e-10
%         t(:,i) = [1;-nx(i)/(ny(i)-nz(i));+nx(i)/(ny(i)-nz(i))];
%     end
% end
% t = t./vecnorm(t);









% t = zeros(size(n));
% for i=1:N
%     if nz(i)~=0
%         t(:,i) = [0;nz(i);-ny(i)];
%     elseif nz(i)==0 && ny(i)~=0
%         t(:,i) = [-ny(i);nx(i);0];
%     elseif nz(i)==0 && ny(i)==0
%         t(:,i) = [0;0;-1];
%     else
%         error('hmmf');
%     end
% end
% t = t./vecnorm(t);






% t = zeros(size(n));
% for i=1:N
%     if nz(i)~=0
%         t(:,i) = [-nz(i);0;nx(i)];
%     elseif nz(i)==0 && ny(i)~=0
%         t(:,i) = [ny(i);-nx(i);0];
%     elseif nz(i)==0 && ny(i)==0
%         t(:,i) = [0;1;0];
%     else
%         error('hmmf');
%     end
% end
% t = t./vecnorm(t);



% %-- Sathwik' choice
% t = zeros(size(n));
% for i=1:N
%     if ny(i)~=0
%         t(:,i) = [-ny(i);nx(i);0];
%         % t(:,i) = t(:,i)/norm(t(:,i));
%     elseif ny(i)==0 && nz(i)~=0
%         t(:,i) = [nz(i);0;-nx(i)];
%         % t(:,i) = t(:,i)/norm(t(:,i));
%     elseif ny(i)==0 && nz(i)==0
%         t(:,i) = [0;0;1];
%     else
%         error('hmmf');
%     end
% end
% t = t./vecnorm(t);