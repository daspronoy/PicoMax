% returns tangential vectors in certain way...





function [t] = tangential(n)



nx = n(1,:);
ny = n(2,:);
nz = n(3,:);
N = numel(nx);



%% modified (100)
% t = zeros(size(n));
% for i=1:N
%     if nz(i)~=0
%         t(:,i) = [0;nz(i);-ny(i)];
%     elseif nz(i)==0 && ny(i)~=0
%         t(:,i) = [-ny(i);nx(i);0];
%     elseif nz(i)==0 && ny(i)==0
%         t(:,i) = [0;1;0];
%     else
%         error('hmmf');
%     end
% end
% t = t./vecnorm(t);


%% modified (010) ??
t = zeros(size(n));
for i=1:N
    if nz(i)~=0
        t(:,i) = [nz(i);0;-nx(i)];
    elseif nz(i)==0 && ny(i)~=0
        t(:,i) = [-ny(i);nx(i);0];
    elseif nz(i)==0 && ny(i)==0
        t(:,i) = [0;1;0];
    else
        error('hmmf');
    end
end
t = t./vecnorm(t);





% %% 2024.03.09
% t = zeros(size(n));
% for i=1:N
%     if ny(i)==0 && nz(i)==0
%         t(:,i) = [0;1;0];
%     else
%         if nx(i)>0
%             tx = 1;
%             ty = -1/2/nz(i)-nx(i)/2/ny(i);
%             tz = 1/2/ny(i)-nx(i)/2/nz(i);
%         else
%             tx = -1;
%             ty = 1/2/nz(i)+nx(i)/2/ny(i);
%             tz = -1/2/ny(i)+nx(i)/2/nz(i);
%             % tx = -tx;
%             % ty = -ty;
%             % tz = -tz;
% 
%             % tx = -1;
%             % ty = -1/2/nz(i)-nx(i)/2/ny(i);
%             % tz = 1/2/ny(i)-nx(i)/2/nz(i);
%         end
% 
%         t(:,i) = [tx;ty;tz];
%     end
% end
% t = t./vecnorm(t);






% % G-X (100)?
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
%         t(:,i) = [0;1;0];
%     else
%         error('hmmf');
%     end
% end
% t = t./vecnorm(t);