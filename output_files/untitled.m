filename = ['band.dat'];
fid = fopen(filename);
data= fread(fid,'double');
fclose(fid);

npts=31;
eband=reshape(data(1:end,:),[80 npts]);
set(gcf,'renderer','opengl', 'Position',  [500, 100, 900, 800])

plot(1:npts,eband(1:end,:), 'Color','#505050','LineWidth',1.5);

% xticks([1 41 121 161 241 281 361 401]);
% xticklabels({"$\Gamma$","A","L","K","M", "L", "H","A"});
% 
% 
% ylabel("$\epsilon - \epsilon_F$ (eV)","Interpreter","latex");
% xline([41 121 161 241 281 361],'--');
% 
% ylim([-2 1]);
% xlim([1 npts]);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 20, 'LineWidth', 0.5);