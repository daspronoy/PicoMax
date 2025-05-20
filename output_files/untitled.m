filename = ['temp2.dat'];
fid = fopen(filename);
data= fread(fid,'double');
fclose(fid);

npts=41;
eband=reshape(data(1:end,:),[40 npts]);
set(gcf,'renderer','opengl', 'Position',  [500, 100, 900, 800])

plot(1:npts,eband(1:end,:), 'Color','#505050','LineWidth',1.5);

% xticks([1 21 61 81 121 141 181]);
% xticklabels({"A","$\Gamma$", "K", "M", "L", "H","A"});
% % xticklabels({"K", "L", "H", "$\Gamma$", "A", "M"});

ylabel("$\epsilon - \epsilon_F$ (eV)","Interpreter","latex");
% xline([21 61 81 121 141],'--');
ylim([-3 1]);
xlim([1 npts]);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 20, 'LineWidth', 0.5);