filename = ['temp3.dat'];
fid = fopen(filename);
data= fread(fid,'double');
fclose(fid);

eband=reshape(data(1:end,:),[40 121]);
set(gcf,'renderer','opengl', 'Position',  [500, 100, 500, 500])

plot(1:121,eband(1:end,:)-0.895, 'Color','#505050','LineWidth',1.5);

xticks([1 21 61 81 121]);
xticklabels({"$\Gamma$", "A", "H", "K", "$\Gamma$"});

ylabel("Energy (eV)","Interpreter","latex");

ylim([-2 2]);
xlim([1 121]);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 20, 'LineWidth', 0.5);