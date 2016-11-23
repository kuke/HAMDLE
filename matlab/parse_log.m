%%
% parse the log file
% example: parse_log('../log/20161122T203118.log')
% param: log_file
%        output figure

function parse_log(log_file)
    
Data = textread(log_file,'','headerlines',5);
Order = Data(:, 1);
figure 
VD_HAMDLE = Data(:, 5);
axes1 = axes('Position',[0.09 0.11 0.81 0.815]);
box(axes1,'on');
hold(axes1,'on');
plot(Order, VD_HAMDLE, 'b', 'linewidth', 1.5);
xlabel('Expansion order'); ylabel('VD');
set(gca,'FontSize',15);
xlim([Order(1) Order(end)]);
ylim([0 1]);

figure
READ = Data(:, 2);
REF = Data(:, 3);
ALGO = Data(:, 4);
Total = READ+REF+ALGO;
axes1 = axes('Position',[0.12 0.11 0.79 0.815]);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[5 100]);
box(axes1,'on');
hold(axes1,'on');
plot(Order, ALGO, 'r', 'linewidth', 1.5); hold on;
plot(Order, Total, 'b', 'linewidth', 1.5);
legend('ALGO', 'Total', 'Location','North');
xlim([Order(1) Order(end)]);
xlabel('Expansion order');ylabel('Convergence time/s');
set(gca,'FontSize',15);

end