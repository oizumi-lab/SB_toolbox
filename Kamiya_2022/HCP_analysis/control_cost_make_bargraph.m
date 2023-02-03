function [f,lgd] = control_cost_make_bargraph(tasks_included_vec, mvt_btstrp)

clear ctr; clear ydt;

all_task_list =...
    {'REST1', 'REST2', 'EMOTION', 'GAMBLING', 'LANGUAGE',...
    'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM'};
tasks_label = all_task_list(tasks_included_vec);

yy = mvt_btstrp(tasks_included_vec,:,:);
er = std(yy,0,3); 

%% 
f = figure('DefaultTextFontName', 'Helvetica', 'DefaultAxesFontName', 'Helvetica');
yy_mean = mean(yy,3);
rank_vec = tiedrank(yy_mean(:,3));
order_vec = sortrows([ [1:length( tasks_included_vec)]',rank_vec],2);
order_vec = order_vec(:,1);
yy_mean = yy_mean(order_vec,:);
er = sortrows([order_vec, er]);
er(:,1) = [];

tasks_label_fig = tasks_label(order_vec)';
X = categorical(tasks_label_fig);
X = reordercats(X, tasks_label_fig);
hBar = bar(yy_mean);
for k1 = 1:size(yy,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, er', '.k');
hold off

lgd = legend('Mean control','Covariance control','Total control','location','southoutside','Orientation','horizontal');
lgd.FontSize = 8;

set(gca,'XTick', [1:length( tasks_included_vec)], 'XTickLabel', tasks_label_fig, 'fontsize', 8)
if size(tasks_included_vec,2) > 5
    set(gca, 'XTickLabelRotation',330)
end

set(gcf,'Position',[30 100 750 600])

end