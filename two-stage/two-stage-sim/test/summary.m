%fileID = fopen('test7_optObjVal_kappa_14.6207.txt';'r');
%formatSpec = '%f';
%sizeA = [150 1];
%A = fscanf(fileID;formatSpec;sizeA);
%mean(A);
%std(A);
kappaList = [ 7,7.5862,8.1724,8.7586,9.3448,9.931,10.5172,11.1034,11.6897,12.2759,...
                   12.8621,13.4483,14.0345,14.6207,15.2069,15.7931,16.3793,16.9655,...
                   17.5517,18.1379,18.7241,19.3103,19.8966,20.4828,21.069,21.6552,...
                   22.2414,22.8276,23.4138,24 ];
sizeA = [150 1];
optVal = zeros(length(kappaList), 2);
constVal = zeros(length(kappaList), 2);
k = 1;
for kappa = kappaList
    inName= strcat('test7_optObjVal_kappa_', num2str(kappa), '.txt');
    fileID = fopen(inName, 'r');
    formatSpec = '%f';
    A = fscanf(fileID, formatSpec, sizeA);
    optVal(k, :) = [nanmean(A), nanstd(A)];
    inName= strcat('test7_optContVal_kappa_', num2str(kappa), '.txt');
    fileID = fopen(inName, 'r');
    A = fscanf(fileID, formatSpec, sizeA);
    formatSpec = '%f';
    constVal(k, :) = [nanmean(A), nanstd(A)];
    k = k + 1;
end

sumtab = horzcat(kappaList', optVal, constVal);
% txt file                     
outName= strcat('val_summary.txt');
dlmwrite(outName, sumtab, '-append');
%     obj_list_mean= mean(obj_list; 2);
%     obj_list_std= std(obj_list; 0; 2);
%     const_list_mean = mean(const_list; 2);
%     const_list_std = std(const_list; 0; 2);
%     theta1_list_mean = mean(theta1_list;  2);
%     theta1_list_std = std(theta1_list; 0; 2);
%     theta2_list_mean = mean(theta2_list; 2);
%     theta2_list_std = std(theta2_list; 0; 2);
%     set_list = repmat(setting; 22; 1);
% tex file
outName2 = strcat('result.tex');
FID = fopen(outName2, 'w');
% fprintf(FID, '\\begin{tabular}{rrrrrrrrrr}\\hline \n');
fprintf(FID, '\\begin{tabular}{rrrrr}\\hline \n');
% fprintf(FID, 'setting & $\\nu$  & $\\wh{V}_1(\\wh{\\bs{\\theta}}_{\\nu})$ & $std(\\wh{V}_1)$ & $\\wh{V}_2(\\wh{\\bs{\\theta}}_{\\nu})$ & $std(\\wh{V}_2)$ & $\\wh{\\theta}_{\\nu,1}$ & $std(\\wh{\\theta}_{\\nu,1})$ & $\\wh{\\theta}_{\\nu,2}$ & $std(\\wh{\\theta}_{\\nu,2})$ \\\\ \\hline \n');
fprintf(FID, '$\\nu$  & $\\wh{V}_1(\\wh{\\bs{\\theta}}_{\\nu})$ & $std(\\wh{V}_1)$ & $\\wh{V}_2(\\wh{\\bs{\\theta}}_{\\nu})$ & $std(\\wh{V}_2)$ \\\\ \\hline \n');
printtab = sumtab;
for k=11:size(printtab,1)
    printline = printtab(k, :);    
    fprintf(FID, '%8.2f & %8.2f & %8.2f  & %8.2f &  %8.2f  \\\\ ', printline);
    if k==size(printtab,1)
        fprintf(FID, '\\hline ');
    end
    fprintf(FID, '\n');
end
fprintf(FID, '\\end{tabular}\n');
fclose(FID);
  
%   % plot
%   width=10;
%   height=8;
%   x0=1;
%   y0=1;
%   figure('Units','inches', 'Position',[x0 y0 width height],...
%   'PaperPositionMode','auto');
% % (x0,y0) = position of the lower left side of the figure
%    c = printtab(:,2);
%    y1 = printtab(:,3);
%    std1 = printtab(:,4);
%    y1U = y1 + 1.96*std1/sqrt(n);
%    y1L = y1 - 1.96*std1/sqrt(n);
%    y2 = printtab(:,5);
%    std2 = printtab(:,6);
%    y2U = y2 + 1.96*std2/sqrt(n);
%    y2L = y2 - 1.96*std2/sqrt(n);
%    plot(c,y1,'--ro',c,y2,'-.bo','LineWidth',1.5);
% %    hold on;
% %    plot(c,y1U,':', c,y1L,':', c,y2U,':',c,y2L,':','LineWidth',1.5);
% %    hold off;
%    hline = refline(1,0);
%    set(hline,'LineStyle',':', 'LineWidth',1.5);
%    xlabel({'$\nu$ bound on secondary potential outcome'}, ...
%              'interpreter' ,'latex', 'FontSize',15 )
%    ylabel({'$\widehat{V}$ values of estimated constrained optimal regimes'},...
%             'interpreter' ,'latex', 'FontSize',15 )
%    title({'Efficient Frontier Plot $\widehat{V}_1$ / $\widehat{V}_2$ vs. $\nu$'},...
%           'interpreter' ,'latex', 'FontSize',15);
%    legend({'$\widehat{V}_1$ vs. $\nu$',  '$\widehat{V}_2$ vs. $\nu$'}, ...
%              'interpreter' ,'latex', 'Location','SouthEast','FontSize',15);
%   % axis([0 t(end) -1.5 1.5]);
%    set(gca, 'Units','normalized', ...
%        'FontUnits','points',... 
%        'FontWeight','normal',... 
%        'FontSize',15);
%    print(strcat('efficient_plot', num2str(setting)), '-dpdf', '-bestfit' ) ;