function cz_plot_cbp_model5(inputFolder,filename)

% needs statC_all, allDataC,EGI_layout, GAc
addpath('\\cin-storage\oganian_data\code\matlab_code\fieldtrip-master')
addpath('..\util\')
ft_defaults

%     inputFolder = 'results124_080524_nfold5_filt0.1to40_01sOmE_newData_win-1to600_generic';
%     filename = 'cbp_uni_10mtrEl_080524';
    load(fullfile(inputFolder,filename));
    load('EGI_layout129.lay.mat')
    chan = statC_all{1}{1}.cfg.channel;
%      ley = cellfun(@(x) x{1}.condicion, allDataC,'UniformOutput',false);    
    clustalfa = [0.05];     
    pval = 1; % do you want to plot p values?
    nclust = 2; % how many clusters do you want?
    scale = 6; % scale cluster number
    ch = find(ismember(EGI_layout129.label,statC_all{1}{1}.cfg.channel)); % average over channels for trf

j = 1;
for cond = [2,4,6,8]
    
        subplot(2,2,j)
        data  = GAc{cond};
        std_error = (sqrt(mean(data.var(ch,:))) / sqrt(26));
        y = mean(data.avg(ch,:));
        time = data.time;
        
        % Calculate the upper and lower bounds for the shaded area
        upper_bound = y + std_error;
        lower_bound = y - std_error;
    
        hold on;
        plot(time, y);   % TRFs grand averages data cond 1, filtered
        fill([time fliplr(time)], [upper_bound fliplr(lower_bound)], 'k',...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');

        % Add labels and title
        xlabel('Time');
        ylabel('weight');
        if nclust < 3
            ylim([-0.6 0.8])
        end
         % Add legend only for the y variable
%           legend(ley{cond}, 'Location', 'southoutside', 'Orientation', 'horizontal')
         grid on
clear data
        data  = GAc{cond+1};
        std_error = (sqrt(mean(data.var(ch,:))) / sqrt(26));
        y = mean(data.avg(ch,:));
        time = data.time;
        
        % Calculate the upper and lower bounds for the shaded area
        upper_bound = y + std_error;
        lower_bound = y - std_error;
    
        hold on;grid on;
        plot(time, y);   % TRFs grand averages data cond 1, filtered
        fill([time fliplr(time)], [upper_bound fliplr(lower_bound)], 'k',...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');

        % Add labels and title
        xlabel('Time');
        ylabel('weight');
        if nclust < 3
            ylim([-0.6 0.8])
        end
%           legend(ley{cond}, 'Location', 'southoutside', 'Orientation', 'horizontal')
        
       

    stati = statC_all{1}{j}; %change 1 to 2 if you have more than one alpha value in the cell
    if isfield(stati,'posclusters')
        disp('hay clusters')
    pc = stati.posclusterslabelmat.*ismember(stati.posclusterslabelmat,1:nclust);
    else 
        pc = 0;
    end
    if isfield(stati,'negclusters')
    nc = stati.negclusterslabelmat.*ismember(stati.negclusterslabelmat,1:nclust)*(-1);
    else
        nc = 0;
    end
    
    clust = pc + nc;
    
    non_zero_indices = clust ~= 0;
    plot(time(non_zero_indices), clust(non_zero_indices) / scale, '*');

    if pval == 1
        % add p values < alfa
        alfa = 0.025;
        signif_pvalues_i = stati.prob < alfa;

        % Mask for non-zero and significant clusters
        signif_non_zero_indices = non_zero_indices & signif_pvalues_i;
        plot(time(signif_non_zero_indices), clust(signif_non_zero_indices) / scale, '*');
    end
    
    title(statC_all{1}{j}.comparison)
    grid on
    
    j = j+1
    %     % Save figure as PNG
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10*1.5 6*1.5])
%     if pval == 1
%         filename_png = ['pval_TRFwClust_clustalfa' num2str(alfavalue) '-' num2str(clustalfa(alfavalue)) '.png'];
%     else 
%         filename_png = ['TRFwClust_clustalfa' num2str(alfavalue) '-' num2str(clustalfa(alfavalue)) '.png'];
%     end
%     print(filename_png,'-dpng','-r0')
end
hold off