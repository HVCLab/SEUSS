%% plot GA for each condition, avg across top 10 chan
addpath 'C:\Users\Camille\OneDrive\CNSP_workshop\fieldtrip-20231025\fieldtrip-20231025'
ft_defaults

% archivos = dir(fullfile(processed_datapath,'preprocessed','clusterinput','int','*avg*'))
GAint5 = load('grandAverages_iepc_5cond_100424_filtwidth0.5.mat');
GAint6 = load('grandAverages_iepc_6cond_100724_filtwidth0.5.mat');
load channel_ind_mtr.mat
load frequencies.mat
load EGI_layout129.lay.mat

%%

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';

GAfig = {};
GAfig{1} = GAint6.GA{1};
GAfig{2} = GAint6.GA{2};
GAfig{3} = ft_math(cfg,GAint6.GA{1},GAint6.GA{2});
GAfig{4} = GAint6.GA{3};
GAfig{5} = GAint6.GA{4};
GAfig{6} = ft_math(cfg,GAint6.GA{3},GAint6.GA{4});
GAfig{7} = ft_math(cfg,GAint6.GA{1},GAint6.GA{3});

GAfig{8} = ft_math(cfg,GAint6.GA{2},GAint6.GA{4});
GAfig{9} = ft_math(cfg,GAfig{7},GAfig{8});


figure(1)
titulos = {'UI','UR','irr vs reg UNS',...
    'SI','SR','irr vs reg STR',...
    'uns vs str IRR','uns vs str REG','int'};

% Main effects
cfg = [];
cfg.baseline = 'no';
% cfg.zlim         = [-0.02 0.02];
cfg.channel      = chanindmtr.ind(1:10);
cfg.layout       = EGI_layout129;
cfg.parameter = 'powspctrm';

j = 1
for i = 1:9
    cfg.figure = subplot(3,3,i);
    if i == 1 || i == 2 || i == 4 || i == 5
     cfg.zlim         = [0.03 0.09];
    else
     cfg.zlim         = [-0.02 0.025];
    end
    
    ft_singleplotTFR(cfg, GAfig{i}); 
    title(titulos{i});
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    hold on
%     plot(0,2,'*')
     if  i > 6
     contour(time,frq,cluster2d{j},3,'color','w','LineWidth',2)
     j = j+1
     end
%     
%     
%     ax = gca(); % this Gets the Current Axis so we can set properties
% %     ax.XAxisLocation = 'origin';
%     ax.YAxisLocation = 'origin';
%     ax.TickDir = 'out';
%     
%     set(gca, 'YColor', 'white');

    % Remove the box around the plot, while we're at it:
%     box off;

    % And move the x-axis label to underneath the axis:
    % ax.XLabel.Position(2) = -30;
end
