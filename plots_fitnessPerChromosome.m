%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          CIN and B-ALL progression       %
% 
%   Authors
%
%       Carmen Ortega Sabater - PhD Student
%           carmen.ortegasabater@uclm.es
%
%       Víctor M. Pérez García  - PI   victor.perezgarcia@uclm.es             
%       Gabriel Fernández Calvo - PI   gabriel.fernandez@uclm.es           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary file 4. This script allow us to generate plots.
% We also study the change in average rho with time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              PLOTS              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters summary
parametersSummary = [replicates dt simSteps rho_init probCIN ProbDeath s] 
simulationSteps = (1:200)';

% preprocess the data
% replace zeros in matrix with NaN
save_cells_karyotypes = zeros( size(cells_karyotype,1), size(cells_karyotype,2), size(saveKaryotypes,2), size(cells_karyotype,4));
for i = 1:replicates
    save_cells_karyotypes(:,:,:,i) = cells_karyotype(:,:,saveKaryotypes(i,:),i);
end 
save_cells_karyotypes(save_cells_karyotypes == 0) = NaN; 
hist_fitness(hist_fitness == 0) = NaN; 
rho_average(rho_average == 0) = NaN;
fitness_average(fitness_average == 0) = NaN; 
Mass(Mass==0) = NaN; 
net_karyotypes(net_karyotypes == 0 ) = NaN; 

    % find the mean karyotypes for replicates
    net_karyotypes_replicates = mean(net_karyotypes, 3, 'omitnan');
        % extract some timepoints
        savekar = saveKaryotypes(1,:);
        net_karyotypes_save = net_karyotypes_replicates(:,savekar);
    mean_netKar = mean(net_karyotypes_replicates, 1 ,'omitnan')';
        mean_netKar_presave = mean(net_karyotypes, 1, 'omitnan');
        mean_netKar_save = zeros(up_simSteps, replicates);
            for i = 1:replicates
                 mean_netKar_save(:,i) = mean_netKar_presave(:,:,i)
            end 

            figure()
tiledlayout(1,3)

% 1. Population growth plot
    % Draw all graphs in the same figure
    nexttile
    hold on
    box on
    ax = gca;
    ax.FontSize = 14; 
    ax.LineWidth = 3;
    set(gca,'FontName','Helvetica Neue' );
    %plot_totpop = plot( fit_x, fit_y, '--', 'LineWidth', 1.5, 'Color', [171 174 171]/255 ); %grey color
    for i = 1:replicates
        plot( Mass(:,i), 'Color', [182/255 207/255 208/255] );
    end
    plot( mean(Mass,2), 'LineWidth', 2, 'Color', 'r' );
    xlabel( '$\mathrm{Time \; [days]}$', 'FontSize', 18, 'Interpreter','latex' ); 
    ylabel('$\log N(t)$', 'FontName','Arial','FontSize', 18, 'Interpreter','latex');
    xlim([0 150]);
    %saveas(gcf, 'case10_totpop_6years_initNumber1e9_pdeath2_10e-1.png')
    hold off

% 2. Average rho in time
    nexttile
    hold on
    box on
    ax = gca;
    ax.FontSize = 14; 
    ax.LineWidth = 3;
    set(gca,'FontName','Helvetica Neue' );
    %plot_totpop = plot( fit_x, fit_y, '--', 'LineWidth', 1.5, 'Color', [171 174 171]/255 ); %grey color
    for i = 1:replicates
        plot( rho_average(:,i), 'Color', [182/255 207/255 208/255] );
    end
    plot( mean(rho_average,2), 'LineWidth', 2, 'Color', 'r' );
    xlabel( '$\mathrm{Time \; [days]}$', 'FontSize', 18, 'Interpreter','latex' ); 
    ylabel('$\bar{\rho}$','FontName','Arial','FontSize', 18, 'Interpreter','latex');
    xlim([0 150]);

    %saveas(gcf, 'case10_totpop_6years_initNumber1e9_pdeath2_10e-1.png')
    hold off

% 3. Average fitness in time
    nexttile
    hold on
    box on
    ax = gca;
    ax.FontSize = 14; 
    ax.LineWidth = 3;
    set(gca,'FontName','Helvetica Neue' );
    %plot_totpop = plot( fit_x, fit_y, '--', 'LineWidth', 1.5, 'Color', [171 174 171]/255 ); %grey color
    for i = 1:replicates
        plot( fitness_average(:,i), 'Color', [182/255 207/255 208/255] );
    end
    plot(mean(fitness_average,2), 'LineWidth', 2, 'Color', 'r' );
    xlabel( '$\mathrm{Time \; [days]}$', 'FontSize', 18, 'Interpreter','latex' ); 
    ylabel('$\mathrm{Fitness}$', 'FontName','Arial', 'FontSize', 18, 'Interpreter','latex');
    xlim([0 150]);

    %saveas(gcf, 'case10_totpop_6years_initNumber1e9_pdeath2_10e-1.png')
    hold off

  % Export data
    mean_fitness = mean( fitness_average,2);
    std_fitness = std(fitness_average, [], 2, 'omitnan');

%     % 1b. Shaded error bars
%     Mass(Mass == NaN ) = 500;
%     logMass = log(Mass);
%     std_Mass = std(logMass, [], 2, 'omitnan');
%     meanLogMass = mean(logMass,2); 
%     logMassSDplus = meanLogMass+std_Mass;
%     logMassSDminus = meanLogMass-std_Mass;
%     logMassSD = [logMassSDplus; logMassSDminus];
% 
%     % Draw all graphs in the same figure
%     figure
%     hold on
%     box on
%     ax = gca;
%     ax.FontSize = 14; 
%     ax.LineWidth = 3;
%     set(gca,'FontName','Helvetica Neue' );
%     fill([simulationSteps; flipud(simulationSteps)], [meanLogMass-std_Mass; flipud(meanLogMass+std_Mass)], [.9,.9,.9], 'LineStyle','none')
%     line(simulationSteps, meanLogMass);
%     %plot( log(mean(Mass,2)), 'LineWidth', 2, 'Color', 'r' );
%     xlabel( '$\mathrm{Time \; [days]}$', 'FontSize', 18, 'Interpreter','latex' ); 
%     ylabel('$\log N(t)$', 'FontName','Arial','FontSize', 18, 'Interpreter','latex');
%     xlim([0 200]);
%     %saveas(gcf, 'case10_totpop_6years_initNumber1e9_pdeath2_10e-1.png')
%     hold off


% 4. Evolution of net karyotype in time 
 

    figure
    hold on
    box on
    ax = gca;
    ax.FontSize = 14; 
    ax.LineWidth = 3;
     for i = 1:replicates
         data = net_karyotypes(:,:,i);
        plot( mean(data,1,'omitnan'), 'Color', [182/255 207/255 208/255] );
    end
    plot( mean(net_karyotypes_replicates, 1 ,'omitnan'), 'LineWidth', 2, 'Color', 'r' );
    xlabel( 'Time (days)', 'FontSize', 18 ); 
    ylabel('Karyotype','FontName','Arial','FontSize', 18);
    %xlim([0 150]);
    %ylim([91 94]);

    %saveas(gcf, 'case10_totpop_6years_initNumber1e9_pdeath2_10e-1.png')
    hold off

mean_netKar = mean(net_karyotypes_replicates, 1 ,'omitnan')';    

% % 4. Evolution of fitness distribution in time     
%     % Figure 4. Evolution of phenotypic distribution with time. 
% 
% % [f1,x1] = ksdensity(1/hist_fitness(:,1)); 
% % [f2,x2] = ksdensity(1/hist_fitness(:,floor(s/2))); 
% % [f3, x3] = ksdensity(1/hist_fitness(:,s))
% 
% %timepoint 1
% numIntervals = 20;
% intervalWidth = (max(1/hist_fitness(:)) - min(1/hist_fitness(:)))/numIntervals;
% x = 1:intervalWidth:6;
% 
% ncount_h1 = histc(1./hist_fitness(:,1),x);
% ncount_h2 = histc(1./hist_fitness(:,floor(s/2)), x);
% ncount_hend = histc(1./hist_fitness(:,s),x);
% 
% relativeFreq_h1 = ncount_h1/length(hist_fitness(:,1));
% relativeFreq_h2 = ncount_h2/length(hist_fitness(:,floor(s/2)));
% relativeFreq_hend = ncount_hend/length(hist_fitness(:,s));
% 
% 
% figure();
% tiledlayout(3,1)
% 
% nexttile
% box on 
% bar(x-intervalWidth/2, relativeFreq_hend,1)
% xlim([min(x) max(x)])
% %set(gca, 'xtick', x)
% 
% nexttile 
% bar(x-intervalWidth/2, relativeFreq_h2,1)
% xlim([min(x) max(x)])
% 
% nexttile
% bar(x-intervalWidth/2, relativeFreq_h1,1)
% xlim([min(x) max(x)])
% 
% ax = gca;
% ax.FontSize = 26; 
% h1 = histogram(hist_fitness(:,1)/nnz(~isnan(hist_fitness(:,1))), 20)
% hold on
% hmid = histogram(hist_fitness(:,floor(s/2))/nnz(~isnan(hist_fitness(:,floor(s/2)))), 20)
% hend = histogram(hist_fitness(:,s)/nnz(~isnan(hist_fitness(:,s))), 20)
% hold off
% 
% 
% a1 = area(hist_fitness, pF_1, 'Linewidth', 2, 'FaceAlpha', 0.3, 'FaceColor',[0 0.25 0.25]);
% %a2 = area(rho_day, pF_midtime, 'Linewidth', 2, 'FaceAlpha', 0.7, 'FaceColor', [187/255 161/255 79/255]);
% %a3 = area(rho_day, pF_end, 'Linewidth', 2, 'FaceAlpha', 0.6, 'FaceColor',[155/255 20/255 14/255]);
% 
% a1.EdgeColor = [0 0.25 0.25];
% a2.EdgeColor = [187/255 161/255 79/255];
% a3.EdgeColor = [155/255 20/255 14/255];
% ylim([0 0.40])
% xlabel('$\mathrm{Proliferation \; rate} \; (\mathrm{days}^{-1})$', 'FontSize', 22, 'Interpreter','latex' );
% ylabel('$\mathrm{Phenotypic \; frequency}$', 'FontSize', 22, 'Interpreter','latex');
% legend({'$t = 0 \; \mathrm{days}$', '$t = 15 \; \mathrm{days}$', '$t = 30 \; \mathrm{days}$' },'Interpreter','latex')
% legend boxoff  
% hold off


% 5. Multiplot chromosomes/karyotype distribution. Evolution of MEAN karyotype
% in time
    % find mean karyotypes for t=x1, t=x2, t=x3?
    mean_kar_replicates = mean( save_cells_karyotypes, 4, 'omitnan');
    mean_karyotypes_time = zeros(size(save_cells_karyotypes,1), size(save_cells_karyotypes,3) );
    for i = 1:size(mean_kar_replicates, 3)
        temp_kar = mean_kar_replicates(:,:,i);
        mean_karyotypes_time(:,i) = round(mean(temp_kar, 2, 'omitnan'));
    end 
    
    grayColor = [.7 .7 .7];

    fig = figure;
    hold on
    box on 
%     ax = gca;
%     ax.FontSize = 14; 
%     ax.Layer = "top"; 
%     ax.LineWidth = 3;
        for i = 1:size(mean_karyotypes_time, 2)
            subplot(size(mean_karyotypes_time,2), 1, i)
            b = bar(mean_karyotypes_time(:,i), 'BarWidth', 1, 'FaceColor', grayColor)
            b.FaceAlpha = 0.6;
            b.EdgeColor = 'none';
%             xlabel( 'Chr ID', 'FontSize', 16 )
%             ylabel( 'Copies', 'FontSize', 18 );
            title(['t=',num2str(saveKaryotypes(i)),' days'], 'FontSize', 18 );
            yline(2,'r--', LineWidth = 2 )
        end
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.YLim = [1 4];
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'Copies');
        xlabel(han, 'Chr ID');
        title( han, 'Mean Karyotype - CIN = 0.1' )
        set(findobj(gcf,'type','axes'), 'ylim', [1 max(mean_karyotypes_time(:))+1], 'FontName','Arial','FontSize',12, 'LineWidth', 1.5);
        

% 6. Evolution of MODE karyotype in time
    % find mode karyotypes for t=x1, t=x2, t=x3?
    grayColor = [.7 .7 .7];
    mode_kar_replicates = mean( save_cells_karyotypes, 4);
    mode_karyotypes_time = zeros(size(mode_kar_replicates,1), size(mode_kar_replicates,3) );
    time = floor(mean(saveKaryotypes,1)); %find average time vector
  
    for i = 1:size(save_cells_karyotypes, 3)
        mode_karyotypes_time(:,i) = mode(save_cells_karyotypes(:,:,i), 2);
    end 
    
    fig = figure;
    hold on
    box on 
%     ax = gca;
%     ax.FontSize = 14; 
%     ax.Layer = "top"; 
%     ax.LineWidth = 3;
        for i = 1:size(mode_karyotypes_time, 2)
            subplot(size(mode_karyotypes_time,2), 1, i)
            b = bar(mode_karyotypes_time(:,i), 'BarWidth', 1, 'FaceColor', grayColor)
            b.FaceAlpha = 0.6;
            b.EdgeColor = 'none';
%             xlabel( 'Chr ID', '', 16 )
%         FontSize    ylabel( 'Copies', 'FontSize', 18 );
            title(['t=',num2str(time(i)),' days'], 'FontSize', 18 );
            yline(2,'r--', LineWidth= 2 )
        end

        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        %han.YLim = [0 6];
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'Copies');
        xlabel(han, 'Chr ID');
        title( han, 'Mode Karyotype ' )
        set(findobj(gcf,'type','axes'), 'ylim', [0 6], 'FontName','Arial','FontSize',12, 'LineWidth', 1.5);
        modeKaryotypes = gcf;
        

% 7. heatmaps 
% initial time point 
% final time point

figure
box on
ax = gca 
ax.LineWidth = 6
colormap('turbo')
imagesc(mode_kar_replicates(:,1:500,1)');
colorbar
lims = clim 
clim([1 6])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.5)
heatmap1 = gcf; 

% final time point
figure
box on
ax = gca 
ax.LineWidth = 6
colormap('turbo')
imagesc(mode_kar_replicates(:,:,6)')
colorbar
lims = clim 
clim([1 6])
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',16, 'LineWidth', 1.5)
heatmap2 = gcf; 


% 7. Cumulative barplot - B-ALL subtypes 
% x axis - time
% y axis - % cells with a certain B-ALL subtype
% near-haploid: 24-30 chr
% low-hypo: 31-39 chr
% high-hyper: 51-67 chr
% diploid: 44-46 chr

% seed = net_karyotypes
% timepoints = save_karyotypes
net_karyotypes_timePoints = zeros( 1200000, size(saveKaryotypes,2),replicates);
for i = 1:size(net_karyotypes, 3)
    net_karyotypes_timePoints(:,:,i) = net_karyotypes(:,saveKaryotypes(i,:), i); 
end
  
bALLsubtypes = zeros( 5, 6, replicates); %subtypes; time points; replicates





cellcount = zeros(replicates, 6);
for k = 1:replicates
    for i = 1:size(net_karyotypes_timePoints,2)
        cellcount(k,i) = nnz(~isnan(net_karyotypes_timePoints(:,i,k)));
    end
end

for k = 1:replicates
    for i = 1:size(net_karyotypes_timePoints,2)
        for j = 1:5 %b-all subtypes
            %select the data
            data = net_karyotypes_timePoints(:,i,k);
            data = data(~isnan(data));
            bALLsubtypes(1,i,k) = length(data(data>=24 & data<=30))%/length(data); %nhaploid
            bALLsubtypes(2,i,k) = length(data(data>=31 & data<=39))%/length(data); %Low-Hypodiploid
            bALLsubtypes(3,i,k) = length(data(data>=45 & data<=47))%/length(data); % euploid
            bALLsubtypes(4,i,k) = length(data(data>=51 & data<=67))%/length(data); % high-hyperdiploid
            bALLsubtypes(5,i,k) = (length(data)-sum(bALLsubtypes(1:4,i,k)))%/length(data); %other
            bALLsubtypes(:,i,k) = bALLsubtypes(:,i,k)./length(data);
        end
    end
end

    bALLsubtypes_replicates = mean(bALLsubtypes,3);
    time = floor(mean(saveKaryotypes,1)); %find average time vector
    % bALLsubtypes in case we want to store all replicates
    bALLsubtypes_pre = num2cell(bALLsubtypes_replicates, [1 2]);
    bALLsubtypesConcatenated = vertcat(bALLsubtypes_pre{:});    

figure
box on
hold on
bar(time, bALLsubtypes_replicates.', 'stacked')
legend('nHaploid', 'Low Hypodiploid', 'Euploid', 'High-hyperdiploid', 'Other');
ylabel('Rel Frequency', 'FontSize',24);
xlabel( 'Time (days)', 'FontSize',24);
title( 'Evolution of karyotypes in time - B-ALL subtypes' )
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',18, 'LineWidth', 1.5);


% % Prepare data for MULLER PLOT
% bALLsubtypes = zeros( 5, s, replicates); %subtypes; time points; replicates
% 
% 
% for k = 1:replicates
%     for i = 1:size(net_karyotypes,2)
%         for j = 1:5 %b-all subtypes
%             %select the data
%             data = net_karyotypes(:,i,k);
%             data = data(~isnan(data));
%             bALLsubtypes_muller(1,i,k) = length(data(data>=24 & data<=30))%/length(data); %nhaploid
%             bALLsubtypes_muller(2,i,k) = length(data(data>=31 & data<=39))%/length(data); %Low-Hypodiploid
%             bALLsubtypes_muller(3,i,k) = length(data(data>=45 & data<=47))%/length(data); % euploid
%             bALLsubtypes_muller(4,i,k) = length(data(data>=51 & data<=67))%/length(data); % high-hyperdiploid
%             bALLsubtypes_muller(5,i,k) = (length(data)-sum(bALLsubtypes(1:4,i,k)))%/length(data); %other
%             bALLsubtypes_muller(:,i,k) = bALLsubtypes(:,i,k)./length(data);
%         end
%     end
% end



