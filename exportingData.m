%%%%% EXPORT THE DATA %%%%%

%fileName = '2_HD53unli_CIN04_m20_10'
fileName = 'BESCAMsemilim_CIN005_PWGD0011_m10_10';

    % Mass
        fileNameMass = strcat(fileName, '_Mass.txt');
        writematrix(Mass, fileNameMass);
    % rho_average
        fileNameRhoAver = strcat(fileName, '_rhoAverage.txt');
        writematrix( rho_average, fileNameRhoAver);
    % fitness_average
        fileNamefitnessAverage = strcat(fileName, '_fitnessAverage.txt');
        writematrix(fitness_average, fileNamefitnessAverage);
    % net_karyotypes_save
        fileNameNetKarTime = strcat(fileName, '_netKaryotypeTimePoints.txt');
        writematrix(net_karyotypes_save, fileNameNetKarTime);
    % mean_netKar_save
        fileNameNetKar = strcat(fileName, '_netKaryotype.txt');
        writematrix(mean_netKar_save, fileNameNetKar);
    % mode_kar_replicates EXPORTAR A MATLAB TMB
        fileNameModeKarRep = strcat(fileName, '_modeKar.txt');
        writematrix(mode_kar_replicates, fileNameModeKarRep);
    % mode_karyotypes_time
        fileNameModeKar = strcat( fileName, '_modeKarTimepoints.txt');
        writematrix(mode_karyotypes_time, fileNameModeKar);
    % bALL subtypes 
        fileNameBALLsubtypes = strcat( fileName, '_bALLsubtypes.txt');
        writematrix(bALLsubtypes, fileNameBALLsubtypes);

    % Heatmaps 
        fileNameHeatmapInitial = strcat(fileName, "_heatmapInitial.pdf")
        exportgraphics(heatmap1,fileNameHeatmapInitial)

        fileNameHeatmapFinal = strcat(fileName, "_heatmapFinal.pdf")
        exportgraphics(heatmap2,fileNameHeatmapFinal)

        fileNameModeKaryotypes = strcat(fileName, "_modeKaryotypes.pdf")
        exportgraphics(modeKaryotypes,fileNameModeKaryotypes)

