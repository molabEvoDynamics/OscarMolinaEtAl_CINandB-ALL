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

% Auxiliary file 2. This file contains the core of the model to perform
% the simulations. 

% Initialize counters and memory
clear all; close all;
tic;

% Load parameters and data structures
run( 'parameters.m' );
run( "ngs_initialCondition.m"); % this will load the karyotypes from real patients
    
for m = 1:replicates
    %for c = 1:size(probCIN_mat,2)        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %          Initial condition          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear cells matrix 
        cells = [];


        
 %-------------------PROOF OF CONCEPT ------------------------------------%       
% Starting with an aneuploid cell population sampled from a normal distribution
% 
%         initialKaryotypes = zeros( totalCells, chromosomes); 
%         
%         % Initial cells are distributed across phenotypes with a
%         % Gaussian of width InitialWidth
%         InitialWidth = 1;
%         
%         % fill initial karyotypes
%         for i = 1:size(initialKaryotypes, 1)
%         initialKaryotypes(i,:) = abs(round(basicNumberOfCopies + randn(chromosomes,InitialWidth)));
%         end 
%         
%         % remove cells which show chromosomes with 0 copies
%         initialKaryotypes( ~any(initialKaryotypes,2), : ) = [];  
%         initialOut =  initialKaryotypes(all( initialKaryotypes,2),:);
%         initialOut = initialOut';
%         
%         initialCellNumber = 500; %cells
%         %pick N=initialCellNumber cells to start the simulations from 
%         % initialOut 
%         k = randperm(size(initialOut,2),initialCellNumber);
%         initKaryotypes = initialOut(:,k);
%-----------------------PROOF OF CONCEPT 4N----------------------------%         

%         % STARTING WITH A 4N POPULATION
%         % set initial cell population
%         initialCellNumber = 50; %cells
%         initKaryotypes = ones(chromosomes, initialCellNumber)*4;

%--------------------DATA FROM PATIENTS (NGS) -------------------------%
% % CHOOSE PATIENT
% Modal karyotypes, considering the standard deviation from patients'
% karyotype
% mean_pat_kar = hyperD_hd53_mean; 
% sd_pat_kar = hyperD_hd53_sd; 

% mat_p = round(mean_pat_kar.*randn(totalCells,1)+sd_pat_kar); 
%  mat_p( ~any(mat_p,2), : ) = [];  
%         initialOut_p =  mat_p(all( mat_p,2),:);
%         allPositiveRows = all(initialOut_p>0, 2);
%         initialOut_p = initialOut_p(allPositiveRows, :); 
%         initialOut_p = initialOut_p';
%        Fixed initial karyotypes from patients (clonal population)
         % Hyperdiploid 
         hd53 = [2 2 2 3 2 3 3 3 2 3 2 2 2 4 2 2 3 3 2 2 4 3 3];
         initKaryotypes = repmat(hd53', 1, initialCellNumber);
%        % Near haploid (not endoreduplicated clone)
%          bescam = [1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 1 1 2 1 2];
%          initKaryotypes = repmat(bescam', 1, initialCellNumber);
 
        % aneuploid initial condition
         % aneuploidExample = [2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 3 2 2 2 2 2];
         % initKaryotypes = repmat(aneuploidExample', 1, initialCellNumber);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %          START SIMULATIONS          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        cells(:, 1:initialCellNumber) = initKaryotypes;
        % save initial population in cells_karyotype at s=1
        cells_karyotype(:, 1:initialCellNumber, 1, m) = initKaryotypes; 

        % fill initial net karyotypes
        net_karyotypes(1:size(cells,2), 1, m) = sum(initKaryotypes,1); 


        % find the fitness for the initial karyotypes
       
        fitness = abs(sum((initKaryotypes).*fitness_hybrid,1));
        hist_fitness(1:initialCellNumber, 1) = fitness; 

        rho_average(1,m) = rho_init; 

        for s=2:up_simSteps % Number of time steps   
           
            % 1. Proliferation stage
%             % find the fitness for each karyotype
                 fitness = sum(cells.*fitness_hybrid,1);
%             % 1a. Fitness threshold;
%                 % set a fitness threshold which determines if cells divide
%                 % or not
%                 fitnessThreshold = 2; 
%                 % get newborn indexes 
%                 dividing = find(fitness < fitnessThreshold);           

            % 1c. Matriz de diferencias
             N = length(fitness); % find total population

%              if s == 2
%                 rho_average(s,m) = rho_init;
%              end 


%              euploid = 2*ones(size(cells));
%              cells_f = cells - euploid;
%              cells_f2 = abs(sum(cells_f.*fitness_hybrid,1));
             % approach b 
             
            net_karyotypes_s = sum(cells,1);         % find net karyotypes
            cells_f2 = fitness./net_karyotypes_s;    % adjust F per netkar
          
            rho = rho_average(s-1,m)*(1/dt)*(1-N/K); % cell division
            newborn = binormal(N, rho);
% 
            % SELECTION TOWARDS ANEUPLOIDY - Select based on max fitness
            [maxFitness, dividing] = maxk(cells_f2, newborn);
            
%             % NEUTRAL SELECTION - Random selection once we select cells 
%             dividing = randperm(length(cells_f2), newborn);

%             %MIXED MODEL
%             % Select newborns between SIGMA % with max FITNESS
%             [maxFitness30, dividingCandidates] = maxk(cells_f2, round(sigma*length(cells_f2)));
%             dividingCandidatesLength = length(dividingCandidates);
%             dividingPositions = randperm(dividingCandidatesLength, newborn);
%             dividing = dividingCandidates(dividingPositions);

%                             %----Whole Genome Duplication (WGD)----%
%                 Implementing whole genome duplications (WGD)
                 N_WGD = 0; 
% %                       %find if any dividing cells have less than 40 chr count
%                        candidatesWGD = sum(net_karyotypes_s(1:length(fitness)) < 40);
%                        if candidatesWGD > 0
%                         %Identify who they are (find positions in
%                          %Cells) 
%                             posCandidatesWGD = find( net_karyotypes_s(1:length(fitness)) < 40 & net_karyotypes_s(1:length(fitness)) > 0);                         
%                             %Is WGD linked to cell division? YES
%                             intersectionDividingWGD = intersect(dividing, posCandidatesWGD);
%                          %Identify how many will suffer WGD
%                             N_WGD = binormal(length(intersectionDividingWGD), probWGD);
%                         %within candidatesWGD, find which of them will
%                         %suffer WGD
%                             id_posWGD = randperm(length(intersectionDividingWGD), N_WGD);
%                             posWGD = intersectionDividingWGD(id_posWGD); 
%                             cells(:, posWGD) = cells(:, posWGD).*2; 
%                        end 
                        %----WGD not linked to proliferation----%
                            %find if any dividing cells have less than 40 chr count
                            candidatesWGD = sum(net_karyotypes_s(1:length(fitness)) < 40);
                            if candidatesWGD > 0
                                %Identify who they are (find positions in
                                %Cells)
                                posCandidatesWGD = find( net_karyotypes_s(1:length(fitness)) < 40 & net_karyotypes_s(1:length(fitness)) > 0);
                                %Is WGD linked to cell division? YES
                                %Identify how many will suffer WGD
                                N_WGD = binormal(length(posCandidatesWGD), probWGD);
                                %within candidatesWGD, find which of them will
                                %suffer WGD
                                id_posWGD = randperm(length(posCandidatesWGD), N_WGD);
                                posWGD = posCandidatesWGD(id_posWGD);
                                cells(:, posWGD) = cells(:, posWGD).*2;
                            end
                       %-----------WGD not linked to proliferation-------%


%                 %----WGD----%    

            % CIN - Cells undergoing chromosomal instability
            
            dividingN = length(dividing);
            cinCells = binormal(dividingN, probCIN); %how many cells have an altered karyotype

            if cinCells ~= 0

           changingCells = randperm(dividingN,cinCells); % which cells of the newborn will undergoCIN
% % % %                 %------ LIMITING LAGGING CHROMOSOMES (1-4 MAX)-------%
%                 % initiate the empty matrix
%                 extracopies2 = zeros(23, cinCells);
%                 whichLaggingChr = zeros(23,1);
% 
%                 for p = 1:cinCells
%                     % Iterate over cinCells
%                     % Choose a random number (K) between 1 and 3. Setting 3 as the
%                     % maximum possible number of lagging chromosomes 
%                     nLaggingChr = randi([1 4]);
%     
%                     % Choose K numbers between 1 and 23 to see which
%                     % chromosomes are acquiring extra copies
%                     whichLaggingChr = fix(1+(23-1).*rand(nLaggingChr,1));
% 
%                     % Fill the matrix in the corresponding positions 
%                     extracopies2(whichLaggingChr, 1) = 1;
%                 end 
% 
%                 % Update daughter2
%                 daughter2 = cells(:,dividing); % we extract all the cells that are dividing
%                 % then we can update daughter1 which underwent CIN 
%                 % in matrix of karyotypes
%                 daughter1 = dividing(changingCells);
%                 cells(:,daughter1) = cells(:,daughter1)+extracopies2;
% 
%                 % We need to substract the extracopies from one of the
%                 % daughter cells
%                 daughter2(:,changingCells) = daughter2(:,changingCells) - extracopies2; 
% % %                 %------ LIMITING LAGGING CHROMOSOMES (1-4) -------%

% %                 %------ LIMITING LAGGING CHROMOSOMES (1 MAX)-------%
%                 % initiate the empty matrix
%                 extracopies2 = zeros(23, cinCells);
%                 for p = 1:cinCells
%                     % Iterate over cinCells
%     
%                     nLaggingChr = 1
%     
%                     % Choose K numbers between 1 and 23 to see which
%                     % chromosomes are acquiring extra copies
%                     whichLaggingChr = randi([1 23],1,nLaggingChr);
%     
%                     % Fill the matrix in the corresponding positions 
%                     extracopies2(whichLaggingChr, p) = 1;
%                 end 
% 
%                 % Update daughter2
%                 daughter2 = cells(:,dividing); % we extract all the cells that are dividing
%                 % then we can update daughter1 which underwent CIN 
%                 % in matrix of karyotypes
%                 daughter1 = dividing(changingCells);
%                 cells(:,daughter1) = cells(:,daughter1)+extracopies2;
% 
%                 % We need to substract the extracopies from one of the
%                 % daughter cells
%                 daughter2(:,changingCells) = daughter2(:,changingCells) - extracopies2; 
% %                 %------ LIMITING LAGGING CHROMOSOMES -------%



% 
%                 %------ UNLIMITED LAGGING CHROMOSOMES ------%
                extracopies = randi([0 1], 23, cinCells); % matrix of extracopies (gains)
                copyLosses = randperm(numel(extracopies), floor(numel(extracopies)/3)); %losses
                extracopies(copyLosses) = -1;


                % if a cell missegregates, its sister cell must lose or
                % gain the counterpart chromosome copies
                daughter2 = cells(:,dividing); % we extract all the cells that are dividing
                % then we can update daughter1 which underwent CIN 
                % in matrix of karyotypes
                daughter1 = dividing(changingCells);
                cells(:,daughter1) = cells(:,daughter1)+extracopies;

                % We need to substract the extracopies from one of the
                % daughter cells
                daughter2(:,changingCells) = daughter2(:,changingCells) - extracopies;
%                 %------ UNLIMITED LAGGING CHROMOSOMES ------%   

                % update population matrix. We incorporate daugther cells
                % matrices
                cells = [cells daughter2]; 

            else % what happens if any cells suffer CIN
                nonChangingCells = cells(:,dividing);
                cells = [cells nonChangingCells];
            
            end      
            
            % kill cells with >6 copies of any chromosomes
             indexes_1 = find(cells < 1);
             indexes_2 = find(cells >=6 );
% 
             [row1 col1] = ind2sub(size(cells), indexes_1);
             [row2 col2] = ind2sub(size(cells), indexes_2);
             aberrant = [col1; col2];
             aberrant = unique(aberrant);
             cells(:,aberrant) = [];                                            
        
% %        % Dead cells
           % normal death       
             ProbDeath = (N/K)*ProbDeath;
             dead = binormal( size(cells,2), ProbDeath );
%              if dead < 0 
%                  dead = binormal(size(cells,2), ProbDeath);
%              end 

             cells(:,size(cells,2)-dead:size(cells,2))=0;

%            
            % find the fitness for each karyotype
            fitness = sum(cells.*fitness_hybrid,1);
        
%           % Save results
            Mass(s, m) = size(cells, 2);     % Find total population
%                 % find the fitness for each karyotype
             fitness = sum(cells.*fitness_hybrid,1);    
                % find the fitness for each karyotype normalized by total
                % karyotype
            %fitness = sum(cells.*fitness_hybrid,1)./net_karyotypes_s; 
        disp(s)
            % Save only    
            %         if ismember(s, saveKaryotypes) == 1
            %             for k = 1:length(saveKaryotypes)
            %                 position_k = find(saveKaryotypes == k) 
            %                 cells_karyotype(:,1:size(cells,2), position_k, m) = cells;
            %                 hist_fitness(1:size(cells,2),s,position_k, m) = fitness'; 
            %             end 
            %         end
           
           hist_fitness(1:size(cells,2),s,m) = fitness';
           if size(cells,2) >= 10000
            saveCellsSubset = floor(1 + (size(cells,2)-1) .* rand(10000,1));   
            cells_karyotype(:,:, s,  m) = cells(:,saveCellsSubset);
           end 
           
          
           fitness_average_norm(s,m) = mean(cells_f2, 'omitnan');
           fitness_average(s,m) = mean(fitness); 
           net_karyotypes_s = sum(cells,1);   
           net_karyotypes(1:size(cells,2), s, m) = net_karyotypes_s;

%        % update rho_basal depending on fitness

%           if s > 2
            rho_average(s,m) = rho_average(s-1, m)*(1+(fitness_average_norm(s,m)-fitness_average_norm(s-1,m)));
%           end
            nTrack_wgd(s, m) = N_WGD;

      end
    saveKaryotypes(m,:) = floor([1:s/5:s s]);
end        

run( 'plots_fitnessPerChromosome.m')
 
toc;





