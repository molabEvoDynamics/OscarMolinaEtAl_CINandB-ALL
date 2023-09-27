%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          CIN and B-ALL progression          %
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

% Auxiliary file 1. In this script we predefine
% all parameters and data structures

%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Parameters       %
%%%%%%%%%%%%%%%%%%%%%%%%%%

    simSteps = 100;             % Measured in days (total simulation time)
    dt = 2;
    up_simSteps = dt*simSteps;  % Adjust number of outputs to have the same 
                                % total simulation time   
    time = 0+(1/dt):1/dt:simSteps; % Time vector (for plots)                            
    replicates = 10;            % number of replicates 

   %%%%%% initial condition %%%%%%%
   
   chromosomes = 23; 
   initialCellNumber = 500; 
   totalCells = 10000; % Initial number of cells from which we'll sample
                       % the 
               
    % Carrying capacity for B-ALL cells
    K = 1e6; %[cells/L]
    
    saveKaryotypes = zeros( replicates, 6);

%%%%%%%%%%%%%%%
%    RATES    %
%%%%%%%%%%%%%%%


%     % rhomin and rhomax values taken from 
%     % https://www.sciencedirect.com/science/article/pii/014521269598846P
    rhomin     = 0.2;            % Minimum proliferation probability [day^-1]
    rhomax     = 0.4;            % Maximum proliferation probability [day^-1]          
    M = 500 % Number of different phenotypes [M]
    % Growth rates 
    % rho_init = 0.2161 % fixed
     rho_init = normrnd(0.21,0.001); % Sample from a normal distribution
    % rho_init = rhomin + (rhomax-rhomin).*rand(1,1);

% CIN = probability of missegregation
%     probCIN  = ones(1,length(rho))*0.5*(1/dt);   % Gamma [days^-1]. 
%     ChrDist = 0.1:0.1:1; % Distribution of chromosomes between daughter cells.
%                          % Daughter 2 will get other-daughter1
%     probCIN_val = 0:0.1:0.5;
%     probCIN_mat = (probCIN_val'.*ones(1,length(rho)))';  
      probCIN = 0.2*(1/dt);

% Probability of undergoing whole genome duplication (WGD)
     probWGD = 0.011*(1/dt);                                                                                                                            11; % Mittelmann database

% Death probability. All cells die with the same death rate
% 1/4*rho from initiating phenotype (HalfLatt)
     % 1. CONSTANT PROB DEATH
     % ProbDeath = 0.01*ones(size(rho))*(1/dt); %[days^-1]
     % ProbDeath = ones(size(rho))*rho(HalfLatt)*(1/4); 
     ProbDeath = rho_init/10;     

    
    % 2. SELECTION AGAINST ANEUPLOID ACCORDING TO A POISSON DISTRIBUTION 
    % Generate death probabilities according to a poisson distribution
%     mu_min = 0.0536 %rho_min/4
%     mu_max = 0.223 %rho_max/4
%     lambda = (mu_min*100+mu_max*100)/2; 
%     d = makedist('poisson', 'lambda', lambda); %generate the distribution
%     td = d.truncate(mu_min*100, mu_max*100); % truncate the distribution
%     deathRate = sort(td.random(Nkar,1)/100, 'descend');
%     ProbDeath = deathRate*(1/dt); 
    
    %%% FITNESS/CHROMOSOME %%%
    fitness_hybrid = [0.115324811 0.100896549 0.093934244 0.112698852 0.144263758...
               0.189852198 0.224924724 0.16900021 0.161659944 0.183724525...
               0.306234471 0.30923896 0.262945822 0.484174593 0.469857796...
               0.33827031 0.775702437 0.290919989 0.848513521 0.393448892...
               0.487194168 0.957313303 0.109378302]'; 


    % Generate initial karyotypes 
    basicNumberOfCopies = 2; 
    chromosomes = 23;

    % Selection strength / selection coefficient (sigma) 
    sigma = 0.5; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%    STORING RESULTS    %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Cells matrix (storing karyotypes)
    cells = [];
    %cells = zeros(chromosomes, K, length(saveKaryotypes));
% Tracking total population
    Mass = zeros(up_simSteps, replicates); 
% Tracking average growth rate (whole population)    
    rho_average = zeros( up_simSteps, replicates);
% Tracking average fitness (whole population)
    fitness_average = zeros( up_simSteps, replicates);
% Tracking average fitness (whole population) - NORMALIZED BY TOTAL KAR
    fitness_average_norm = zeros( up_simSteps, replicates);
% Tracking karyotypes in time 
    %cells_karyotype = zeros(chromosomes, K, length(saveKaryotypes), replicates);
    cells_karyotype = zeros(chromosomes, 1e4, up_simSteps, replicates);
% Fitness histograms 
    hist_fitness = zeros(K,up_simSteps,replicates);
% Tracking net karyotype from cells 
    net_karyotypes = zeros(1.2*K, up_simSteps, replicates);
% Tracking whole genome duplication
    nTrack_wgd = zeros( up_simSteps, replicates);   
