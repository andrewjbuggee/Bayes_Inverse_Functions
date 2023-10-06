%% Retreive vertical profiles in a loop
% The loop can be used to change pixels or to change retrieval settings

% This script uses MODIS retrievals rather than running my own TBLUT
% algorithm



% By Andrew John Buggee

%% Load paths

clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;

%% Load MODIS data

% Load modis data and create input structure
clear variables

% Determine which computer you're using
computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';

    % ----------------------------------------
    % ***** Define the VOCALS-REx Folder *****
    % ----------------------------------------

    % ----- November 9th data -----
    vocalsRexFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';
    vocalsRexFile = 'RF11.20081109.125700_213600.PNI.nc';



    % ----- November 11 data -----
%     vocalsRexFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/2008_11_11/';
%     vocalsRexFile = 'RF12.20081111.125000_214500.PNI.nc';






elseif strcmp(computer_name,'andrewbuggee')==true



    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % Define the MODIS folder name

    % ----- November 9th at decimal time 0.611 (14:40) -----
    modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


    % ----- November 11th at decimal time 0.604 (14:30) -----
    %modisFolder = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/'];


    % ----- November 11th at decimal time 0.784 (18:50) -----
    %modisFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';




    % ----------------------------------------
    % ***** Define the VOCALS-REx Folder *****
    % ----------------------------------------



    % ----- November 9 data -----
    vocalsRexFolder = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/';
    vocalsRexFile = 'RF11.20081109.125700_213600.PNI.nc';


end




[modis,L1B_fileName] = retrieveMODIS_data(modisFolder);

% ----- Create a structure defining inputs of the problem -----

modisInputs = create_modis_inputs(modisFolder, L1B_fileName);

%% LOAD VOCALS-REx data

% ---------------------------------------------------------------------
% ------------ Do you want to load VOCALS-REx data? -----------------------
% ---------------------------------------------------------------------


loadVOCALSdata = true;

% % ---------------------------------------------------------------------
% % ----------- Do you want to use VOCALS-REx pixels? -------------------
% % ---------------------------------------------------------------------
% 
% % If flag below is true, only suitable pixels within the VOCALS-REx data set
% % will be used
% useVOCALS_pixelLocations = true;



if loadVOCALSdata==true
    vocalsRex = readVocalsRex([vocalsRexFolder,vocalsRexFile]);

end

%% CROP VOCALS REX DATA

% ------------------------------------------------
% ----- We dont need all of this data ------------
% ------------------------------------------------

% ----- Find all vertical profiles within VOCALS-REx data ------
% find all vertical profiles in the vocals-REx data set. Vertical profiles
% qualify if the total number concentration is above 10^0 and relatively
% flat, the plane is climbing in altitude, and clears the cloud deck from
% bottom to top. If this is all true, save the vocals rex data
lwc_threshold = 0.03;           % g/m^3
stop_at_max_lwc = false;         % truncate profile after the maximum lwc value
Nc_threshold = 1;               % # droplets/cm^3 


% Lets only keep the vocalsRex data we need
% Time is measured in seconds since the startTime

% ---- DO YOU WANT TO USE ADVECTION? -----
modisInputs.flags.useAdvection = true;

tic
vocalsRex = cropVocalsRex2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold, modis, modisInputs);
toc


%% PLOT VOCALS REX DATA WITH MODIS RETRIEVALS

plot_vocalsRex_with_MODIS_retrieved_re(vocalsRex, modis)

%% FIND MODIS PIXELS CLOSEST TO VOCALS

% There may be a time delay, lets subtract a minute from the index found
% above
bufferTime = 0;

% Find the MODIS pixels that overlap with VOCALS-REx

pixels2use = find_MODIS_VOCALS_overlapping_pixels(modis, modisInputs, vocalsRex);


% Tell MODIS_INPUTS to use only these pixels found
modisInputs.pixels.num_2calculate = length(pixels2use.res1km.linearIndex);


%% CREATE GAUSS-NEWTON INPUTS

% We use the estimates calcualted by the TBLUT as our a priori
GN_inputs = create_bayes_inputs(modisInputs);
disp('Dont forget to check the inputs and change if needed!!')

%% CREATE MODEL PRIOR AND COVARIANCE MATRIX AND MEASUREMENT COVARIANCE

% I don't need anything but the covariance matrix and the expected values
%inputs = create_model_prior(inputs,data_inputs);

% -------------------------------------------------------
% do you want to use your estimates or the MODIS estimate?
% -------------------------------------------------------

use_MODIS_estimates = true;

if use_MODIS_estimates==true
    truth_estimate_table = [];
end

GN_inputs = create_model_prior_covariance_andCloudHeight(GN_inputs, pixels2use, truth_estimate_table, use_MODIS_estimates, modis, vocalsRex);
GN_inputs = create_MODIS_measurement_covariance(GN_inputs, modis, modisInputs, pixels2use);


%% CALCULATE RETRIEVAL PARAMETERS

% change the model apriori each loop
r_top_apriori_percentage = 0.1:0.05:0.4;        % percentage of the TBLUT guess
r_bot_apriori_percentage = 0.3:0.05:1;        % percentage of the TBLUT guess
tau_c_apriori_percentage = 0.05:0.05:0.3;        % percentage of the TBLUT guess

tic
for rt = 1:length(r_top_apriori_percentage)
    for rb = 1:length(r_bot_apriori_percentage)
        for tc = 1:length(tau_c_apriori_percentage)
            
            disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])

            % set the new covariance matrix
            % the percentage above multipled by the TBLUT retrieval is the
            % STD. Square it to get the variance
            
            % ----------------------------------------------------------
            % ------- Set the covariance matrix of the FIRST pixel -----
            % ----------------------------------------------------------
            GN_inputs.model.covariance(:,:,1) = diag([(GN_inputs.model.apriori(1,1)*r_top_apriori_percentage(rt))^2,...
                (GN_inputs.model.apriori(1,2)*r_bot_apriori_percentage(rb))^2, (GN_inputs.model.apriori(1,3)*tau_c_apriori_percentage(tc))^2]);


            % ----------------------------------------------------------
            % ------- Set the covariance matrix of the MEDIAN pixel -----
            % ----------------------------------------------------------
            GN_inputs.model.covariance(:,:,2) = diag([(GN_inputs.model.apriori(2,1)*r_top_apriori_percentage(rt))^2,...
                (GN_inputs.model.apriori(2,2)*r_bot_apriori_percentage(rb))^2, (GN_inputs.model.apriori(2,3)*tau_c_apriori_percentage(tc))^2]);


            % ----------------------------------------------------------
            % ------- Set the covariance matrix of the LAST pixel -----
            % ----------------------------------------------------------
            GN_inputs.model.covariance(:,:,3) = diag([(GN_inputs.model.apriori(3,1)*r_top_apriori_percentage(rt))^2,...
                (GN_inputs.model.apriori(3,2)*r_bot_apriori_percentage(rb))^2, (GN_inputs.model.apriori(3,3)*tau_c_apriori_percentage(tc))^2]);
            
            % Compute the retrieval variables
            [GN_outputs, GN_inputs] = calc_retrieval_gauss_newton_4modis(GN_inputs,modis,modisInputs,pixels2use);
            
            save([modisInputs.savedCalculations_folderName,'GN_inputs_outputs_rt_',num2str(rt),'_rb_', num2str(rb),...
                '_tc_', num2str(tc), '_',char(datetime("today")),'.mat'],"GN_outputs","GN_inputs");

        end
    end
end

toc
%% PLOT RETRIEVED VERTICAL PROFILE WITH MODIS RETRIEVAL

plot_vocalsRex_with_MODIS_retrieved_re_and_vertProf_retrieval(vocalsRex, modis, GN_outputs)

