% ---- Compute fowrad model using MODIS spectral channels ----

% this function will compute reflectances using LibRadTran radiative
% transfer solver. It will do this for a water cloud with a droplet
% profile over MODIS spectral channels

% for now we will run this using the first 7 MODIS channels

% By Andrew J. Buggee
%%
function measurement_estimate = compute_forward_model_4modis(modis,current_guess,GN_inputs,pixel_row,pixel_col,modisInputs)

% Define some needed folder and file names
saveCalculations_fileName = GN_inputs.save_calcs_fileName; % where to save the computed data
INP_folderName = modisInputs.INP_folderName; % Where to save the INP files

% --- compute the forward model at our current estimate ---
r_top = current_guess(1);
r_bottom = current_guess(2);
tau_c = current_guess(3);

profile_type = GN_inputs.model.profile.type; % type of water droplet profile

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = modisBands(1);    % nm - Wavelength used for cloud optical depth calculation
% ----------------------------------------------------------

% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------

H = 0.5;                                % km - geometric thickness of cloud
n_layers = 10;                           % number of layers to model within cloud

% Cloud top
z_top = modis.cloud.topHeight(pixel_row, pixel_col)/1e3;        % meters -  cloud top height

%z0 = 0.9;                                 % km - base height of cloud
z = linspace(z_top-H, z_top,n_layers);        % km - altitude above ground vector
indVar = 'altitude';                    % string that tells the code which independent variable we used
constraint = profile_type;              % string that tells the code which physical constraint to use

re = create_droplet_profile2([r_top, r_bottom], z, indVar, constraint);     % microns - effective radius vector

dist = 'mono';                         % droplet distribution
homogenous_str = 'non-homogeneous';     % This tells the function to create a multi-layered cloud
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector. 

parameterization_str = modisInputs.flags.wc_parameterization;

wc_filename = write_wc_file(re, tau_c, z_topBottom, wavelength_tau_c(1,1), dist, homogenous_str, parameterization_str);


% ----- Write an INP file --------
GN_names.inp = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename);
    
% now lets write the output names
    
GN_names.out = writeOutputNames(GN_names.inp);

% ---- Run uvspec for the files created -----
[measurement_estimate,~] = runReflectanceFunction_4gaussNewton(GN_names,INP_folderName,saveCalculations_fileName);






end