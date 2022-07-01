% --- Compute Jacobian matrix for MODIS channels-----



% By Andrew J. Buggee
%%

function jacobian = compute_jacobian_4modis(modis,state_vector,measurement_estimate,GN_inputs, modisInputs, pixel_row,pixel_col, measurement_variance)

% --- Define the filename to save all calculations ---
saveCalculations_fileName = GN_inputs.save_calcs_fileName;

% --- Define the INP Folder location ---
INP_folderName = modisInputs.INP_folderName;

% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% ---- define the incremental change to each variable -----

change_in_state = [0.05 * r_top, -0.05 * r_bottom, 0.05 * tau_c];


% ----------------------------------------------------------

% Set up a few constants for the water cloud
H = 0.50;                                % km - geometric thickness of cloud
n_layers = 10;                          % number of layers to model within cloud
    
z0 = 0.9;                                 % km - base height of cloud
z = linspace(z0, z0+H,n_layers);        % km - altitude above ground vector
indVar = 'altitude';                    % string that tells the code which independent variable we used

profile_type = GN_inputs.model.profile.type; % type of water droplet profile
num_model_parameters = GN_inputs.num_model_parameters;
dist = 'mono';                         % droplet distribution
homogenous_str = 'non-homogeneous';     % This tells the function to create a multi-layered cloud
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector.

parameterization_str = modisInputs.flags.wc_parameterization;

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = modisBands(1);    % nm - Wavelength used for cloud optical depth calculation

% Lets step through each model variable and compute the derivative
jacobian = zeros(length(measurement_estimate),num_model_parameters);
change_in_measurement = zeros(length(measurement_estimate),num_model_parameters);

for xx = 1:num_model_parameters
    
    % We start with the original state vector
    newState_vector = state_vector;
    
    % Then we alter just one of them
    newState_vector(xx) = state_vector(xx) + change_in_state(xx);
    
    new_r_top = newState_vector(1);
    new_r_bottom = newState_vector(2);
    new_tau_c = newState_vector(3);
    % --------------------------------------------
    % create water cloud file with droplet profile
    % --------------------------------------------
    
    re = create_droplet_profile2([new_r_top, new_r_bottom], z, indVar, profile_type);     % microns - effective radius vector
    
    
    wc_filename = write_wc_file(re, new_tau_c, z_topBottom, wavelength_tau_c(1,1), dist, homogenous_str, parameterization_str);
    
    
    % ----- Write an INP file --------
    names.inp = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
    
    % ---- Run uvspec for the files created -----
    [new_measurement_estimate,~] = runReflectanceFunction_4gaussNewton(names,INP_folderName,saveCalculations_fileName);
    
    change_in_measurement(:,xx) = new_measurement_estimate' - measurement_estimate;

    jacobian(:,xx) = change_in_measurement(:,xx)./change_in_state(xx);


    
end


% --- Optional Plot! ---

spectral_bands = modisBands(1:7);
[~, index_sort] = sort(spectral_bands);
string_bands = string(round(spectral_bands(index_sort(:,1),1)));


f = figure; bar([abs(change_in_measurement(index_sort(:,1),:)),sqrt(measurement_variance(index_sort(:,1)))])
hold on
xticklabels(string_bands);
xlabel('Wavelength $(nm)$', 'Interpreter','latex')
ylabel('$\triangle$ Reflectance','Interpreter','latex')
legend('$\triangle r_{top}$','$\triangle r_{bot}$', '$\triangle \tau_{c}$','$\sigma_\lambda$',...
     'interpreter', 'latex', 'Location','best','Fontsize',20); 
grid on; grid minor
set(f, 'Position',[0 0 1000 500])
title('The Jacobian', 'Interpreter','latex')
dim = [.14 0.67 .3 .3];
str = ['$r_{top} = $',num2str(state_vector(1)),', $r_{bot} = $ ',num2str(state_vector(2)),', $\tau_c = $ ',num2str(state_vector(3))];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
    'FontWeight','bold','FontSize',14, 'EdgeColor','w','Interpreter','latex');



%-------------------------------------------------------------
% ---- SPECIAL PLOT FOR r_bot --------------------------------
%-------------------------------------------------------------

% spectral_bands = modisBands(1:7);
% [~, index_sort] = sort(spectral_bands);
% string_bands = string(round(spectral_bands(index_sort(:,1),1)));
% 
% load('jacobian_rt-10_rb-9_tau-20.mat', 'change_in_measurement')
% change_r_bot = change_in_measurement(index_sort(:,1),2);
% load('jacobian_rt-10_rb-9_tau-15.mat', 'change_in_measurement')
% change_r_bot = [change_r_bot, change_in_measurement(index_sort(:,1),2)];
% load('jacobian_rt-10_rb-9_tau-10.mat', 'change_in_measurement','measurement_variance')
% change_r_bot = [change_r_bot, change_in_measurement(index_sort(:,1),2)];
% 
% f = figure; bar([abs(change_r_bot),sqrt(measurement_variance(index_sort(:,1)))])
% hold on
% xticklabels(string_bands);
% xlabel('Wavelength $(nm)$', 'Interpreter','latex')
% ylabel('$\triangle$ Reflectance','Interpreter','latex')
% legend('$\tau_c = 20$','$\tau_c = 15$', '$\tau_c = 10$','$\sigma_\lambda$',...
%      'interpreter', 'latex', 'Location','best','Fontsize',20); 
% grid on; grid minor
% set(f, 'Position',[0 0 1000 500])
% title('$\partial F(\vec{x})/\partial r_{bot}$', 'Interpreter','latex')






end