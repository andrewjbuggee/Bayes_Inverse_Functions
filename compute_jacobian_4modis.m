% --- Compute Jacobian matrix for MODIS channels-----



% By Andrew J. Buggee
%%

function jacobian = compute_jacobian_4modis(modis,state_vector,measurement_estimate,bayes_inputs,pixel_row,pixel_col,INP_folderName,saveCalculations_fileName)


% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% ---- define the incremental change to each variable -----

change_in_state = [0.01 * r_top, 0.01 * r_bottom,0.01 * tau_c];


% ----------------------------------------------------------

profile_type = bayes_inputs.model.profile.type; % type of water droplet profile
num_model_parameters = bayes_inputs.num_model_parameters;


% Lets step through each model variable and compute the derivative
jacobian = zeros(length(measurement_estimate),num_model_parameters);

for xx = 1:num_model_parameters
    
    newState_vector = state_vector;
    newState_vector(xx) = state_vector(xx) + change_in_state(xx);
    
    new_r_top = newState_vector(1);
    new_r_bottom = newState_vector(2);
    new_tau_c = newState_vector(3);
    % --------------------------------------------
    % create water cloud file with droplet profile
    % --------------------------------------------
    
    wc_profile_fileName = write_waterCloud_files_with_profile(new_r_top,new_r_bottom,new_tau_c,profile_type);
    
    % ----- Write an INP file --------
    names.inp = write_INP_4_MODIS_singlePixel_withProfile(wc_profile_fileName,INP_folderName,modis,pixel_row,pixel_col,newState_vector);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
    
    % ---- Run uvspec for the files created -----
    [new_measurement_estimate,~] = runReflectanceFunction_4gaussNewton(names,INP_folderName,saveCalculations_fileName);
    
    jacobian(:,xx) = (new_measurement_estimate' - measurement_estimate)./change_in_state(xx);
    
end





end