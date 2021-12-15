% ---- Compute fowrad model using MODIS spectral channels ----

% this function will compute reflectances using LibRadTran radiative
% transfer solver. It will do this for a water cloud with a droplet
% profile over MODIS spectral channels

% for now we will run this using the first 7 MODIS channels

% By Andrew J. Buggee
%%
function measurement_estimate = compute_forward_model_4modis(modis,retrieval,bayes_inputs,pixel_row,pixel_col,INP_folderName,saveCalculations_fileName)

% --- compute the forward model at our current estimate ---
r_top = retrieval(1);
r_bottom = retrieval(2);
tau_c = retrieval(3);

profile_type = bayes_inputs.model.profile.type; % type of water droplet profile

% ----------------------------------------------------------

% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------

wc_profile_fileName = write_waterCloud_files_with_profile(r_top,r_bottom,tau_c,profile_type);

% ----- Write an INP file --------
names.inp = write_INP_4_MODIS_singlePixel_withProfile(wc_profile_fileName,INP_folderName,modis,pixel_row,pixel_col,retrieval);
    
% now lets write the output names
    
names.out = writeOutputNames(names.inp);

% ---- Run uvspec for the files created -----
[measurement_estimate,~] = runReflectanceFunction_4gaussNewton(names,INP_folderName,saveCalculations_fileName);






end