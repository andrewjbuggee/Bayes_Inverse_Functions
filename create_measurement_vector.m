% --- Create the measurement vector needed for Rodgers solution of gaussian
% model pdf and measurement pdf ---


% By Andrew J. Buggee

%%

function y = create_measurement_vector(modis,data_inputs)

pixels2use = data_inputs.pixels2use;
refl_data = cat(3,modis.EV.m250.reflectance,modis.EV.m500.reflectance);
bands2run = data_inputs.inputs.bands2run;
num_pixels = data_inputs.inputs.pixels.num_2calculate;

% SWITCH TO 1KM DATA
% but for now, lets propose a simple algebraic fix so taht we select
% the proper pixel in our 500 meter data set
pixelIndex_500m = 1:4:((num_pixels-1)*4 +1);


y = zeros(length(bands2run),num_pixels);

for pp = 1:num_pixels
    
    
    row = pixels2use.res500m.row(pixelIndex_500m(pp));
    col = pixels2use.res500m.col(pixelIndex_500m(pp));
    
    y(:,pp) = reshape(refl_data(row,col,bands2run),length(bands2run),1);
    
end


end