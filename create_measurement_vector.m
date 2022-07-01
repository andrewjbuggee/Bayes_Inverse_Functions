% --- Create the measurement vector needed for Rodgers solution of gaussian
% model pdf and measurement pdf ---


% By Andrew J. Buggee

%%

function y = create_measurement_vector(modis, GN_inputs, pixels2use)


% We will retireve cloud properties for n pixels, designated by the
% Gauss_Newton input structure
num_pixels = GN_inputs.numPixels2Calculate;

% Right now, use all 7 bands
y = zeros(size(modis.EV1km.reflectance,3),num_pixels);

% For now we use all 7 bands
for pp = 1:num_pixels
    
    
    row = pixels2use.res1km.row(pp);
    col = pixels2use.res1km.col(pp);
    
    y(:,pp) = reshape(modis.EV1km.reflectance(row,col,:),[],1);
    
end


end