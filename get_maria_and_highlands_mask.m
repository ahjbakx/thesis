function [maria_mask, highlands_mask] = get_maria_and_highlands_mask()
    latitude = fitsread("/Users/aaron/thesis/Data/polarisation/latitude.fits");
    longitude = fitsread("/Users/aaron/thesis/Data/polarisation/longitude.fits");
    albedo = load_fits("Av");
    maria_mask = load("maria_mask.mat").maria_mask;
    highlands_mask = (latitude>-99 & albedo>-99) & not(maria_mask); % filter NaN and infs
end