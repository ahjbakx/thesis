function fits = load_fits(name)

    fits = fitsread( strcat("/Users/aaron/thesis/Data/polarisation_data/", name, ".fits") );

    fits_size = size(fits);
    for i = 1:fits_size(1)
        for j = 1:fits_size(2)
           if fits(i,j) == -99
               fits(i,j) = NaN;
           end
        end
    end
    
   fits = fits;

end