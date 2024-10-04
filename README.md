## Usage:

Clone the environment and cd into the directory and run:

```
python example_run.py
```

## Notes:
example_run.py can be edited with the partial pressures of each gas, the final spectral resolution, the temperature, and wavelength range. Note that the code has been validated against other gas cell simulation codes, but enough edits have happened now that the output should be re-validated using https://hitran.iao.ru/gasmixture/simlaunch or http://spectraplot.com/absorption 

## Details
Since hapi assumes pressure broadening is just due to air broadening, the code makes some simplifications in generating the combined gas cell. It works by generating 1D transmittance models for individual gas species at the total pressure of the final gas cell and a cell length of L_0=1cm (this is done using the HITRAN hapi.py code, documentation [here](https://hitran.org/static/hapi/hapi_manual.pdf).  These files are saved and then reloaded and are multiplied together after raising each molecule's transmittance values to the power of (p_fraction * L_tot/L_0), where p_fraction is the ratio of the gas's partial pressure in the cell to the total cell pressure, and L_tot is the total gas cell length in cm. This basically scales to the amount of the gas in the cell.


3 - run code - this saves it to out_path with option of plotting
	run(species,pres,temp,path_length,out_path,ploton=True)


* I also built in compare_sims which compares HAPI results to SpectralPlot online app and shows that the results are slightly different
* if you want multiple in one cell at total pressure P, evaluate each individual spectrum at P, but wth path_length equal to the fraction of the gas desired e.g. 1/3 of each gas at P=2, use P=2, with path_length = 2 * total_length/3 for each gas
* will probably want to save the hapi sims to file..
*if you change the wavelength range then you should delete the hitran files and rerun the code 
so that it initiates a download again with the correct 'nu' data. if no transitions exist
in a region then an error will occur


4/1/21 Note: I took the gh version on april 1 2021 to make plotspec_gh then updated the plot_wl function so that it convolves the high resolution spectrum to R~30000 or whatever res is wanted. plotspec2 may be pointless, plotspec and plotspec2 should be similar and have the code to open the gas cell. data without an isotopologue [1] number in the name are wrong at low resolution. I must have changed the resolution. So the lines looked deeper if using those files. plotspec_gh.py is most reliable although the change of resolution fix isn't a big deal. its just more accurate to convolve after adjusting the high res spectrum, but the errors for the purposes here were small.

To do: should consolidate the files/functions. then reupload to github and make a development branch and be happy w/ master branch

