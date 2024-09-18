## Usage:

Clone the environment and cd into the directory and run:

```
python example_run.py
```

## Notes:
example_run.py can be edited with the partial pressures of each gas, the final spectral resolution, the temperature, and wavelength range. Note that the code has been validated against other gas cell simulation codes, but enough edits have happened now that the output should be re-validated using https://hitran.iao.ru/gasmixture/simlaunch or http://spectraplot.com/absorption 

## Details
Since hapi assumes pressure broadening is just due to air broadening, the code makes some simplifications in generating the combined gas cell. It works by generating 1D transmittance models for individual gas species at the total pressure of the final gas cell and a cell length of L_0=1cm (this is done using the HITRAN hapi.py code, documentation [here](https://hitran.org/static/hapi/hapi_manual.pdf).  These files are saved and then reloaded and are multiplied together after raising each molecule's transmittance values to the power of (p_fraction * L_tot/L_0), where p_fraction is the ratio of the gas's partial pressure in the cell to the total cell pressure, and L_tot is the total gas cell length in cm. This basically scales to the amount of the gas in the cell.




