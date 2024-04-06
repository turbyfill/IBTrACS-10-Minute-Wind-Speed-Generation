# NA-Wind-Speed-Generation

## Description

The script included in this repository is a master program that complements the `holland_model.py`. It processes tropical cyclone track data to estimate wind speeds across different geographic regions. Specifically, the script accepts data from the International Best Track Archive for Climate Stewardship (IBTrACS) and applies the Holland model to estimate 10 minute wind speeds.

### What the Script Does

- **Data Preprocessing**: The script starts by loading and preprocessing IBTrACS data, focusing on relevant attributes like storm identifiers, locations, pressures, and wind speeds.
- **Storm Movement Speed Calculation**: It calculates the movement speed of each storm, which is crucial for accurate wind field estimation.
- **Holland Model Application**: Using the `holland_model.py` module, the script applies the Holland wind model to estimate the radial wind field around the tropical cyclone's center.
- **Wind Field Estimation**: It generates a spatial grid over the area of interest and calculates the wind speed at each grid point based on the cyclone's characteristics and the Holland model output.
- **Data Aggregation**: The script aggregates the wind field data, providing insights into the maximum wind speeds experienced in different regions during the lifetime of a tropical cyclone.

This tool is essential for researchers, meteorologists, and climate scientists who require detailed estimates of wind speeds from historical tropical cyclone tracks for various applications, including risk assessment and climate studies.

## Background and Acknowledgements

This project is inspired by and derived from the research conducted by Bloemendaal, N., Haigh, I.D., de Moel, H., and others, as detailed in their work:

- Bloemendaal, N., Haigh, I.D., de Moel, H. et al. "Generation of a global synthetic tropical cyclone hazard dataset using STORM." Sci Data 7, 40 (2020). [https://doi.org/10.1038/s41597-020-0381-2](https://doi.org/10.1038/s41597-020-0381-2)

Additionally, the methodology adopts the principles from Lin, N., and Chavas, D.'s research:

- Lin, N., and Chavas, D. (2012), "On hurricane parametric wind and applications in storm surge modeling," J. Geophys. Res., 117, D09120, [doi:10.1029/2011JD017126](https://doi.org/10.1029/2011JD017126).

The original script that this work is based on can be found at:

- [https://github.com/NBloemendaal/STORM-return-periods/blob/master/masterprogram.py](https://github.com/NBloemendaal/STORM-return-periods/blob/master/masterprogram.py)

## Dependencies

- Python
- Pandas
- NumPy
- SciPy
- holland_model.py: https://github.com/NBloemendaal/STORM-return-periods/blob/master/holland_model.py 

## Usage

Place the `holland_model.py` in the same directory as the main script. Ensure that the necessary data files are in the correct format and location as expected by the script. Run the main script using Python:

```bash
MASTER.py
