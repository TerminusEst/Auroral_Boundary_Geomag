# Auroral_Boundary_Geomag
Python scripts which were used to calculate auroral boundaries from geomagnetic data, as appears in Blake et al. (2020) *Estimating Location of Auroral Equatorward Boundary using Historical and Simulated Surface Magnetic Field Data*, JGR Space Physics, *under review*, DOI: XXXXXXXXXXXX

There are two scripts:
1. **EField_Calc.py** calculates the E-field for a single [INTERMAGNET](https://intermagnet.org/) site using the [Quebec 1-D resistivity model](https://doi.org/10.1046/j.1365-246x.1998.00388.x).
2. **Boundary_Calc.py** calculates the extent of the auroral boundary using magnetic latitudes and maximum calculated Eh values from multiple INTERMAGNET sites.

## Dependencies and Versions used
- [gcvspline](https://github.com/charlesll/gcvspline) (0.4)
- numpy (1.18.4)
- matplotlib (3.1.1)
- math
- datetime (4.3)

## Script Examples
**EField_Calc.py**
Takes INTERMAGNET 60s B-Field data, and calculates E-Field data using a 1-D resistivity profile. This example is for 3 days of data from the Eskdalemuir site (found in the Data/ folder).
![Figure_2](https://user-images.githubusercontent.com/20742138/91243809-d232d700-e718-11ea-89df-23bcc4f6e114.png)

**Boundary_Calc.py**
This script takes in magnetic latitudes and maximum E-field values found in the Data/ folder (2003_10_30.txt and 2009_10_07.txt).
Natural cubic splines are fit to these data. The boundary is taken as the point on this fit that has the largest gradient:

![Figure_1](https://user-images.githubusercontent.com/20742138/91241931-0f489a80-e714-11ea-8c8b-d77a3d3cdff2.png)

And here are those boundaries for the non-log10 versions of the above plot:
![Figure_3](https://user-images.githubusercontent.com/20742138/91332316-00a4c680-e79a-11ea-87e3-ed9f24dadd9c.png)

## Author
Written by Sean Blake in NASA GSFC, 2018-2020

Email: blakese@tcd.ie
