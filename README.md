# bandstructureforgpaw
This script allows you to interpolate the different bands in eigenenergies vs. k-point graphs. You must have a previous calculation with GPAW for using it.

## Instructions

Provide a .gpw file and run the program. 

$ python band_structure.py file.gpw

Then, the program will ask you for the desired method of interpolation: LCAO, PAW, PW and Finite Difference.

```
#########################################################
################## Band Structure Plotter ###############
#########################################################

>>> Which method should I use for calculating your band structure?

      a. Finite difference(Real Space)
      b. Plane Wave Expansion
      c. Linear Combination of Atomic Orbitals (LCAO)
      d. Projector Augmented Wave

>>>
```
Then, it lets you display the plot, save the data to a file or run again the calculations.

```
>>>>     What do you want to do now?


        a. Display the plot of your band structure
        b. Export your band structure data into a .txt file
        c. Redo the calculations
        d. Close the program

>>> 

```




