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

      a. Finite Difference (Real Space)
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
If you want to download a specific file from Github using your command line you may:

* Open the file in Github.

* Copy the specific file's raw link from GitHub. On the top right corner you can see the option to open the file in raw mode.

* Open it in raw mode and copy the URL.

* Now use curl command in command line to download the file:

      curl -o filename raw-link-to-file
      

But, if you want to download the entire repository, you may:

* Copy the URL where appear all the files. Make sure you are not in any branch.

* Once you have chosen a directory, in the command line type:

      git clone <Repository URL>
      
Now you have in your directory all the files from your Github repository.
      
 




