# Cell-adhesion-on-fluid-membranes
This public repository contains the MATLAB codes to make theoretical predictions and to generate the figures in the Supplemental Material for our article:

**Cell adhesion and spreading on fluid membranes through microtubules-dependent mechanotransduction**

Oleg Mikhajlov, Ram M. Adar, Maria Tătulea-Codrean, Anne-Sophie Macé, John Manzi, Fanny Tabarin, Aude Battistella, Fahima di Federico, Jean-François Joanny, Guy Tran van Nhieu, Patricia Bassereau

doi: https://doi.org/10.1101/2022.09.12.507658

## Contents

The folder contains these MATLAB routines:

* **figF2.m** - plot Figure F2 from Supplemental Material (corresponds to Supp. Fig. S13 from main manuscript)
* **figF3.m** - plot Figure F3 from Supplemental Material (corresponds to Fig. 4F from main manuscript)
* **FigF5_X.m** - plot Figure F5 panel X (X being from A to H) (corresponds to Supp. Fig. S14 from main manuscript)

The folder contains these MATLAB functions:

* **findturningpoints.m**
This function takes as input the parameters of the theoretical model and a given cytoskeletal force f (measured in pN) and returns:
(1) the intervals of phib between 0 and 1 where the effective chemical potential is monotonically increasing or decreasing, and
(2) the values of the chemical potential at the ends of those intervals.


* **findsteadystates.m**
This function takes as input the parameters of the theoretical model, 
a given cytoskeletal force f (measured in pN) and chemical potential
(measured in kBT), and it returns the solution for the bond density
within the interval indicated by the input parameter phibinterval.

## How to use

Run each file **figFX.m** to plot Figure FX from the Supplemental Material file. These routines call the functions **findturningpoints.m** and **findsteadystates.m** for the theoretical calculations. 
