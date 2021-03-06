# CompartmentalAnalysis
Matlab Graphic User Interface (GUI) for compartmental analysis of various compartmental systems describing FDG uptake in murine models.

The code is based on the following publications:

*2018 F. Delbary and S. Garbarino "Compartmental analysis of dynamic nuclear medicine data: regularization procedure and application to physiology", Inverse Problems in Science and Engineering 1--19.

*2016 F. Delbary, S. Garbarino and V. Vivaldi "Compartmental analysis of dynamic nuclear medicine data: models and identifiability", Inverse Problems 32 125010


A GUI is provided for the following compartmental systems:

1) a standard 2D compartmental model for reversible FDG uptake. See for instance [2002 Gunn, R.N., Gunn, S.R., Turkheimer, F.E., Aston, J.A. and Cunningham, V.J. "Positron emission tomography compartmental models: a basis pursuit strategy for kinetic modeling", Journal of Cerebral Blood Flow & Metabolism, 22(12), pp.1425-1439] for oncological applications.

2) a (3+1)D compartmental model for FDG kinetic in the kidneys [2014 Garbarino S, Caviglia G, Sambuceti G, Benvenuto F and Piana M "A novel description of FDG excretion in the renal system: application to metformin-treated models" Physics in Medicine and Biology, 59, 2469-2484]

3) a 2D-2input compartmental model for FDG kinetics in the liver [2015 Garbarino S, Vivaldi V, Delbary F, Caviglia G, Piana M, Marini C, Capitanio S, Calamia I, Buschiazzo A and Sambuceti G, "A new compartmental method for the analysis of liver FDG kinetics in small animal models" European Journal of Nuclear Medicine and Molecular Imaging Research, 2015, 5-35].


# Usage:
Code is written in Matlab R2015b and tested with versions up to R2019b.

run main.m file launches a Matlab GUI for tumor, kidneys or liver compartmental modelling. It requires standard .voistat data. 
Test data can be found in the "data test" folder. Pictorial representation of the 3 compartmental models implemented (2D, (3+1)D and 2D2input) are available in the "models" folder, and accessible from the GUI itself.
