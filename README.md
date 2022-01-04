# CLSM spectral analysis

Software package with graphical user interface for plotting and analysis of spectral image (stacks).

## Documentation

The documetation for the graphical user interface is work in progress.

The documentation for scripting can be called by typing in the MATLAB command line
```
doc CLSMspectraLab
```
for the docuementation of the background correction, type
```
doc bgCorrection
```
and for the documentation of the mask, type
```
doc imageMask
```

The CLSMspectraLab class heavily relies on the imageStack class, which documentation can be called when in the MATLAB path via
```
doc imageStack
```

## Getting started

To load an image via the BioFormats reader into the CLSMspectraLab class 'cl', type 
```
cl = CLSMspectraLab.import('path/to/my/file.nd2')
```
Make sure the BioFormats package is in the MATLAB path and the 'bioformats_package.jar' is in the java path, see https://nl.mathworks.com/help/matlab/ref/javaaddpath.html. A list of the supported formats can be found here: https://docs.openmicroscopy.org/bio-formats/5.8.2/supported-formats.html and is not limited to .nd2 files.
The image is stored in an instance of the imageStack class in cl.input. Now plot the stored image directly via the imageStack functionality
```
figure; cl.input.plotImage;
```
and plot the channels (spectrum)
```
figure; cl.plotChannels();
```

The background correction is stored in cl.bgCorrection and can be changed to a different child of the bgCorrection class (see src/+bgCorrection) and its properties can be modified. To run the default settings
```
cl = cl.computeBgCorrection;
```
Similarly for the mask, make sure a child of the bgCorrection class (see src/+imageMask) is loaded in cl.mask and compute. Please note that the background substraction only will be considered in the computation of the mask when the cl.bgCorrectionFlag = true. Compute
```
cl = cl.computeMask;
```
To make sure the mask is applied, make sure that cl.maskFLag = true. Now, the summed intensity along x,y can be plotted with the mask applied
```
figure; cl.plotSumIntensity2D
```

For more advanced use of the class, please check the documentation and source code of the GUI (see: src/GUI/CLSMspectraLabGUI.mlapp).

## Project organization
- PG = project-generated
- HW = human-writable
- RO = read only
```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)
