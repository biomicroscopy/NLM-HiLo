# NLM-HiLo

> This content is based on our bioRxiv preprint titled ["Multiplane HiLo microscopy with speckle illumination and non-local means denoising"](https://www.biorxiv.org/content/10.1101/2023.09.08.556851v1).


### Background

[HiLo microscopy](https://sites.bu.edu/biomicroscopy/research/hilo/) is a variant of structured illumination microscopy. It makes use of an uniform-illumination image and a structured illumination image (with a grid or speckle pattern) to achieve widefield optical sectioning.

Implementation of basic speckle-based HiLo reconstruction and non-local means (NLM) denoising to reduce residual speckle noise in the processed image is provided here.

### Instructions

 - Set up environment:

MATLAB 2021a 

Visual Studio 2019 

Check [here](https://www.mathworks.com/support/requirements/previous-releases.html) to see if your version of Visual C++ is supported in MATLAB.

 - Run HiLo example
```Matlab
main_test.m
```
A single image $(571\times 571)$ with a search window width of $31\times 31$ should take ~9s to process.

### Contact
Please reach out to <sqzheng@bu.edu> with any questions.
