# NLM-HiLo
HiLo microscopy reconstruction with NLM denoising

### Background

Implementation of basic speckle-based HiLo reconstruction and non-local means denosing to reduce speckle artifact.

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
