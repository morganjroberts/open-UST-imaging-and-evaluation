# open-UST-imaging-and-evaluation
Evaluating the acoustic performance and imaging capability of the [open-UST](https://github.com/morganjroberts/open-UST) ultrasound tomography system.

This repo contains code for analysing experimental data from the open-UST system.

## Accessing Experimental Data
The experimental data files are too large to store on GitHub. If you are external from the Biomedical Ultrasound Group, please get in contact to request access to the data.

Internal BUG users: the data is stored in `ScanData\open-UST-imaging-and-evaluation`. When running these scripts in Matlab, first map the `ScanData` network drive to your PC and add the data subfolder to the path, for example:
```
addpath(genpath("Y:\open-UST-imaging-and-evaluation"))
```

## Requirements

- Matlab 2022a or more recent
- Matlab toolboxes:
  - Signal processing
  - Curve fitting
- [k-Wave](https://github.com/ucl-bug/k-wave)
- [ust-sart](https://github.com/ucl-bug/ust-sart)