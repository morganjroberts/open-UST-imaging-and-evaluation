# open-UST-imaging-and-evaluation
Evaluating the acoustic performance and imaging capability of the [open-UST](https://github.com/morganjroberts/open-UST) ultrasound tomography system.

This repo contains code for analysing experimental data from the open-UST system.

## Accessing Experimental Data
The experimental data files are too large to store on GitHub. If you are external from the Biomedical Ultrasound Group, please get in contact to request access to the data.

Internal BUG users: the data is stored in `Datasets\open-UST-imaging-and-evaluation`. 

**Note:** Do not copy the data folder into your local copy of this repository - you risk accidentally adding them to git.

## Getting Started

Open a terminal and change directories to an appropriate folder:
```
cd C:\Users\Morgan\Documents\PhD Git Repos
```
Clone this repository:
```
git clone https://github.com/morganjroberts/open-UST-imaging-and-evaluation.git
```
Open MATLAB, and open the function `addRepoDataPath.m`. Edit the absolute path to the repository and the data directory:
```
function [repo_dir, data_dir] = getRepoDataPath()

repo_dir = 'C:\Users\Morgan\Documents\PhD Git Repos\open-UST-imaging-and-evaluation';
data_dir = 'Z:\open-UST-imaging-and-evaluation';

end
```
**Note:** Do not include a file separator (`/` or `\`) at the end of either path.

This function is used frequently throughout the repository to define absolute paths to scripts and data files that are outside of the current working directory.

Add this function to the Matlab path by typing `pathtool` at the MATLAB command line. This dialog box can also be accessed using the `Set Path` button on the ribbon bar. Once the dialog box is open, click `Add with Subfolders`, select the `open-UST-imaging-and-evaluation` folder, and click `Save`.

## Requirements

- Matlab 2022a or more recent
- Matlab toolboxes:
  - Signal processing
  - Curve fitting
- [k-Wave](https://github.com/ucl-bug/k-wave)
- [ust-sart](https://github.com/ucl-bug/ust-sart)