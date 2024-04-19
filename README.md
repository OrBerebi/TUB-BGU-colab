# TUB-BGU-colab
Enhancing spatial cues for low-order HRTF SH representation

## Matlab
TUB-BGU-colab/main.m script generates a .mat file that contains SH coefficients (and other parameters) from a .sofa file

## python
TUB-BGU-colab/PyTorch/main.ipynb notebook is the iMagLS implementation for the HRTF SH coefficients solution. its loads the .mat file from TUB-BGU-colab/main.m and calculate the iMagLS solution, also can save static and dynamic source .wav binaural audio files and save the SH coefficients as a space/time .sofa files

TUB-BGU-colab/PyTorch/baumgartner2014/run_demo.ipynb a python implemention of the amtoolbox baumgartner2014 sound-source localization in sagittal plane model
