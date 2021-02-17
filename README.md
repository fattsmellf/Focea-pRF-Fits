# Focea-pRF-Fits
MATLAB code to fit wide-field pRF mapping data

This example code will load one session from a bilaterally imaged
clear-skull mouse. The evoked response (DeltaF/F) has been calculated for
the 800x800 pixel image.

This mouse was tested with 31 bar stimuli, the predicted responses to the
31 stimuli from the set of Gaussians used in the paper is loaded from the
file PredResp_WF31.

The pRF model is then fit to the data using linear regression, this takes
some time (~30 mins) due to the large number of Gaussians we take a number of steps to improve processing speed;
i)     We only fit every other pixek
ii)    The regression model is fit for all pixels simultaneously using lscov
iii)   The data is split into 12 slices and processed in parallel on a multi-core machine.
iv)    Only the minimum value of SS is stored and updated

Tested in MATAB 2019b, but should work in earlier versions that support
parallel processing.

     Copyright (C) 2021 M.W.Self

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
     http://www.gnu.org/licenses/
