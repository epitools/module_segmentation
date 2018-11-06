# GUI-Independent segmentation module from EpiTools

this repository contains the segmentation module from EpiTools-Matlab ([/src/src_analysis/module_segmentation](https://github.com/epitools/epitools-matlab/tree/master/src/src_analysis/module_segmentation)) with all the GUI dependencies removed.

The main segmentation function can be found in `SegmentIm.m` which calls the c routines `findcellsfromregiongrowing.cc` and `growcellsfromseeds3.cc` for which all the compiled MEX files (win/lin/mac) are included as well.

To know more about the segmentation routine check out our wiki page: [epitools.github.io/wiki/Analysis_Modules/03_segmentation/](https://epitools.github.io/wiki/Analysis_Modules/03_segmentation/) or the epitools article [http://dx.doi.org/10.1016/j.devcel.2015.12.012](http://dx.doi.org/10.1016/j.devcel.2015.12.012)

# example execution

```matlab
%% load example image (Projected Drosophila Wing Disc - Ecad:GFP)
load('ProjIm.mat')

%% crop image for testing (right click -> "Crop Image")
[crop, rect] = imcrop(ProjIm,[]);
close();

%% compute segmentation and output segmentation feedback
s1=0.5;     % first smoothing (seeding)
minA=2;     % minimum cell area
minI=20;    % minimum membrane intensity
mergeT=0.8; % minimum border intensity ratio for merging cells
s2=0.5;     % second smoothing (segmentation)
maxA=3000;  % maximum cell area
boundT=0.1; % mimimal seed/boundary ratio
verbose=1;  % segmentation feedback plot

[segmentation, seeds, labels] = SegmentIm(crop,s1,minA,minI,mergeT,s2,maxA,boundT,verbose);
```

# Segmentation feedback

with `verbose=1` you can obtain a feedback visualization like this

![Segmentation feedback](ProjIm_feedback.png)
