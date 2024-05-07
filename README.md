MATLAB and Python programs to take spot movies as input and try to register all frames to a common center. Versions should work similarly. Couldn't make MiniWave-capable movies using python+ffmpeg/imageio.

This is a challenging problem since normal computer vision techniques (cross correlation, phase correlation, etc.) don't seem to work well for these spot movies.

Solution used here, which is fairly successful:
1) Threshold image so pixels below some value are zero'd (~20 out of 255) 
2) Try to guess spot/lenslet spacing using mean vectors in horizontal and vertical directions, with FindPeaks (conservative params like peaks must be >50% of maximum peak). Use average spacing between peaks.
3) Assume brightest spot in average image is the center of a good spot. Form a grid, using spacing from #2, centered on this bright pixel.
4) For each frame:
5) Sum over each grid cell (lenslet), to get total intensity in each cell
6) Threshold this number (i.e. 100) to binarize a grid of suspected valid spots
7) Fit a circle to the binarized grid using LSQ. Initial params are the center of mass and r=sqrt(total/pi). In practice, these two methods give very similar results.
8) Shift (circular) each image so that the center of the fitted circle (presumed pupil) is in the center-most cell.
9) Repeat for all frames
10) Write shifted frames as output AVI file.

Other ideas/TBD:
1) Not truly registering each frame to a reference. Could try to align binarized circle images with standard techniques (xcorr).
2) A guess is that failures arise at edges of circle. Instead of binarizing, could weight cell totals at outside of circle, since cell not completely in circle.
3) Could use fancy techniques like priors since don't expect a lot of (x,y) motion or changes in size/radius.
