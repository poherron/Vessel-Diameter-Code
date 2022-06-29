# Vessel-Diameter-Code

Code used for computing the diameter of blood vessels either by:
1. Drawing a line across the vessel and computing the diameter along the line 
(works for lumen label or 'wall' label where the lumen is dark like the background outside the vessel but the vessel wall is labeled with dye).
2. Fitting a circle to pixels that cross a luminance threshold (when the lumen is labeled a sobel filter is used to extract the edge of the 
thresholded mask to fit the circle only to the edge of the vessel rather than to all the pixels within the lumen).

The cross section or region mask is drawn in an earlier function using "imline" or "impoly". Apparently, Matlab no longer recommends using them
but instead says we should use "drawline", "drawpolygon" or "drawpolyline". 
