function outim = scaleim(im)
%This function scales the image to get good contrast, see autocontrast
%in segment_main for details.

%Einar Heiberg, inspired by Janes autocontrast implementation

%Set percentile values, see autocontast
lowerpercentile=0.02;
upperpercentile=0.99;

sortim = sort(im(:));
lowerthresh = sortim(round(lowerpercentile*length(sortim)));
upperthresh = sortim(round(upperpercentile*length(sortim)));

%Cut above/below
outim = max(im,lowerthresh);
outim = min(outim,upperthresh);

%Normalize between 0..1
outim = (outim-lowerthresh)/(upperthresh-lowerthresh);
