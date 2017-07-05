function newobjectlabel=labelunknown(objectlabel,unknownindex)

xsize=size(objectlabel,1);
ysize=size(objectlabel,2);
zsize=size(objectlabel,3);

[x,y,z]=ind2sub([xsize ysize zsize],unknownindex);
valid=(x>1 &x<xsize &y>1 &y<ysize &z>1 &z<zsize);
unknownindex=unknownindex(valid);

newobjectlabel=fastlabelunknown(uint32(objectlabel),uint32(unknownindex));

newobjectlabel=double(newobjectlabel);


