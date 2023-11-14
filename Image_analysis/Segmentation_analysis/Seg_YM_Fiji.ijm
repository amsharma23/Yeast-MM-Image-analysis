fovs_dir = getDirectory("/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p/FOVs/");
fovs = getFileList(fovs_dir);

for(i=0;i<fovs.length;i++)
{
fov = fovs_dir+fovs[i];
fv = fovs[i].substring(0, fovs[i].length() - 1);
fov_ts = fov+fv+"_ts/";
fov_seg = fov+fv+"_seg_im/";
fov_zip = fov+fv+"_zip/";
File.openSequence(fov_ts, "virtual");
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=1 frames=18 display=Color");
run("YeastMate", "scorethresholdsingle=0.9 scorethresholdmating=1.0 scorethresholdbudding=1.0 minnormalizationqualtile=0.015 maxnormalizationqualtile=0.985 addsinglerois=true addmatingrois=false addbuddingrois=false showsegmentation=true onlyselectedclassesinmask=false processeveryframe=true mintrackingoverlap=0.1 ipadress=127.0.0.1:11005");
selectWindow("segmentation of "+fv+"_ts");
run("Image Sequence... ", "select="+fov_seg+" dir="+fov_seg+" format=TIFF");
roiManager("Deselect");
roiManager("Save", fov_zip+fv+".zip");

close();
close();
close("ROI Manager");

}
