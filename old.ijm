selectWindow("FOV4_ts");
File.openSequence("/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p/FOVs/FOV4/FOV4_ts/", "virtual");
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=1 frames=18 display=Color");
run("YeastMate", "scorethresholdsingle=0.9 scorethresholdmating=1.0 scorethresholdbudding=1.0 minnormalizationqualtile=0.015 maxnormalizationqualtile=0.985 addsinglerois=true addmatingrois=false addbuddingrois=false showsegmentation=true onlyselectedclassesinmask=false processeveryframe=true mintrackingoverlap=0.1 ipadress=127.0.0.1:11005");
run("Image Sequence... ", "select=/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p/FOVs/FOV4/FOV4_seg_im dir=/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p/FOVs/FOV4/FOV4_seg_im/ format=TIFF");
roiManager("Deselect");
roiManager("Save", "/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p/FOVs/FOV4/FOV4_zip/RoiSet.zip");
