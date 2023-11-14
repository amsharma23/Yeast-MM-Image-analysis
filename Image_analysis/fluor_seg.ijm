fovs_dir = getDirectory("/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/FOVs/");
fovs = getFileList(fovs_dir);
print("Fov length: "+fovs.length);
for(i=22;i==22;i++)
{
	fov = fovs_dir+fovs[i];
	fv = fovs[i].substring(0, fovs[i].length() - 1);
	crops = fov+"/Crops/";
	traps = getFileList(crops);
	fov_fl_ts = fov+fv+"_ts_fl/";
	fov_ph_ts = fov+fv+"_ts_ph/";
	fov_seg = fov+fv+"_seg_im/";
	fov_zip = fov+fv+"_zip/";
	print("FoV: "+fv);
	if(i==22){ji = 39;}
	else{ji =4;}
	for(j=ji; j<traps.length-1;j++){
			trap_f = crops+traps[j];
			trp_f = fov_seg+traps[j];
			trp_z = fov_zip+traps[j]+"/";
			File.openSequence(trap_f, "virtual");
			selectWindow(traps[j].substring(0, traps[j].length() - 1));
			wait(1000);
			run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=1 frames=18 display=Color");
			wait(3000);
			selectWindow(traps[j].substring(0, traps[j].length() - 1));
			run("YeastMate", "scorethresholdsingle=0.9 scorethresholdmating=1.0 scorethresholdbudding=1.0 minnormalizationqualtile=0.015 maxnormalizationqualtile=0.985 addsinglerois=true addmatingrois=false addbuddingrois=false showsegmentation=true onlyselectedclassesinmask=false processeveryframe=true mintrackingoverlap=0.1 ipadress=127.0.0.1:11005");
			wait(1000);
			selectWindow("segmentation of "+traps[j].substring(0, traps[j].length() - 1));
			run("Image Sequence... ", "select="+trp_f+" dir="+trp_f+" format=TIFF");
			roiManager("Deselect");
			str = traps[j].substring(0, traps[j].length() - 1);
			roiManager("Save", trp_z+str+".zip");
			print("j: "+j);
			close();
			close();
			close("ROI Manager");
		}
	}
