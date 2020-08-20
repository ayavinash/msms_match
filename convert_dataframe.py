import pandas as pd
import numpy as np
from pyteomics import mass,parser
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os,sys,itertools,re



def convert_dataframe():
	
	data=pd.read_table("D:/ownclowd/Projects/Spectrum_annotator_gui/scripts/tmp/uvpd_20ppm_all_loss.txt",comment="#")
	out_file_path="D:/ownclowd/Projects/Spectrum_annotator_gui/scripts/tmp/uvpd_20ppm_all_loss_new.txt"
	out_file=open(out_file_path,"w")
	out_file.write("peak_no\tion_type\tmz\tintensity\tnoise_level\tion_series\tion_index\tcharge\tneutral_loss_quantity\tneutral_losses\tisotope\tcalc_mz\tdelta_mz\tdelta_mz_ppm\n")
	peak_dict={}
	for peak_no in data.index:
		peak_dict[peak_no]=tuple([data.mz[peak_no],data.intensity[peak_no]])
	#print peak_dict
	theo_data=data.Theoretical.dropna().str.split("|")
	theo_data_index=theo_data.index
	noise=data.Noise_level[0]

	for index in peak_dict:
		if index not in theo_data_index:
			continue
		#print index,peak_dict[index],theo_data[index]
		matches=theo_data[index]
		for match in matches:
			out_file.write(str(index)+"\t"+"major_ion_series"+"\t"+str(peak_dict[index][0])+"\t"+str(peak_dict[index][1])+"\t"+\
				str(noise)+"\t"+
				match.replace("/","\t")+"\n")


convert_dataframe()