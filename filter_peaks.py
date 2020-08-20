import pandas as pd
import numpy as np

def generate_df_with_matches():
	match_file_path="D:/Data/MCF7_SNPs_FINAL/FINAL_SEARCH/All_spectra/OwnData_FDR1/filtered_result/pyteomics/tmp/uvpd.txt"
	data=pd.read_table(match_file_path,comment="#",names=["mz","intensity","Noise_level","Mascot","Theoretical","Precursors","Immonium_ions","Internal_fragments"],dtype={"mz":np.float64,"intensity":np.float64,
	"Noise_level":np.float64,"Mascot":np.str,"Theoretical":np.str,"Precursors":np.str,"Immonium_ions":np.str,"Internal_fragments":np.str})

	#data=all_settings["data"]
	annotation_set=set(["Theoretical","Internal_fragments","Immonium_ions","Mascot","Precursors"])

	all_matches={}
	row_count=0
	for column in annotation_set:
		column_data=data[column][data[column].notnull()]
		for index in column_data.index:
			#print index,type(index)

			row_data=data[column][index]
			for each_match in row_data.split("|"):
				match_dict={}
				ion_type,ion_index,ion_charge,ion_loss_quantity\
					,ion_neutral_loss,ion_isotope,ion_calc_mz,ion_mass_dev,ion_mass_dev_ppm\
				 		=each_match.split("/")


				match_dict["ion_type"]=ion_type
				match_dict["ion_loss_quantity"]=int(ion_loss_quantity)
				match_dict["ion_neutral_loss"]=ion_neutral_loss

				if ion_loss_quantity=="0":
					match_dict["ion_neutral_loss"]=None

				match_dict["ion_calc_mz"]=float(ion_calc_mz)
				match_dict["ion_charge"]=int(ion_charge.rstrip("+"))
				match_dict["ion_exp_mz"]=data.get_value(index,"mz")
				match_dict["ion_intensity"]=data.get_value(index,"intensity")
				match_dict["ion_source"]=column
				match_dict["ion_old_index"]=index
				match_dict["ion_mass_dev"]=float(ion_mass_dev)
				match_dict["ion_mass_dev_ppm"]=float(ion_mass_dev_ppm)
				match_dict["ion_isotope"]=ion_isotope

				


				if column=="Internal_fragments":
					match_dict["ion_sequence"],start_end_index=ion_index.split(":")
					ion_start,ion_end=start_end_index.split("-")
					match_dict["ion_start"],match_dict["ion_end"]=int(ion_start),int(ion_end)

				else:
					#match_dict["ion_sequence"]=np.nan
					match_dict["ion_start"]=0
					match_dict["ion_end"]=int(ion_index)
				
				all_matches[row_count]=match_dict

				row_count+=1
	
	filter_data=pd.DataFrame.from_dict(all_matches,orient="index")
	#filter_data[filter_data["ion_neutral_loss"]==""]=np.nan

	#writer=pd.ExcelWriter("te")
	filter_data.to_excel("text.xlsx")
	return filter_data

generate_df_with_matches()
sys.exit()



def filter_by_annotation(all_settings):

    data=all_settings["data"]

    annotation_types=["Mascot","Theoretical","Precursors","Immonium_ions","Internal_fragments"]
    annotation_set=set()

    for annotation in annotation_types:
        annotation_state=all_settings[annotation]
        if annotation_state:
            annotation_set.add(annotation)
        else:
            data[annotation]=np.nan

    #if not all_settings["match_file_path"]:
     #   annotation_set.remove("Mascot")
    filter_data=generate_df_with_matches(data,annotation_set)


    if all_settings["Precursors"]:
    	pre_fil_data=filter_data

        if all_settings["Pre_max_loss"]!="All":
        	pre_fil_data=filter_data[(filter_data.source=="Precursors") & (filter_data.ion_loss_quantity<=all_settings["Pre_max_loss"])]
        if all_settings["Pre_min_charge"]!="All":
        	pre_fil_data=pre_fil_data[filter_data.ion_charge<=all_settings["Pre_min_loss"]]
        if all_settings["Pre_remove_peaks"]:
            pre_filter_data=pre_filter_data.drop(pre_filter_data.index)


    if all_settings["Internal_fragments"]:
    	int_fil_data=filter_data

        if not all_settings["annotate_c13_peaks"]:
        	int_fil_data=filter_data[(filter_data.source=="Internal_fragments") & (filter_data.ion_isotope=="0C13")]

        if all_settings["ab_int"]:
        	int_fil_data=int_fil_data[(filter_data.source=="Internal_fragments") & (filter_data.ion_type.str() in all_settings["ab_int"])]


        else:
        	int_fil_data=int_fil_data.drop(int_fil_data.index)


        
        if all_settings["Int_nterm_P"]:
        	#int_fil_data=int_file_data[]
            data=filter_by_fragment_nterminus(data,"Internal_fragments",int(all_settings["Int_nterm_P"]))

        if all_settings["Int_cterm_DEQWH"]:
            data=filter_by_fragment_cterminus(data,"Internal_fragments",int(all_settings["Int_cterm_DEQWH"]))



        if all_settings["Int_max_len"]!="All":
            data=filter_by_fragment_length(data,"Internal_fragments",int(all_settings["Int_max_len"]))

        if all_settings["Int_max_loss"]!="All":
            data=filter_by_loss(data,"Internal_fragments",int(all_settings["Int_max_loss"]))

        if all_settings["Int_max_charge"]!="All":
            data=filter_by_charge(data,"Internal_fragments",int(all_settings["Int_max_charge"]))            


    if all_settings["Theoretical"]:
        if not all_settings["annotate_c13_peaks"]:
            data=filter_by_isotopes(data,"Theoretical")

        all_ions=all_settings["abc"]+all_settings["xyz"]
        if all_ions:
            data=filter_by_ion_type(data,"Theoretical",set(all_ions))
        else:
            data["theoretical"]=np.nan
 
        if all_settings["Theo_max_charge"]!="All":
            data=filter_by_charge(data,"Theoretical",int(all_settings["Theo_max_charge"]))

        if all_settings["Theo_max_loss"]!="All":
            data=filter_by_loss(data,"Theoretical",int(all_settings["Theo_max_loss"]))

    all_settings["data"]=data
    all_settings["annotation_set"]=annotation_set
    return all_settings





















