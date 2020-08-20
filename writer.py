import os,sys


def write_match_to_file(all_settings):
    print all_settings.keys()


    title=all_settings["match_file_name"]
    mz=all_settings["mz"]
    score=all_settings["score"]
    interference=all_settings["interference"]

    file_name=title
    dir_path=os.getcwd()+os.sep+"tmp"+os.sep
    print dir_path
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    match_file_path=dir_path+file_name

    match_file=open(match_file_path,"w")
    print match_file_path
    peptide=all_settings["peptide"]
    pep_charge=all_settings["pep_charge"]
    modX_seq=all_settings["pep_seq"]

    noise_level=all_settings["noise_level"]
    data=all_settings["data"]
    data["Noise_level"]=noise_level


     
    #match_file.write("####\n")
    #match_file.write("####")
    match_file.write("####"+"Title:"+title+"\n")
    match_file.write("####"+"Source:"+"NA"+"\n")
    match_file.write("####"+"Peptide:"+peptide+"\n")
    match_file.write("####"+"Modx_seq:"+modX_seq+"\n")
    match_file.write("####"+"Charge:"+str(pep_charge)+"\n")
    match_file.write("####"+"ExpMz:"+str(mz)+"\n")
    match_file.write("####"+"Score:"+"NA"+"\n")
    match_file.write("####"+"Interference:"+"NA"+"\n")
    match_file.write("####"+"Noise_type:"+"DNL"+"\n")
    match_file.write("####"+"Noise_Level:"+str(noise_level)+"\n")

    match_file.write("\tmz\tintensity\tNoise_level\tMascot\tTheoretical\tPrecursors\tImmonium_ions\tInternal_fragments\tSide_chain_losses\n")
    index_set=set(data.index)

    
    for n in data.index:
        mz,intensity,noise_level,mascot,theo,pre,imm_ion,internal,side_chain=data.mz[n],data.intensity[n],data.Noise_level[n],data.Mascot[n],data.Theoretical[n],data.Precursors[n],data.Immonium_ions[n],data.Internal_fragments[n],data.Side_chain_losses[n]
        data_line=str(n)+"\t"+str(data.mz[n])+"\t"+str(data.intensity[n])+"\t"+str(data.Noise_level[n])+"\t"+str(data.Mascot[n])+"\t"+str(data.Theoretical[n])+"\t"+str(data.Precursors[n])+"\t"+str(data.Immonium_ions[n])+"\t"+str(data.Internal_fragments[n])+"\t"+str(data.Side_chain_losses[n])+"\n"
        match_file.write(data_line)
    #match_file.write("####")    
    match_file.close()
