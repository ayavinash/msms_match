# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from pyteomics import mass,parser
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os,sys,itertools,re



def internal_fragments_old(pep_list):

    trimmed_pep_list=pep_list[1:-1]
    len_trimmed_pep_list=len(trimmed_pep_list)
   
    for k in xrange(2, len_trimmed_pep_list+2):
        for i in xrange(len_trimmed_pep_list-k+1):
            end=i+k
            yield trimmed_pep_list[i:end],i+2,end+1


def check_nterm_cterm(pep_list,nterm_is_P,cterm_is_DEQWH):
    has_p,has_deqwh,ok=0,0,0

    if not (nterm_is_P and cterm_is_DEQWH):
        ok=1

    elif nterm_is_P and cterm_is_DEQWH:
        print pep_list
        if (pep_list[0][-1]=="P") and (pep_list[-1][0] in "DEQWH"):
            ok=1

    elif nterm_is_P:
        if pep_list[0][-1]=="P":
            ok=1
        else:
            ok=0
    elif cterm_is_DEQWH:
        if pep_list[-1][0] in "DEQWH":
            ok=1
        else:
            ok=0
    
    return ok



def internal_fragments(pep_list,max_len,nterm_is_P,cterm_is_DEQWH):

    trimmed_pep_list=pep_list[1:-1]
    len_trimmed_pep_list=len(trimmed_pep_list)
    if (max_len=="All") or (int(max_len)>len(trimmed_pep_list)):
        max_len=len_trimmed_pep_list
    else:
        max_len=int(max_len)
        
 
    for k in xrange(2, len_trimmed_pep_list+2):
        for i in xrange(len_trimmed_pep_list-k+1):
            end=i+k
            fragment=trimmed_pep_list[i:end]
            pos_start,pos_end=i+2,end+1
            if (pos_end-pos_start+1)<=max_len:
                ok=check_nterm_cterm(pep_list,nterm_is_P,cterm_is_DEQWH)
                if ok:
                    yield fragment,pos_start,pos_end




def annotate_internal_fragments(all_settings):
    pep_list=parser.parse(all_settings["pep_seq"],show_unmodified_termini=True,split=True)
    max_len=all_settings["Int_max_len"]

    print "Internal fragments"
    print max_len,type(max_len)
    nterm_is_P=all_settings["Int_nterm_P"]
    cterm_is_DEQWH=all_settings["Int_cterm_DEQWH"]


    fragment_list=internal_fragments(pep_list,max_len,nterm_is_P,cterm_is_DEQWH)    
    

    data=all_settings["data"]
    column="Internal_fragments"
    ab_int=all_settings["ab_int"]
    for fragment in fragment_list:
        frag_aa_comp=parser.amino_acid_composition(fragment[0],labels=aa_comp_keys)
        fragment_modx=parser.tostring(fragment[0])
        for ion_type in ab_int:
            data=generate_ion_mass(ion_type,fragment_modx+":"+str(fragment[1])+"-"+str(fragment[2]),frag_aa_comp,fragment_modx,column,all_settings)
    data.Internal_fragments[data.Internal_fragments[data.Internal_fragments==""].index]=np.nan
            
    all_settings["data"]=data

    return all_settings  

def generate_immonium_ions():
    immonium_dict={}
    aa_dict=mass.std_aa_comp
    carbonyl_comp=mass.Composition({"C":1,"O":1})
    for aa in aa_dict:
        if "-" in aa:
            continue
        immonium_comp=aa_dict[aa]-carbonyl_comp
        if aa=="C":
            immonium_mass=mass.calculate_mass(immonium_comp+aa_comp["cam"],charge=1)
        else:
            immonium_mass=mass.calculate_mass(immonium_comp,charge=1)
        immonium_dict[aa]=immonium_mass

    return immonium_dict





def initiate_globals(all_settings):
    global proton_mass
    global frag_mass_tol
    global loss_mass_dict
    global aa_comp_keys
    global aa_comp
    global immonium_dict
    global frag_mass_tol_unit

    global show_mass
    global show_sequence
    global show_mass_dev
    global show_mass_dev_ppm

    show_mass=all_settings["show_mass"]
    show_sequence=all_settings["show_sequence"]
    show_mass_dev=all_settings["show_mass_dev"]
    show_mass_dev_ppm=all_settings["show_mass_dev_ppm"]


    #proton_mass=mass.calculate_mass(mass.std__comp())
    aa_comp=get_comp_from_unimod()
    aa_comp_keys=aa_comp.keys()

    frag_mass_tol=float(all_settings["mass_tol"])
    all_settings["mass_tol"]=float(frag_mass_tol)


    immonium_dict=generate_immonium_ions()

    deviation_ppm=lambda x,y:((x-y)/y)*1000000
    global check_and_apply_mass
    global check_and_apply_pre_mass
    global annotate_fragment_isotopes


    if all_settings["mass_tol_unit"]=="ppm":
        from check_and_apply_mass_todf import check_and_apply_mass_ppm as check_and_apply_mass
        from check_and_apply_mass_todf import check_and_apply_pre_mass_ppm as check_and_apply_pre_mass
        
    else:
        from check_and_apply_mass_todf import check_and_apply_mass_dalton as check_and_apply_mass
        from check_and_apply_mass_todf import check_and_apply_pre_mass_dalton as check_and_apply_pre_mass

    if all_settings["annotate_c13_peaks"]:
        from isotopes import annotate_fragment_isotopes_all as annotate_fragment_isotopes
    else:
        from isotopes import annotate_fragment_isotopes_none as annotate_fragment_isotopes



    proton_mass=mass.calculate_mass(mass.std_aa_comp["H-"])
    ammonia_mass=mass.calculate_mass(mass.Composition({"N":1,"H":3}))
    water_mass=mass.calculate_mass(mass.Composition({"H":2,"O":1}))
    metloss_mass=mass.calculate_mass(mass.Composition({"C":1,"H":4,"S":1,"O":1}))
    all_settings["proton_mass"]=proton_mass
    std_ion_comp=mass.std_ion_comp


    std_ion_comp["a+1"]=mass.Composition({'H': -2, 'C': -1, 'O': -2})+mass.Composition({'H':1})
    std_ion_comp["x+1"]=mass.Composition({'H': -2, 'C': 1, 'O': 1})+  mass.Composition({'H':1})


    std_ion_comp["c+1"]=mass.Composition({'H': 1, 'O': -1, 'N': 1})+mass.Composition({'H': 1})
    std_ion_comp["c+2"]=mass.Composition({'H': 1, 'O': -1, 'N': 1})+mass.Composition({'H': 2})
    std_ion_comp["c-1"]=mass.Composition({'H': 1, 'O': -1, 'N': 1})+mass.Composition({'H': -1})

    std_ion_comp["z+1"]=mass.Composition({'H': -3, 'N': -1})+mass.Composition({'H': 1})
    std_ion_comp["z+2"]=mass.Composition({'H': -3, 'N': -1})+mass.Composition({'H': 2})
    std_ion_comp["z+3"]=mass.Composition({'H': -3, 'N': -1})+mass.Composition({'H': 3})



    loss_mass_dict={}
    loss_mass_dict["NH3"]=ammonia_mass
    loss_mass_dict["H2O"]=water_mass
    loss_mass_dict["CH4S"]=metloss_mass



    #return loss_mass_dict,proton_mass

def get_mascot_ann(ann):
    ann_split=ann.split("/")
    #print ann_split
    
    ion_type,ion_index,charge=ann_split[0],ann_split[1],ann_split[2]
    if charge=="1+":
        charge=""
  
    ann=(ion_type+"("+ion_index+")"+charge).replace("frag: ","").replace("ion","").replace(" ","").replace("NH3","*").replace("H2O","$").replace("-","")
    return ann
 

def assign_annotation_by_loss(all_settings):
    data=all_settings["data"]
    annotation_set=all_settings["annotation_set"]
    #data["Annotation"]=(xrange(0,len(data)))
    data["Annotation"]=""
    data["From"]=""
    mz_check_set=set()
    #data["From"]=
    annotation_order_list=["Mascot","Immonium_ions","Precursors","Theoretical","Internal_fragments","Side_chain_losses"]
    #annotation_function_list[get_mascot_ann,get_]
    annotation_set=all_settings["annotation_set"]

    assigned_index_set=set()
    for annotation in annotation_order_list:
        if annotation not in annotation_set:
            continue

        annotation_series=data[annotation][data[annotation].notnull()]
        if annotation_series.empty:
            continue
        annotation_series_index_set=set(annotation_series.index)
        not_assigned_index=annotation_series_index_set.difference(assigned_index_set)
            
        for index in not_assigned_index:
            ann=annotation_series[index]

            if annotation=="Mascot":
                final_ann=get_mascot_ann(ann)
                source="Mascot"
            elif annotation=="Immonium_ions":
                final_ann=ann.split("/")[0]
                source="Immonium_ions"
            elif annotation=="Precursors":
                best_ann=select_ann_with_least_neutral_loss(ann.split("|"))
                final_ann=get_pre_ann(best_ann)
                source="Precursors"
            elif annotation=="Theoretical":
                best_ann=select_ann_with_least_neutral_loss(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Theoretical")
                source="Theoretical"
                #new_ann=select_ann_with_least_mass_deviation(ann.split("|"))
            elif annotation=="Internal_fragments":
                best_ann=select_ann_with_least_neutral_loss(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Internal_fragments")
                source="Internal_fragments"
            elif annotation=="Side_chain_losses":
                best_ann=select_ann_with_least_neutral_loss(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Side_chain_losses")
                source="Side_chain_losses"


            else:
                final_ann=np.nan
                source=np.nan

            data.set_value(index,"Annotation",final_ann)
            data.set_value(index,"From",source)
        assigned_index_set.update(not_assigned_index)

    empty_index=data.From[data.From==""].index
#    data["From"][empty_index]=np.nan
#    data["Annotation"][empty_index]=np.nan


    for index in empty_index:
        data.set_value(index,"From",np.nan)
        data.set_value(index,"Annotation",np.nan)



    all_settings["data"]=data
    return all_settings

def assign_annotation_to_all(all_settings):
    data=all_settings["data"]
    annotation_set=all_settings["annotation_set"]
    #data["Annotation"]=(xrange(0,len(data)))
    data["Annotation"]=""
    data["From"]=""
    mz_check_set=set()
    #data["From"]=
    annotation_order_list=["Mascot","Immonium_ions","Precursors","Theoretical","Internal_fragments","Side_chain_losses"]
    #annotation_function_list[get_mascot_ann,get_]
    annotation_set=all_settings["annotation_set"]

    dvw_ions=all_settings["dvw_ions"]


    for index in data.index:
        mascot_ann=data.get_value(index,"Mascot")
        immo_ann=data.get_value(index,"Immonium_ions")
        pre_ann=data.get_value(index,"Precursors")
        theo_ann=data.get_value(index,"Theoretical")
        int_ann=data.get_value(index,"Internal_fragments")
        dvw_ann=data.get_value(index,"Side_chain_losses")

        all_anns=[]
        source=""
        if not pd.isnull(mascot_ann):
            final_ann=get_mascot_ann(mascot_ann)
            source="Mascot"
        elif not pd.isnull(immo_ann):
            final_ann=immo_ann.split("/")[0]
            source="Immonium_ions"
        else:
            
            if not pd.isnull(theo_ann):
                all_anns+=theo_ann.split("|")
            if not pd.isnull(int_ann):
                all_anns+=int_ann.split("|")
            if not pd.isnull(pre_ann):
                all_anns+=pre_ann.split("|")
            if not pd.isnull(dvw_ann):
                all_anns+=dvw_ann.split("|")

            if all_anns:
                print all_anns
                if all_settings["match_type"]==1:
                    best_ann=select_ann_with_least_neutral_loss(all_anns)
                else:
                    best_ann=select_ann_with_least_mass_deviation(all_anns)
                print best_ann

                if best_ann.startswith("M/"):
                    final_ann=get_pre_ann(best_ann)
                    source="Precursors"
                elif ":" in best_ann.split("/")[1]:
                    final_ann=extract_annotation_theo(best_ann,"Internal_fragments")
                    source="Internal_fragments"
                elif best_ann.split("/")[0] in dvw_ions:
                    final_ann=extract_annotation_theo(best_ann,"Side_chain_losses")
                    source="Side_chain_losses"

                else:
                    final_ann=extract_annotation_theo(best_ann,"Theoretical")
                    source="Theoretical"

        if (not all_anns) and (not source):
            final_ann=np.nan
            source=np.nan
 

        data.set_value(index,"Annotation",final_ann)
        data.set_value(index,"From",source)
     

    all_settings["data"]=data
    return all_settings



def assign_annotation_by_deviation(all_settings):
    data=all_settings["data"]
    annotation_set=all_settings["annotation_set"]
    #data["Annotation"]=(xrange(0,len(data)))
    data["Annotation"]=""
    data["From"]=""
    mz_check_set=set()
    #data["From"]=
    annotation_order_list=["Mascot","Immonium_ions","Precursors","Theoretical","Internal_fragments","Side_chain_losses"]

    #annotation_function_list[get_mascot_ann,get_]
    annotation_set=all_settings["annotation_set"]

    assigned_index_set=set()
    for annotation in annotation_order_list:
        if annotation not in annotation_set:
            continue

        annotation_series=data[annotation][data[annotation].notnull()]
        if annotation_series.empty:
            continue
        annotation_series_index_set=set(annotation_series.index)
        not_assigned_index=annotation_series_index_set.difference(assigned_index_set)
            
        for index in not_assigned_index:
            ann=annotation_series[index]

            if annotation=="Mascot":
                final_ann=get_mascot_ann(ann)
                source="Mascot"
            elif annotation=="Immonium_ions":
                final_ann=ann.split("/")[0]
                source="Immonium_ions"
            elif annotation=="Precursors":
                best_ann=select_ann_with_least_mass_deviation(ann.split("|"))
                final_ann=get_pre_ann(best_ann)
                source="Precursors"
            elif annotation=="Theoretical":
                best_ann=select_ann_with_least_mass_deviation(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Theoretical")
                source="Theoretical"
                #new_ann=select_ann_with_least_mass_deviation(ann.split("|"))
            elif annotation=="Internal_fragments":
                best_ann=select_ann_with_least_mass_deviation(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Internal_fragments")
                source="Internal_fragments"
            elif annotation=="Side_chain_losses":
                best_ann=select_ann_with_least_mass_deviation(ann.split("|"))
                final_ann=extract_annotation_theo(best_ann,"Side_chain_losses")
                source="Side_chain_losses"

            else:
                final_ann=np.nan
                source=np.nan

            data.set_value(index,"Annotation",final_ann)
            data.set_value(index,"From",source)
        assigned_index_set.update(not_assigned_index)

    empty_index=data.From[data.From==""].index
#    data["From"][empty_index]=np.nan
#    data["Annotation"][empty_index]=np.nan


    for index in empty_index:
        data.set_value(index,"From",np.nan)
        data.set_value(index,"Annotation",np.nan)



    all_settings["data"]=data
    return all_settings



def extract_annotation_theo_old(ann,ann_type):
    
    ann_split=ann.split("/")
    #print ann_split
    #print ann_split
    # ann
    if ann_type=="Theoretical":
        ion_type,ion_index,charge,loss_quantity,losses,isotope=ann_split[0],ann_split[1],ann_split[2],ann_split[3],ann_split[4],ann_split[5]
    else:
        ion_type,ion_index,charge,loss_quantity,losses,isotope=ann_split[0],ann_split[1].split(":")[1],ann_split[2],ann_split[3],ann_split[4],ann_split[5]
    #if len(ann_split)!=8:
     #   print ann_split
    
    
    if charge=="1+":
        charge=""

    final_ann=ion_type+"("+ion_index+")"+charge
    if loss_quantity=="0":
        return final_ann
    else:
        neutral_string=""

        if loss_quantity=="1":
            neutral_string=losses.replace("1-NH3","*").replace("1-H2O","$").replace("1-CH4S","@")
        else:
            neutral_string=losses.replace("-NH3","*").replace("-H2O","$").replace("-CH4S","@")
        if isotope=="0C13":
            isotope=""
            
        final_ann=final_ann+"("+neutral_string+")"+isotope
        return final_ann

def extract_annotation_theo(ann,ann_type):
    
    ann_split=ann.split("/")
    ion_type,ion_index,charge,loss_quantity,losses,isotope,cal_mz,mass_dev,mass_dev_ppm=ann_split[0],ann_split[1],ann_split[2],ann_split[3],ann_split[4],ann_split[5],ann_split[6],ann_split[7],ann_split[8]


    if ann_type=="Internal_fragments":
        sequence,ion_index=ion_index.split(":")
        sequence=","+sequence
    else:
        sequence=""


    if charge=="1+":
        charge_string=""
    else:
        charge_string=charge

    if loss_quantity=="0":
        neutral_string=""
    else:
        if loss_quantity=="1":
            neutral_string=losses.replace("1-NH3","*").replace("1-H2O","$").replace("1-CH4S","@")
        else:
            neutral_string=losses.replace("-NH3","*").replace("-H2O","$").replace("-CH4S","@")
            neutral_string="("+neutral_string+")"

    if isotope=="0C13":
        isotope_string=""
    else:
        isotope_string=isotope        
    #print final_ann
    final_ann=ion_type+"("+ion_index+")"+charge_string+neutral_string+isotope_string
    if show_mass:
        final_ann+=","+str(round(float(cal_mz),2))
    if show_mass_dev:
        final_ann+=","+str(round(float(mass_dev),2))
    if show_mass_dev_ppm:
        final_ann+=","+str(round(float(mass_dev_ppm),2))
    if show_sequence:
        final_ann+=sequence

    return final_ann





def select_ann_with_least_neutral_loss(all_anns):
    #print all_anns
    #(all_anns.split("|")[0].split("/"))
    min_neutral_loss=sorted(map(lambda i: int(i.split("/")[3]),all_anns))[0]
    min_loss_list=[]
    for ann in all_anns:
        #print "here",ann
        temp=ann.split("/")[3]
        if temp==str(min_neutral_loss):
            min_loss_list.append(ann)
            
            #best_ann=ann
    #print "before return",best_ann
    sep_count=sorted(map(lambda i: i.split("/")[4].count(":"),min_loss_list))[0]
   # print "#######"
    #print sep_count
    
    for ann in min_loss_list:
        sep=ann.split("/")[4].count(":")
        #print sep,sep_count
        if sep==sep_count:
            best_ann=ann
            return best_ann
            
        
    return best_ann

def select_ann_with_least_mass_deviation(all_anns):
    #print all_anns
    #(all_anns.split("|")[0].split("/"))
    if len(all_anns)==1:
        return all_anns[0]

    min_mass_deviation=sorted(map(lambda i: abs(float(i.split("/")[7])),all_anns))[0]
    #print min_mass_deviation
    #print all_anns
    min_deviation_list=[]
    for ann in all_anns:
        temp=abs(float(ann.split("/")[7]))
        if temp==min_mass_deviation:
            min_deviation_list.append(ann)
    #print min_deviation_list    
    return min_deviation_list[0]




 

def get_pre_ann(best_ann):
    #all_anns=ann.split("|")
    #all_anns=select_pre_with_no_C13(all_anns)
    #best_ann=select_ann_with_least_neutral_loss(all_anns)
    ann_split=best_ann.split("/")
    cal_mz,mass_dev,mass_dev_ppm=ann_split[6],ann_split[7],ann_split[8]
    
    #ann_split=all_anns[0].split("/")
    final_ann=ann_split[0]+"+"+ann_split[2].split("+")[0]+"H"
    neutral_loss=""
    loss_quantity,loss_string=ann_split[3],ann_split[4]
    if loss_quantity=="0":
        loss_str=""
    else:
        #print loss_string
        loss_str=loss_string.replace(",","").replace("-NH3","*").replace("-H2O","$").replace("-CH4S","@")
        #print loss_str
    if loss_string:
        final_ann=final_ann+"("+loss_str+")"
    isotope=ann_split[5]
    if isotope!="0C13":
        final_ann=final_ann+"("+isotope+")"       
    #print final_ann

    if show_mass:
        final_ann+=","+str(round(float(cal_mz),2))
    if show_mass_dev:
        final_ann+=","+str(round(float(mass_dev),2))
    if show_mass_dev_ppm:
        final_ann+=","+str(round(float(mass_dev_ppm),2))

    return final_ann



def filter_by_loss(data,column,max_loss_val):
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            loss=int(item_split[3])
            if loss>max_loss_val:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data

def filter_by_isotopes(data,column):
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            if "/0C13/" not in item:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data


def filter_by_isotope_envelope(data,column):
    ion_isotope_dict=generate_isotope_peaks(data,column)
    
    
    ion_isotope_dict=remove_if_no_monoisotopic(ion_isotope_dict)
    

    ion_isotope_dict=remove_ion_isotope_peaks(ion_isotope_dict,3)



    

    ion_isotope_dict=remove_if_broken_envelope(ion_isotope_dict,3)
    #print ion_isotope_dict


    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            ion="/".join(item_split[:5])
            if ion not in ion_isotope_dict:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data


def filter_by_isotope_envelope_consider_no_monoisotopic(data,column,all_settings):
    #strictly_monosiotopic,minimum_peaks_for_no_mono_to_consider):
    ion_isotope_dict=generate_isotope_peaks(data,column)
    
    minimum_isotope_peaks=int(all_settings["minimum_isotopes"])
    strictly_monoisotopic=all_settings["strictly_monoisotopic"]
    minimum_peaks_for_no_mono_to_consider=int(all_settings["minimum_peaks_to_consider_no_0C13"])
    print minimum_isotope_peaks,strictly_monoisotopic,minimum_peaks_for_no_mono_to_consider
    if strictly_monoisotopic:
        ion_isotope_dict=remove_if_no_monoisotopic(ion_isotope_dict)
    else:
        ion_isotope_dict=remove_if_no_monoisotopic_or_1C13(ion_isotope_dict)
    ion_isotope_dict=remove_ion_isotope_peaks(ion_isotope_dict,minimum_isotope_peaks)
    if strictly_monoisotopic:
        ion_isotope_dict=remove_if_broken_envelope(ion_isotope_dict,minimum_isotope_peaks)
    else:
        ion_isotope_dict=remove_if_broken_envelope_consider_no_monosisotopic(ion_isotope_dict,minimum_isotope_peaks,minimum_peaks_for_no_mono_to_consider)

    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            ion="/".join(item_split[:5])
            if ion not in ion_isotope_dict:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data



def generate_isotope_peaks(data,column):
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    ion_isotope_dict={}    

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            ion="/".join(item_split[:5])
            isotope=item_split[5]
            if ion not in ion_isotope_dict:
                isotope_list=[]
            else:
                isotope_list=ion_isotope_dict[ion]
            isotope_list.append(isotope)
            ion_isotope_dict[ion]=isotope_list
    return ion_isotope_dict

def remove_ion_isotope_peaks(ion_isotope_dict,min_isotopes):
    edited_ion_isotope_dict=dict(ion_isotope_dict)
    for ion in ion_isotope_dict:
        if len(ion_isotope_dict[ion])<min_isotopes:
            del edited_ion_isotope_dict[ion]
    return edited_ion_isotope_dict

def remove_if_no_monoisotopic(ion_isotope_dict):
    edited_ion_isotope_dict=dict(ion_isotope_dict)
    for ion in ion_isotope_dict:
        if ion_isotope_dict[ion][0]!="0C13":
            del edited_ion_isotope_dict[ion]
    return edited_ion_isotope_dict

def remove_if_no_monoisotopic_or_1C13(ion_isotope_dict):
    edited_ion_isotope_dict=dict(ion_isotope_dict)
    allowed_isotopes=tuple(["0C13","1C13"])
    for ion in ion_isotope_dict:
        if ion_isotope_dict[ion][0] not in allowed_isotopes:
            del edited_ion_isotope_dict[ion]
    return edited_ion_isotope_dict



def remove_if_broken_envelope(ion_isotope_dict,min_envelope_to_check):
    edited_ion_isotope_dict=dict(ion_isotope_dict)
    isotope_series_to_check=tuple(range(min_envelope_to_check))
    
    for ion in ion_isotope_dict:
        isotope_list=ion_isotope_dict[ion]
        for isotope_number in isotope_series_to_check:
            current_number=int(isotope_list[isotope_number].split("C13")[0])
            if current_number!=isotope_number:
                del edited_ion_isotope_dict[ion]
                break
    return edited_ion_isotope_dict

def remove_if_broken_envelope_consider_no_monosisotopic(ion_isotope_dict,min_envelope_to_check,allowed_peaks):
    edited_ion_isotope_dict=dict(ion_isotope_dict)
    isotope_series_to_check=tuple(range(min_envelope_to_check))
    isotope_series_to_check_without_mono=[]
    if len(isotope_series_to_check)>=allowed_peaks:
        isotope_series_to_check_without_mono=isotope_series_to_check[1:]
    if isotope_series_to_check_without_mono:
        isotope_series_to_check=isotope_series_to_check_without_mono

    #print ion_isotope_dict
    for ion in ion_isotope_dict:
        isotope_list=ion_isotope_dict[ion]
        for isotope_number in isotope_series_to_check:
            current_number=int(isotope_list[isotope_number].split("C13")[0])

            if current_number!=isotope_number:
                del edited_ion_isotope_dict[ion]
                break
    #print edited_ion_isotope_dict
    return edited_ion_isotope_dict



def filter_by_charge(data,column,max_charge_val):
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            charge=int(item_split[2].rstrip("+"))
            if charge>max_charge_val:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data

def filter_by_charge_ifmz_is_less(data,column,max_charge_val,min_mz_level):
    col_view=data[column][(data[column].notnull()) & (data.mz<min_mz_level)].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            charge=int(item_split[2].rstrip("+"))
            if charge>max_charge_val:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data


def filter_by_charge_precursor(data,column,min_charge_val):
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
    

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            charge=int(item_split[2].rstrip("+"))
            if charge<min_charge_val:
                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    return data

def remove_peaks(data,annotation):

    print len(data),annotation
    #print data[annotation].notnull()
    data=data.drop(data[annotation][data[annotation].notnull()].index)
    print len(data)
    return data


 
def filter_by_fragment_length(data,column,max_int_len):
    print len(data.Internal_fragments[data.Internal_fragments.notnull()])
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            item_split=item.split("/")
            start,end=item_split[1].split(":")[1].split("-")
            length=int(end)-int(start)+1
            if length>max_int_len:
               list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    print len(data.Internal_fragments[data.Internal_fragments.notnull()])
    return data   

def filter_by_fragment_nterminus(data,column,amino_acids):
    print amino_acids,len(data[column][data[column].notnull()])
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
     
    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            nterm_aa=item.split("/")[1][0]
            if nterm_aa not in "P":

                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    print len(data[column][data[column].notnull()])
    return data



def filter_by_fragment_cterminus(data,column,amino_acids):
    print amino_acids,len(data[column][data[column].notnull()])
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)
     
    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            cterm_aa=item.split("/")[1].split(":")[0][-1]

            if cterm_aa not in "DEQWH":

                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
    print len(data[column][data[column].notnull()])
    return data   



def filter_by_ion_type(data,column,ion_types):
   
    col_view=data[column][data[column].notnull()].str.split("|")
    col_view_copy=pd.Series.copy(col_view)

    for index,item_list in col_view.iteritems():
        list_to_modify=list(col_view_copy[index])
        for item in item_list:
            ion_type=item.split("/")[0]
            if ion_type not in ion_types:

                list_to_modify.remove(item)

        if not list_to_modify:
            new_ann=np.nan
        else:
            new_ann="|".join(list_to_modify)
        data.set_value(index,column,new_ann)
    
   
    return data   





def spectrum_plotter_old(all_settings):
    print all_settings.keys()
    print "-----------"

    data=all_settings["data"]
    modx_seq=all_settings["pep_seq"]
    title=all_settings["title"]
    mz=all_settings["mz"]
    charge=all_settings["pep_charge"]
    score=all_settings["score"]
    interference=all_settings["interference"]
    noise=all_settings["noise_level"]
    noise=float(noise)
    print type(noise)
    #print data
    within=float(all_settings["within"])
    theo_color=all_settings["theo_color"]
    pre_color=all_settings["pre_color"]
    int_color=all_settings["int_color"]
    immo_color=all_settings["immo_color"]
    mascot_color=all_settings["mascot_color"]
    original_data=all_settings["original_data"]
    print data.mz


    fig,ax=plt.subplots(figsize=(8,4),dpi=1000)
    ax2=ax.twinx()

    #fig,ax=plt.subplots(dpi=100)
    #ax.axis([min(data.mz),max(data.mz),min(data.intensity),max(data.intensity)])
    base_peak_intensity=max(original_data.intensity)
    print base_peak_intensity
    intensity_percent=(original_data.intensity/base_peak_intensity)*100.0
    print intensity_percent


    plt.vlines(data.mz,[0],data.intensity,"grey")
    plt.hlines(noise,min(data.mz),max(data.mz),color="blue")
    plt.title("Peptide: "+modx_seq+", "+str(charge)+"+, m/z:"+str(round(float(mz),2))+", Score:"+score+"\n"+title+", Interference:"+str(round(float(interference),2))+"%, Noise:"+str(round(noise,2))+"\n",fontsize=20,fontweight="bold")
    plt.xlabel("m/z",fontsize=15,fontweight="bold")
    plt.ylabel("intensity",fontsize=15,fontweight="bold")
    plt.tight_layout()
    data=data[data.intensity>noise]
    plt.vlines(data.mz[data.From=="Precursors"],[0],data.intensity[data.From=="Precursors"],pre_color)
    plt.vlines(data.mz[data.From=="Internal_fragments"],[0],data.intensity[data.From=="Internal_fragments"],int_color)
    plt.vlines(data.mz[data.From=="Immonium_ions"],[0],data.intensity[data.From=="Immonium_ions"],immo_color)
    plt.vlines(data.mz[data.From=="Theoretical"],[0],data.intensity[data.From=="Theoretical"],theo_color)

    plt.vlines(data.mz[data.From=="Mascot"],[0],data.intensity[data.From=="Mascot"],"red")
    #data=data.sort_values(by="intensity",ascending=False)
    data=data.sort("intensity",ascending=False)
    r=fig.canvas.get_renderer()
    text_list=[]
    col_dict={"Mascot":mascot_color,"Immonium_ions":immo_color,"Theoretical":theo_color,"Internal_fragments":int_color,"Precursors":pre_color}

    for index in data.From[(data.From.notnull()) & (data.From!="")].index:
        source=data.get_value(index,"From")
        color=col_dict[source]
        text_obj=plt.text(data.get_value(index,"mz"),data.get_value(index,"intensity"),data.get_value(index,"Annotation"),rotation=90,color=color,horizontalalignment="center",verticalalignment="bottom",size="smaller",visible=False,fontsize=10)
        text_data=[text_obj,data.get_value(index,"mz"),data.get_value(index,"intensity")]
        text_list.append(text_data)
    if not text_list:
        plt.tight_layout()
        plt.show()
    else:        
        text_list[0][0].set_visible(True)
        visible_list=[]
        visible_list.append(text_list[0])
        for n,each_text in enumerate(text_list[1:]):
            for k in range(0,len(visible_list)):
                if abs(each_text[1]-visible_list[k][1])<within:
                    too_close=1
                    break
                else:
                    too_close=0    
            if too_close:
                continue
            each_text[0].set_visible(True)
            visible_list.append(each_text)
    plt.axis([min(original_data.mz)-20.0,max(data.mz)+20.0,min(original_data.intensity),max(original_data.intensity)+0.1*max(original_data.intensity)],labelsize=10,size=10)
    
    #ax.plot(original_data.mz,original_data.intensity)



    plt.tight_layout()
    print "---------"
    plt.show()



def spectrum_plotter(all_settings):
    print all_settings.keys()
    print "-----------"

    data=all_settings["data"]
    modx_seq=all_settings["pep_seq"]
    title=all_settings["title"]
    mz=all_settings["mz"]
    charge=all_settings["pep_charge"]
    score=all_settings["score"]
    interference=all_settings["interference"]
    noise=all_settings["noise_level"]
    noise=float(noise)
    print type(noise)
    print mz
    print type(mz)

    
    #print data
    within=float(all_settings["within"])
    theo_color=all_settings["theo_color"]
    pre_color=all_settings["pre_color"]
    int_color=all_settings["int_color"]
    immo_color=all_settings["immo_color"]
    mascot_color=all_settings["mascot_color"]
    original_data=all_settings["original_data"]
    count_basic_residues=all_settings["count_basic_residues"]
    dvw_color="brown"

    #print data.mz
    
    fig,ax=plt.subplots(figsize=(14,7),dpi=100)
    ax2=ax.twinx()

    #fig,ax=plt.subplots(dpi=100)
    #ax.axis([min(data.mz),max(data.mz),min(data.intensity),max(data.intensity)])
    base_peak_intensity=max(original_data.intensity)
    ax.vlines(data.mz,[0],data.intensity,"grey")
    #if noise >base_peak_intensity:
    ax.hlines(noise,min(data.mz),max(data.mz),color="blue")
    try:
        ax.set_title("Peptide: "+modx_seq+", "+str(charge)+"+, m/z:"+str(round(float(mz),2))+", Score:"+score+"\n"+title+", Interference:"+str(round(float(interference),2))+"%, Noise level:"+str(round(noise,2))+"\n",fontsize=15)
    except:
        pass
   
    ax.set_xlabel("m/z",fontsize=15)
    ax.set_ylabel("ion current",fontsize=15)
    ax2.set_ylabel("% of base peak",fontsize=15)

    plt.tight_layout()
    if noise!=base_peak_intensity and len(data.intensity[data.intensity>noise])>=8:
        data=data[data.intensity>noise]

    #len(data.intensity[data.intensity>noise])
    

    ax.vlines(data.mz[data.From=="Precursors"],[0],data.intensity[data.From=="Precursors"],pre_color)
    ax.vlines(data.mz[data.From=="Internal_fragments"],[0],data.intensity[data.From=="Internal_fragments"],int_color)
    ax.vlines(data.mz[data.From=="Immonium_ions"],[0],data.intensity[data.From=="Immonium_ions"],immo_color)
    ax.vlines(data.mz[data.From=="Theoretical"],[0],data.intensity[data.From=="Theoretical"],theo_color)
    ax.vlines(data.mz[data.From=="Side_chain_losses"],[0],data.intensity[data.From=="Side_chain_losses"],"brown")




           
    ax.vlines(data.mz[data.From=="Mascot"],[0],data.intensity[data.From=="Mascot"],"red")
    data=data.sort_values(by="intensity",ascending=False)
    #data=data.sort("intensity",ascending=False)
    r=fig.canvas.get_renderer()
    text_list=[]
    col_dict={"Mascot":mascot_color,"Immonium_ions":immo_color,"Theoretical":theo_color,"Internal_fragments":int_color,"Precursors":pre_color,"Side_chain_losses":"brown"}

    for index in data.From[(data.From.notnull()) & (data.From!="")].index:
        source=data.get_value(index,"From")
        color=col_dict[source]
        ann_val=data.get_value(index,"Annotation")
        #print ann_val
        ann_val=re.sub(r"\$","**",ann_val)
        ann_val=re.sub(r"\+","$\mathregular{^{+}}$",ann_val)
        #print ann_val
        #ann_val=re.sub(("+","$\mathregular{^{+}}$")


        text_obj=ax.text(data.get_value(index,"mz"),data.get_value(index,"intensity"),ann_val,rotation=90,color=color,horizontalalignment="center",verticalalignment="bottom",size="large",visible=False)
        text_data=[text_obj,data.get_value(index,"mz"),data.get_value(index,"intensity")]
        text_list.append(text_data)
    
    ax.set_xlim(min(original_data.mz),max(original_data.mz))
    ax.set_ylim(0,base_peak_intensity)
    if noise==base_peak_intensity:
        ax.set_ylim(0,base_peak_intensity+0.1*base_peak_intensity)
       
    ax2.set_ylim(0,100)
    ax2.set_yticks([20,40,60,80,100])
    if text_list:    

        transf = ax.transData.inverted()
        text_list[0][0].set_visible(True)
        visible_list=[]
        bb=text_list[0][0].get_window_extent(renderer=r)
        text_loc=text_list[0][2]
        new_bb=bb.transformed(transf)
        width=new_bb.width
        height=new_bb.height        
        print "Width,height",width,height
        visible_list.append(text_list[0])
        ylim_up=max(original_data.intensity)+0.1*max(original_data.intensity)
        ylim_up_percent=110.0
        actual_height=text_loc+height
        print actual_height

        if actual_height>ylim_up:
            ylim_up_percent=ylim_up_percent+(((actual_height-ylim_up)/ylim_up)*100.0)
            ylim_up=actual_height

        print ylim_up,ylim_up_percent
        ax.set_ylim(0,ylim_up)
        max_point=int(ylim_up)
        point=max_point/5
        print range(point,max_point+point,point)

        #ax.set_yticks(range(point,max_point,point))

        ax2.set_ylim(0,ylim_up_percent)
        
        ax.set_xlim(min(original_data.mz)-2*width,max(original_data.mz)+width)
        #ax.set_xlim(min(original_data.mz)-2*width,max(data.mz)+2*width)

        #ax.hlines(max(data.intensity),min(data.mz),max(data.mz),color="magenta")
        for n,each_text in enumerate(text_list[1:]):
            for k in range(0,len(visible_list)):
                #print abs(each_text[1]-visible_list[k][1]),width
                if abs(each_text[1]-visible_list[k][1])<within:
                    too_close=1
                    break
                else:
                    too_close=0    
            if too_close:
                continue
            each_text[0].set_visible(True)
            visible_list.append(each_text)
    fontsize=15
    fontweight="bold"
    tick_list=[ax.xaxis.get_majorticklabels(),ax.yaxis.get_majorticklabels(),ax2.yaxis.get_majorticklabels()]
    for each_type in tick_list:
        for tick in each_type:
            tick.set_fontsize(fontsize)
            #tick.set_fontweight(fontweight)

    plt.tight_layout()
    print "---------"
    plt.show()




def filter_by_annotation(all_settings):

    data=all_settings["data"]

    annotation_types=["Mascot","Theoretical","Precursors","Immonium_ions","Internal_fragments","Side_chain_losses"]
    annotation_set=set()

    for annotation in annotation_types:
        annotation_state=all_settings[annotation]
        if annotation_state:
            annotation_set.add(annotation)
        else:
            data[annotation]=np.nan

    #if not all_settings["match_file_path"]:
     #   annotation_set.remove("Mascot")


    if all_settings["Precursors"]:

        if all_settings["Pre_max_loss"]!="All":
            data=filter_by_loss(data,"Precursors",int(all_settings["Pre_max_loss"]))

        if all_settings["Pre_min_charge"]!="All":
            data=filter_by_charge_precursor(data,"Precursors",int(all_settings["Pre_min_charge"]))
        if all_settings["Pre_remove_peaks"]:
            data=remove_peaks(data,"Precursors")




    if all_settings["Internal_fragments"]:
        if not all_settings["annotate_c13_peaks"]:
            data=filter_by_isotopes(data,"Internal_fragments")
        else:
            #data=filter_by_isotope_envelope(data,"Internal_fragments")
            data=filter_by_isotope_envelope_consider_no_monoisotopic(data,"Internal_fragments",all_settings)



        if all_settings["ab_int"]:
            data=filter_by_ion_type(data,"Internal_fragments",set(all_settings["ab_int"]))
        else:
            data["Internal_fragments"]=np.nan

        
        if all_settings["Int_nterm_P"]:
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
        else:
            #data=filter_by_isotope_envelope(data,"Theoretical")
            data=filter_by_isotope_envelope_consider_no_monoisotopic(data,"Theoretical",all_settings)




        all_ions=all_settings["abc"]+all_settings["xyz"]
        if all_ions:
            data=filter_by_ion_type(data,"Theoretical",set(all_ions))
        else:
            data["theoretical"]=np.nan
 
        if all_settings["Theo_max_charge"]!="All":
            data=filter_by_charge(data,"Theoretical",int(all_settings["Theo_max_charge"]))
        if not all_settings["count_basic_residues"]:
            data=filter_by_charge_ifmz_is_less(data,"Theoretical",2,500.0)
            


        if all_settings["Theo_max_loss"]!="All":
            data=filter_by_loss(data,"Theoretical",int(all_settings["Theo_max_loss"]))


    if all_settings["Side_chain_losses"]:
        if not all_settings["annotate_c13_peaks"]:
            data=filter_by_isotopes(data,"Side_chain_losses")
        else:
            data=filter_by_isotope_envelope_consider_no_monoisotopic(data,"Side_chain_losses",all_settings)

        dvw_ions=all_settings["dvw_ions"]    
        if dvw_ions:
            data=filter_by_ion_type(data,"Side_chain_losses",set(dvw_ions))
        else:
            data["Side_chain_losses"]=np.nan

    all_settings["data"]=data
    all_settings["annotation_set"]=annotation_set
    return all_settings


def annotate_mascot(data,ions):
    
    for ion in ions:
        ion_name,charge,index=ion["name"],ion["charge"],ion["index"]
        frag_array=ion["FragmentArray"]
        for array in frag_array:
            #print "-------"
            #print array
            if array["measure_ref"]=="m_mz":
                mascot_mz_array=array["values"]
                for k,mascot_mz in enumerate(mascot_mz_array):
                    i,=np.where(data["mz"]==mascot_mz)
                    #if charge>1:
                    peak=ion_name+"/"+str(index[k])+"/"+str(charge)+"+"+"/"
                    #peak=ion_name+"("+str(index[k])+")"+str(charge)+"+"
                   
                    #else:
                     #   peak=ion_name+"("+str(index[k])+")"
                        
                        
                    ann=data.get_value(i[0],"Mascot")
                    if ann.startswith(""):
                        data.set_value(i[0],"Mascot",peak+"0C13")
                    else:
                        data.set_value(i[0],"Mascot",ann+"|"+peak+"0C13")
                    neutral_mass=charge*(mascot_mz-proton_mass)
                    
                    data=annotate_fragment_isotopes(data,charge,peak,neutral_mass,"Mascot")
                    
    return data





def get_isotope_peaks(pep_exp_mr):
    isotope_peaks=3
    pep_exp_mr=float(pep_exp_mr)
    if pep_exp_mr<1000.0:
        isotope_peaks=3
    elif 1000.0<=pep_exp_mr<=1999.99:
        isotope_peaks=4
    elif 2000.0<=pep_exp_mr<=2999.99:
        isotope_peaks=5
    elif 3000.0<=pep_exp_mr<=3999.99:
        isotope_peaks=6
    elif 4000.0<=pep_exp_mr<=4999.99:
        isotope_peaks=7
    elif 5000.0<=pep_exp_mr<=5999.99:
        isotope_peaks=8
    elif 6000.0<=pep_exp_mr<=6999.99:
        isotope_peaks=9
    elif pep_exp_mr>=7000.0:
        isotope_peaks=10
    return isotope_peaks    
            



   
def apply_isotope_precursors(isotope_peaks,pep_mono_mass,delta,min_mz,max_mz,data,charge,loss_string,ion_index):
    for isotopic_peak in range(1,isotope_peaks+1):
        isotope_mass=pep_mono_mass+isotopic_peak*delta
        if (isotope_mass > min_mz) and (isotope_mass < max_mz):
            #isotope_peak="M-"+str(isotopic_peak)+"-C13/"+str(charge)+"+/"+loss_string+"/"+str(isotope_mass)
            isotope_peak="M/"+ion_index+"/"+str(charge)+loss_string+str(isotopic_peak)+"C13/"+str(isotope_mass)
            data=check_and_apply_pre_mass(data,isotope_mass,isotope_peak,frag_mass_tol)
    return data
    

    
def apply_all_losses_pre(data,neutral_mass,all_losses,charge,base,column,min_mz,max_mz,isotope_peaks,delta,ion_index):
    for each_set in all_losses:
        ion_mass=neutral_mass
        total_loss,loss_string,loss_string_list=0,"",[]
        for each_loss in each_set:
            loss_quantity,loss_type=each_loss[0],each_loss[1]
            ion_mass-=loss_quantity*loss_mass_dict[loss_type]
            total_loss+=loss_quantity
            loss_string_list.append(str(loss_quantity)+"-"+loss_type)
            
        loss_mono_mass=(ion_mass+proton_mass*charge)/charge
        loss_string=":".join(loss_string_list)
        #loss_peak="+/"+str(total_loss)+"/"+loss_string+
        if (loss_mono_mass > min_mz) and (loss_mono_mass < max_mz):
            #loss_peak=base+"/"+str(charge)+"+/"+loss_string+"/"+str(loss_mono_mass)
            
            loss_peak=base+str(total_loss)+"/"+loss_string+"/0C13/"+str(loss_mono_mass)
            data=check_and_apply_pre_mass(data,loss_mono_mass,loss_peak,frag_mass_tol)
        loss_string="+/"+str(total_loss)+"/"+loss_string+"/"   
        data=apply_isotope_precursors(isotope_peaks,loss_mono_mass,delta,min_mz,max_mz,data,charge,loss_string,ion_index)   

    return data            




def annotate_precorsor(all_settings):
    column="Precursors"
    data=all_settings["data"]
    pep_charge=all_settings["pep_charge"]
    modX_seq=all_settings["pep_seq"]
    pep_list=parser.parse(all_settings["pep_seq"],show_unmodified_termini=True,split=True)
    max_mz=all_settings["max_mz"]
    min_mz=all_settings["min_mz"]


    
    pep_comp=parser.amino_acid_composition(modX_seq,aa_comp=aa_comp)
    pep_neutral_mass=mass.fast_mass2(modX_seq,aa_comp=aa_comp)
    isotope_peaks=get_isotope_peaks(pep_neutral_mass)
    max_global_ammonia_loss=all_settings["max_global_ammonia_loss"]
    max_global_water_loss=all_settings["max_global_water_loss"]


    #print "pep_charge",pep_charge
    
    all_losses=get_multiple_neutral_losses(pep_comp,max_global_ammonia_loss,max_global_water_loss)

    if all_losses:
        loss_charges=range(1,pep_charge+1)
               
    else:
        loss_charges=None
    
    #ion_index=str(len(peptide))
    ion_index=str(len(pep_list))

    for charge in range(1,pep_charge+1):
        delta=1.0/charge
        pep_mono_mass=(pep_neutral_mass+charge*proton_mass)/charge
        if (pep_mono_mass > min_mz) and (pep_mono_mass < max_mz):
            #generate_peak()
            mono_peak="M/"+str(ion_index)+"/"+str(charge)+"+/0//0C13/"+str(pep_mono_mass)
            data=check_and_apply_pre_mass(data,pep_mono_mass,mono_peak,frag_mass_tol)
            
        data=apply_isotope_precursors(isotope_peaks,pep_mono_mass,delta,min_mz,max_mz,data,charge,"+/0//",ion_index)    

    if loss_charges:
        for charge in loss_charges:
            delta=1.0/charge
            base="M/"+ion_index+"/"+str(charge)+"+/"
            apply_all_losses_pre(data,pep_neutral_mass,all_losses,charge,base,column,min_mz,max_mz,isotope_peaks,delta,ion_index)
    
    empty_index=data.Precursors[data.Precursors==""].index
    for index in empty_index:
        data.set_value(index,"Precursors",np.nan)      
    data["Precursors"][empty_index]=np.nan

    all_settings["data"]=data        
    return all_settings




def get_peptide_from_pep_seq(pep_seq):
    only_peptide=pep_seq.split("-")[1]
    pep_list=[]
    for aa in only_peptide:
        if aa in "abcdefghijklmnopqrstuvwxyz":
            continue
        pep_list.append(aa)
    
    return "".join(pep_list)



def annotate_immonium_ions(all_settings):
    data=all_settings["data"]
    peptide=all_settings["peptide"]


    immo_in_pep={}
    available_aa=set(immonium_dict.keys()).intersection(set(list(peptide)))
    for aa in available_aa:
        immo_in_pep[aa]=immonium_dict[aa]
        
    for immo_aa in immo_in_pep:
        immo_mass=immo_in_pep[immo_aa]

        data=check_and_apply_mass(data,immo_mass,immo_aa,"Immonium_ions",frag_mass_tol)
    empty_index=data.Immonium_ions[data.Immonium_ions==""].index
    for index in empty_index:
       data.set_value(index,"Immonium_ions",np.nan)
    all_settings["data"]=data
    return all_settings


def annotate_signal_noise_paper(data):
    #print mass.std_aa_mass

    #print proton_mass
    #sys.exit()

    data_index=data.index
    sorted_data=data.copy(deep=True)
    #print sorted_data
    #sorted_data=sorted_data.sort_values("intensity",ascending=True)
    sorted_data=sorted_data.sort_values("intensity",ascending=True)
    noise=sorted_data.intensity.iloc[0]
    noise=sorted_data.intensity.iloc[0]
    second_int=sorted_data.intensity.iloc[1]
    predicted_noise_for_second=1.5*noise
    snr_min=1.5
    len_sorted_data=len(sorted_data)
    k_list=range(1,len(sorted_data)+1)
    print predicted_noise_for_second
    snr_est=second_int/float(predicted_noise_for_second)
    if snr_est>snr_min:
        noise_level=predicted_noise_for_second
        print "Second peak is signal"
    else:
        noise_level=sorted_data.intensity.iloc[1]
        for k in k_list[2:-1]:
            coeff=np.polyfit(k_list[:k],sorted_data.intensity.iloc[:k],1)
            predicted_noise_for_k=coeff[0]*k_list[k]+coeff[1]
            snr_est=sorted_data.intensity.iloc[k]/float(predicted_noise_for_k)
            if snr_est>snr_min:
                noise_level=predicted_noise_for_k
                print "Signal--",sorted_data.intensity.iloc[k],"Noise--",noise_level,"found at k--",k,"total_peaks--",len_sorted_data
                break
            else:
                noise_level=sorted_data.intensity.iloc[k]    
    if not noise_level:
        noise_level=sorted_data.intensity.iloc[k]
    data["Noise_level"]=noise_level
    return data,noise_level

#annotate_signal_noise_paper()



#def annotate_theoretical(data,peptide,pep_charge,modX_seq,pep_list,column,min_mz,max_mz,abc,xyz,max_global_water_loss,max_global_ammonia_loss,proton_mass):


def annotate_theoretical(all_settings):
    #loss_mass_dict,proton_mass=some_masses()
    

    pep_list=parser.parse(all_settings["pep_seq"],show_unmodified_termini=True,split=True)
    abc=all_settings["abc"]
    xyz=all_settings["xyz"]
    pep_charge=int(all_settings["pep_charge"])
    #max_global_water_loss=all_settings["Theo_max_loss"]
    #max_global_ammonia_loss=all_settings["Theo_max_loss"]
    #all_settings["proton_mass"]=proton_mass
    #all_settings["loss_mass_dict"]=loss_mass_dict
    column="Theoretical"
   
#    if pep_charge>3:
 #       pep_charge=pep_charge-1

 #   all_settings["pep_charge"]=pep_charge





    for i in xrange(1,len(pep_list)):
    #for i in xrange(1,len(peptide)):
        abc_fragment=pep_list[:i]
        abc_ion_index=i
        abc_frag_comp=parser.amino_acid_composition(abc_fragment,labels=aa_comp_keys)
     
        xyz_fragment=pep_list[i:]
        xyz_ion_index=len(pep_list)-i
        xyz_frag_comp=parser.amino_acid_composition(xyz_fragment,labels=aa_comp_keys)

        for ion_type in xyz:
            data=generate_ion_mass(ion_type,str(xyz_ion_index),xyz_frag_comp,xyz_fragment,column,all_settings)
        for ion_type in abc:
            data=generate_ion_mass(ion_type,str(abc_ion_index),abc_frag_comp,abc_fragment,column,all_settings)




    #empty_inde        
    data["Theoretical"][data["Theoretical"][data["Theoretical"]==""].index]=np.nan
    
    data["Mascot"][data["Mascot"][data.Mascot==""].index]=np.nan
    all_settings["data"]=data

    return all_settings


def create_side_chain_loss_mass_dict():
    da_dict={}
    da_dict["S"]=mass.calculate_mass(mass.Composition({"O":1}))
    da_dict["V"]=mass.calculate_mass(mass.Composition({"C":1,"H":2}))
    da_dict["T"]=mass.calculate_mass(mass.Composition({"O":1}))
    da_dict["C"]=mass.calculate_mass(mass.Composition({"S":1}))
    da_dict["I"]=mass.calculate_mass(mass.Composition({"C":2,"H":4}))
    da_dict["L"]=mass.calculate_mass(mass.Composition({"C":2,"H":6}))
    da_dict["N"]=mass.calculate_mass(mass.Composition({"C":1,"H":1,"N":1,"O":1}))
    da_dict["D"]=mass.calculate_mass(mass.Composition({"C":1,"O":2}))
    da_dict["Q"]=mass.calculate_mass(mass.Composition({"C":2,"H":3,"N":1,"O":1}))
    da_dict["K"]=mass.calculate_mass(mass.Composition({"C":3,"H":7,"N":1}))
    da_dict["E"]=mass.calculate_mass(mass.Composition({"C":2,"H":2,"O":2}))
    da_dict["M"]=mass.calculate_mass(mass.Composition({"C":2,"H":4,"S":1}))
    da_dict["F"]=mass.calculate_mass(mass.Composition({"C":6,"H":4}))
    da_dict["R"]=mass.calculate_mass(mass.Composition({"C":3,"H":7,"N":3}))

    db_dict={}
    db_dict["T"]=mass.calculate_mass(mass.Composition({"C":1,"H":2}))
    db_dict["I"]=mass.calculate_mass(mass.Composition({"C":1,"H":2}))

    v_dict={}
    v_dict["A"]=mass.calculate_mass(mass.Composition({"C":1,"H":4}))
    v_dict["S"]=mass.calculate_mass(mass.Composition({"C":1,"H":4,"O":1}))
    v_dict["V"]=mass.calculate_mass(mass.Composition({"C":3,"H":8}))
    v_dict["C"]=mass.calculate_mass(mass.Composition({"C":1,"H":4,"S":1}))
    v_dict["I"]=mass.calculate_mass(mass.Composition({"C":4,"H":10}))
    v_dict["L"]=mass.calculate_mass(mass.Composition({"C":4,"H":10}))
    v_dict["N"]=mass.calculate_mass(mass.Composition({"C":2,"H":5,"N":1,"O":1}))
    v_dict["D"]=mass.calculate_mass(mass.Composition({"C":2,"H":4,"O":2}))
    v_dict["Q"]=mass.calculate_mass(mass.Composition({"C":3,"H":7,"N":1,"O":1}))
    v_dict["K"]=mass.calculate_mass(mass.Composition({"C":4,"H":11,"N":1}))
    v_dict["E"]=mass.calculate_mass(mass.Composition({"C":3,"H":6,"O":2}))
    v_dict["M"]=mass.calculate_mass(mass.Composition({"C":3,"H":8,"S":1}))
    v_dict["H"]=mass.calculate_mass(mass.Composition({"C":4,"H":6,"N":2}))
    v_dict["F"]=mass.calculate_mass(mass.Composition({"C":7,"H":8}))
    v_dict["R"]=mass.calculate_mass(mass.Composition({"C":4,"H":11,"N":3}))
    v_dict["Y"]=mass.calculate_mass(mass.Composition({"C":7,"H":8,"O":1}))
    v_dict["W"]=mass.calculate_mass(mass.Composition({"C":9,"H":9,"N":1}))

    wa_dict={}
    wa_dict["S"]=mass.calculate_mass(mass.Composition({"O":1,"H":1}))
    wa_dict["V"]=mass.calculate_mass(mass.Composition({"C":1,"H":3}))
    wa_dict["T"]=mass.calculate_mass(mass.Composition({"O":1,"H":1}))
    wa_dict["C"]=mass.calculate_mass(mass.Composition({"S":1,"H":1}))
    wa_dict["I"]=mass.calculate_mass(mass.Composition({"C":2,"H":5}))
    wa_dict["L"]=mass.calculate_mass(mass.Composition({"C":3,"H":7}))
    wa_dict["N"]=mass.calculate_mass(mass.Composition({"C":1,"H":2,"N":1,"O":1}))
    wa_dict["D"]=mass.calculate_mass(mass.Composition({"C":1,"H":1,"O":2}))
    wa_dict["Q"]=mass.calculate_mass(mass.Composition({"C":2,"H":4,"N":1,"O":1}))
    wa_dict["K"]=mass.calculate_mass(mass.Composition({"C":3,"H":8,"N":1}))
    wa_dict["E"]=mass.calculate_mass(mass.Composition({"C":2,"H":3,"O":2}))
    wa_dict["M"]=mass.calculate_mass(mass.Composition({"C":2,"H":5,"S":1}))
    wa_dict["R"]=mass.calculate_mass(mass.Composition({"C":3,"H":8,"N":3}))

    wb_dict={}
    wb_dict["T"]=mass.calculate_mass(mass.Composition({"C":1,"H":3}))
    wb_dict["I"]=mass.calculate_mass(mass.Composition({"C":1,"H":3}))

    aa_string="ACDEFGHIKLMNPQRSTVWY"
    dict_list=[da_dict,db_dict,v_dict,wa_dict,wb_dict]
    for each_dict in dict_list:
        for aa in aa_string:
            if aa not in each_dict:
                each_dict[aa]=mass.calculate_mass(mass.Composition({}))
    return da_dict,db_dict,v_dict,wa_dict,wb_dict




def annotate_side_chain_losses(all_settings):
    data=all_settings["data"]
    proton_mass=all_settings["proton_mass"]
    frag_mass_tol=all_settings["mass_tol"]

    da_dict,db_dict,v_dict,wa_dict,wb_dict=create_side_chain_loss_mass_dict()
    pep_list=parser.parse(all_settings["pep_seq"],show_unmodified_termini=True,split=True)
    pep_charge=int(all_settings["pep_charge"])
    
    column="Side_chain_losses"
    data[column]=""

    v_ion=all_settings["v"]
    da_ion=all_settings["da"]
    db_ion=all_settings["db"]
    wa_ion=all_settings["wa"]
    wb_ion=all_settings["wb"]

    print all_settings["dvw_ions"]
    
    for i in xrange(1,len(pep_list)):
        abc_fragment=pep_list[:i]
        abc_ion_index=i
        xyz_fragment=pep_list[i:]
        xyz_ion_index=len(pep_list)-i        
        #print pep_list

        cterm_aa_at_break_point=abc_fragment[-1][-1]
        #print abc_fragment
        #print "-------------------------"
        #print xyz_fragment

        nterm_aa_at_break_point=xyz_fragment[0][-1]
        if "-" in nterm_aa_at_break_point:
            nterm_aa_at_break_point=xyz_fragment[0][-2]


        da_loss_mass=da_dict[cterm_aa_at_break_point]
        db_loss_mass=db_dict[cterm_aa_at_break_point]
        v_loss_mass=v_dict[nterm_aa_at_break_point]
        wa_loss_mass=wa_dict[nterm_aa_at_break_point]
        wb_loss_mass=wb_dict[nterm_aa_at_break_point]

        
        a_neutral_mass= mass.fast_mass2(abc_fragment,ion_type="a")
        y_neutral_mass= mass.fast_mass2(xyz_fragment,ion_type="y")
        z_neutral_mass= mass.fast_mass2(xyz_fragment,ion_type="z")

        if da_ion and da_loss_mass:
            base="da/"+str(abc_ion_index)
            da_neutral_mass=a_neutral_mass-da_loss_mass
            da_mono_mass=(da_neutral_mass+proton_mass)/1.0
            peak1=base+"/1+/0//0C13/"+str(da_mono_mass)
            data=check_and_apply_mass(data,da_mono_mass,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,da_neutral_mass,column,check_and_apply_mass,frag_mass_tol)

        if db_ion and db_loss_mass:
            base="db/"+str(abc_ion_index)
            db_neutral_mass=a_neutral_mass-db_loss_mass
            db_mono_mass=(db_neutral_mass+proton_mass)/1.0
            peak1=base+"/1+/0//0C13/"+str(db_mono_mass)
            data=check_and_apply_mass(data,db_mono_mass,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,db_neutral_mass,column,check_and_apply_mass,frag_mass_tol)

        if v_ion and v_loss_mass:

            base="v/"+str(xyz_ion_index)
            v_neutral_mass=y_neutral_mass-v_loss_mass
            v_mono_mass=(v_neutral_mass+proton_mass)/1.0
            peak1=base+"/1+/0//0C13/"+str(v_mono_mass)
            data=check_and_apply_mass(data,v_mono_mass,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,v_neutral_mass,column,check_and_apply_mass,frag_mass_tol)

        if wa_ion and wa_loss_mass:

            base="wa/"+str(xyz_ion_index)
            wa_neutral_mass=z_neutral_mass-wa_loss_mass
            wa_mono_mass=(wa_neutral_mass+proton_mass)/1.0
            peak1=base+"/1+/0//0C13/"+str(wa_mono_mass)
            data=check_and_apply_mass(data,wa_mono_mass,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,wa_neutral_mass,column,check_and_apply_mass,frag_mass_tol)


        if wb_ion and wb_loss_mass:

            base="wb/"+str(xyz_ion_index)
            wb_neutral_mass=z_neutral_mass-wb_loss_mass
            wb_mono_mass=(wb_neutral_mass+proton_mass)/1.0
            peak1=base+"/1+/0//0C13/"+str(wb_mono_mass)
            data=check_and_apply_mass(data,wb_mono_mass,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,wb_neutral_mass,column,check_and_apply_mass,frag_mass_tol)
  
    #empty_inde        
    data["Side_chain_losses"][data["Side_chain_losses"][data["Side_chain_losses"]==""].index]=np.nan
    
    all_settings["data"]=data
    print data.Side_chain_losses


    return all_settings




def get_comp_from_unimod():
    #db=mass.Unimod()
    aa_comp=dict(mass.std_aa_comp)
    #mass_h=mass.std_aa_mass["H-"]
    #aa_comp["p"]=db.by_title("Phospho")["composition"]
    aa_comp["p"]=mass.Composition({"H":1,"O":3,"P":1})
    
    #aa_comp["me"]=db.by_title("Methyl")["composition"]
    aa_comp["me"]=mass.Composition({"H":2,"C":1})
    
    
    #aa_comp["Me-"]=db.by_title("Methyl")["composition"]+mass.std_aa_comp["H-"]
    aa_comp["Me-"]=mass.Composition({"H":3,"C":1})
    
    #aa_comp["pyq"]=db.by_title("Gln->pyro-Glu")["composition"]
    aa_comp["pyq"]=mass.Composition({"H":-3,"N":-1})
    
    
    #aa_comp["Pyq-"]=db.by_title("Gln->pyro-Glu")["composition"]+mass.std_aa_comp["H-"]
    aa_comp["Pyq-"]=mass.Composition({"H":-2,"N":-1})

    #aa_comp["pye"]=db.by_title("Glu->pyro-Glu")["composition"]    
    aa_comp["pye"]=mass.Composition({"H":-2,"O":-1})
    
    #aa_comp["Pye-"]=db.by_title("Glu->pyro-Glu")["composition"]+mass.std_aa_comp["H-"]    
    aa_comp["Pye-"]=mass.Composition({"H":-1,"O":-1})
    
    #aa_comp["ox"]=db.by_title("Oxidation")["composition"]
    aa_comp["ox"]=mass.Composition({"O":1})
       
    
    #aa_comp["de"]=db.by_title("Deamidated")["composition"]
    aa_comp["de"]=mass.Composition({"H":-1,"N":-1,"O":1})
    
    
    #[aa_comp["am"]=db.by_title("Ammonia-loss")["composition"]
    aa_comp["am"]=mass.Composition({"N":-1,"H":-3})
    
    #aa_comp["Am-"]=db.by_title("Ammonia-loss")["composition"]+mass.std_aa_comp["H-"]
    aa_comp["Am-"]=mass.Composition({"N":-1,"H":-2})
    
    
    #aa_comp["ca"]=db.by_title("Carbamyl")["composition"]
    aa_comp["ca"]=mass.Composition({"H":1,"C":1,"N":1,"O":1})
    
    #aa_comp["Ca-"]=db.by_title("Carbamyl")["composition"]+mass.std_aa_comp["H-"]
    aa_comp["Ca-"]=mass.Composition({"H":2,"C":1,"N":1,"O":1})
    

    #aa_comp["ac"]=db.by_title("Acetyl")["composition"]
    aa_comp["ac"]=mass.Composition({"H":2,"C":2,"O":1})
     
    #aa_comp["Ac"]=db.by_title("Acetyl")["composition"]
    #aa_comp["Ac-"]=db.by_title("Acetyl")["composition"]+mass.std_aa_comp["H-"]
    aa_comp["Ac-"]=mass.Composition({"H":3,"C":2,"O":1})
    
    #aa_comp["cam"]=db.by_title("Carbamidomethyl")["composition"]
    aa_comp["cam"]=mass.Composition({"H":3,"C":2,"N":1,"O":1})
    
    #aa_comp["cam"]=db.by_title("Carbamidomethyl")["composition"]
    
    #aa_comp["oxn"]=db.by_title("Dethiomethyl")["composition"]
    aa_comp["oxn"]=mass.Composition({"H":-4,"C":-1,"S":-1})
    

    aa_comp["heavyk"]=mass.Composition({"H":8})

    aa_comp["heavyr"]=mass.Composition({"H":10})


    #print aa_comp
    for key in aa_comp:
        mass.std_aa_mass[key]=mass.calculate_mass(aa_comp[key])
    
    return aa_comp



def generate_ion_mass(ion_type,ion_index_str,frag_aa_comp,pep_fragment,column,all_settings):
    max_global_water_loss=all_settings["max_global_ammonia_loss"]
    max_global_ammonia_loss=all_settings["max_global_ammonia_loss"]
    proton_mass=all_settings["proton_mass"]
    data=all_settings["data"]
    min_mz=all_settings["min_mz"]
    max_mz=all_settings["max_mz"]
    pep_charge=all_settings["pep_charge"]
    frag_mass_tol=all_settings["mass_tol"]
    max_frag_charge=all_settings["max_frag_charge"]
    count_basic_residues=all_settings["count_basic_residues"]




    all_losses=get_multiple_neutral_losses(frag_aa_comp,max_global_water_loss,max_global_ammonia_loss)
    if all_losses:
        loss_charges=[1,2]
    else:
        loss_charges=None

    neutral_mass=mass.fast_mass2(pep_fragment,ion_type=ion_type)
    base=ion_type+"/"+ion_index_str

    if max_frag_charge>=1:

        ion_mass_charge1=(neutral_mass+proton_mass)/1.0
        if (ion_mass_charge1 > min_mz) and (ion_mass_charge1 < max_mz):
            
            peak1=base+"/1+/0//0C13/"+str(ion_mass_charge1)
            data=check_and_apply_mass(data,ion_mass_charge1,peak1,column,frag_mass_tol)
            peak1_iso=base+"/1+/0//"
            data=annotate_fragment_isotopes(data,1.0,peak1_iso,neutral_mass,column,check_and_apply_mass,frag_mass_tol)
        
        
    if max_frag_charge>=2:

        ion_mass_charge2=(neutral_mass+2*proton_mass)/2.0
        if (ion_mass_charge2 > min_mz) and (ion_mass_charge2 < max_mz):
            peak2=base+"/2+/0//0C13/"+str(ion_mass_charge2)
            data=check_and_apply_mass(data,ion_mass_charge2,peak2,column,frag_mass_tol)
            peak2_iso=base+"/2+/0//"
            data=annotate_fragment_isotopes(data,2.0,peak2_iso,neutral_mass,column,check_and_apply_mass,frag_mass_tol)

    for charge in range(3,max_frag_charge+1):
        #higher_loss_charges=[]
        if count_basic_residues:
            basic_residues=frag_aa_comp["R"]+frag_aa_comp["K"]+frag_aa_comp["H"]
            if charge>basic_residues:
                break

        if loss_charges:
            loss_charges.append(charge)
        ion_mass=(neutral_mass+charge*proton_mass)/charge
        if (ion_mass > min_mz) and (ion_mass < max_mz):
            peak=base+"/"+str(charge)+"+/0//0C13/"+str(ion_mass)
            data=check_and_apply_mass(data,ion_mass,peak,column,frag_mass_tol)

        
    if loss_charges:
        for charge in loss_charges:
            if charge<=max_frag_charge:
                data=apply_all_losses(data,neutral_mass,all_losses,charge,base,column,min_mz,max_mz)
                
                 #data=annotate_fragment_isotopes(data,1.0,peak1,neutral_mass,column)
                
            

    return data
        
def get_multiple_neutral_losses(frag_aa_comp,max_global_water_loss,max_global_ammonia_loss):
    #H-STSFQGGLGSRSoxMAAGoxMAGGLAGoxMGGIQNEKETMQSLNDR-OH
    #aa_comp=get_comp_from_unimod()
    #frag_aa_comp=parser.amino_acid_composition("oxMAAGoxMASTERDGLAGoxM",labels=aa_comp.keys())
    #print frag_aa_comp
    #max_neutral_loss=3
    #print frag_aa_comp
    aa_list=["R","K","N","Q","S","T","E","D","oxM"]
    for aa in aa_list:
        try:
            frag_aa_comp[aa]
        except:
            frag_aa_comp[aa]=0


    
    total_ammonia_loss=frag_aa_comp["R"]+frag_aa_comp["K"]+frag_aa_comp["N"]+frag_aa_comp["Q"]
    total_water_loss=frag_aa_comp["S"]+frag_aa_comp["T"]+frag_aa_comp["E"]+frag_aa_comp["D"]
    total_met_loss=frag_aa_comp["oxM"]
    
    total_neutral_loss=total_ammonia_loss+total_water_loss+total_met_loss
    #print total_ammonia_loss,total_water_loss  
    
    if total_neutral_loss==0:
        return None
        
    if max_global_water_loss=="All":
        max_water_loss=total_water_loss
        
    elif total_water_loss>max_global_water_loss:
        total_water_loss=max_global_water_loss
    else:
        pass
        
    if max_global_ammonia_loss=="All":
        max_ammonia_loss=total_ammonia_loss
    
    elif total_ammonia_loss>max_global_ammonia_loss:
        total_ammonia_loss=max_global_ammonia_loss
    else:
        pass
        
    #print total_ammonia_loss,total_water_loss
               
    ammonia_list,water_list,met_list=[],[],[]
    for k in range(1,total_ammonia_loss+1):
        ammonia_list.append((k,"NH3"))
       # ammonia_list.append(str(k)+":NH3")
    for l in range(1,total_water_loss+1):
        water_list.append((l,"H2O"))
        
    for m in range(1,total_met_loss+1):
        met_list.append((m,"CH4S"))
        
    
    if len(ammonia_list)==0 and len(water_list)==0 and len(met_list)==0:
        return None
        
        print "no neutral_loss"
       
        
        #return data,calculated_set
        
    elif len(ammonia_list)==0 or len(water_list)==0 or len(met_list)==0:
        if len(ammonia_list)==0 and len(met_list)==0:
            all_losses=tuple(itertools.product(water_list))
            
        elif (len(water_list)==0) and (len(met_list)==0):
            all_losses=tuple(itertools.product(ammonia_list))
            
        elif (len(water_list)==0) and (len(ammonia_list)==0):
            all_losses=tuple(itertools.product(met_list))
            
        elif len(met_list)==0:
            all_losses=list(itertools.product(ammonia_list))
            all_losses=all_losses+list(itertools.product(water_list))
            all_losses=tuple(all_losses+list(itertools.product(ammonia_list,water_list)))            
        
        elif len(ammonia_list)==0:
            all_losses=list(itertools.product(met_list))
            all_losses=all_losses+list(itertools.product(water_list))
            all_losses=tuple(all_losses+list(itertools.product(met_list,water_list)))            

        elif len(water_list)==0:
            all_losses=list(itertools.product(met_list))
            all_losses=all_losses+list(itertools.product(ammonia_list))
            all_losses=tuple(all_losses+list(itertools.product(met_list,ammonia_list)))
                        
                        
                
            #print "ammonia",all_losses
    elif len(ammonia_list)>=1 and len(water_list)>=1 and len(met_list)>=1:
        all_losses=list(itertools.product(ammonia_list))
        all_losses=all_losses+list(itertools.product(water_list))
        all_losses=all_losses+list(itertools.product(met_list))
        
        all_losses=all_losses+list(itertools.product(ammonia_list,water_list))
        all_losses=all_losses+list(itertools.product(met_list,water_list))
        all_losses=all_losses+list(itertools.product(met_list,ammonia_list))

                        
        all_losses=tuple(all_losses+list(itertools.product(ammonia_list,water_list,met_list)))
    else:
        print "unexpected!"
        #return data,calculated_set
        sys.exit()
        
    #print all_losses
    all_losses_short=[]
    for loss_comb in all_losses:
        loss_quantity=0
        for loss in loss_comb:
            loss_quantity+=loss[0]
            
        #if loss_quantity>max_neutral_loss:
         #   continue
        all_losses_short.append(loss_comb)
    
    return tuple(all_losses_short)
        
def apply_all_losses(data,neutral_mass,all_losses,charge,base,column,min_mz,max_mz):
    for each_set in all_losses:
        ion_mass=neutral_mass
        total_loss,total_loss_string=0,""
        loss_temp_list=[]
        for each_loss in each_set:
            loss_quantity,loss_type=each_loss[0],each_loss[1]
            total_loss+=loss_quantity
            each_loss_string=str(loss_quantity)+"-"+loss_type
            loss_temp_list.append(each_loss_string)
            ion_mass-=loss_quantity*loss_mass_dict[loss_type]
        #loss_temp_list.insert(0,str(total_loss))
                    
        ion_mz=(ion_mass+proton_mass*charge)/charge
        
        if (ion_mz > min_mz) and (ion_mz < max_mz):
            loss_string=str(total_loss)+"/"+(":".join(loss_temp_list))
            
            loss_peak=base+"/"+str(charge)+"+/"+loss_string+"/0C13/"+str(ion_mz)
            #print loss_peak
            #loss_peak=base+"/"+str(charge)+"+/"+("".join(str(each_set)))+"/"+str(ion_mass)
            data=check_and_apply_mass(data,ion_mz,loss_peak,column,frag_mass_tol)
            if charge<3:
                temp_base=base+"/"+str(charge)+"+/"+loss_string+"/"
                data=annotate_fragment_isotopes(data,charge,temp_base,ion_mass,column,check_and_apply_mass,frag_mass_tol)
            
            
        
    return data


def write_to_excel(all_settings):
    cwd=os.getcwd()
    dir_path=cwd+os.path.sep+"tmp"+os.path.sep
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    match_file_name="all_matches.xlsx"
    match_file_path=dir_path+match_file_name
    print match_file_path
    print os.path.exists(match_file_path)
    if os.path.exists(match_file_path):
        for i in range(1,1000):
            print i
            new_match_file="all_matches_"+str(i)+".xlsx"
            new_match_file_path=dir_path+new_match_file
            if os.path.exists(new_match_file_path):
                continue
            else:
                match_file_path=new_match_file_path
                break
        print match_file_path
    else:
        match_file_path=match_file_path

    print match_file_path
    print "--------"
    out_excel=pd.ExcelWriter(match_file_path)
    #print all_settings
    
    all_settings["data"].to_excel(out_excel,sheet_name="all_matches")
    out_excel.save()


