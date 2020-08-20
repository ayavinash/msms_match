import pandas as pd
import numpy as np




def check_and_apply_mass_dalton(data,ion_mass,peak,column,frag_mass_tol):
    temp=data[abs(data.mz-ion_mass)<frag_mass_tol]
    data=check_and_apply_mass_data(temp,data,ion_mass,peak,column)
    return data
       
def check_and_apply_mass_ppm(data,ion_mass,peak,column,frag_mass_tol):
   
    temp=data[(abs(data.mz-ion_mass)/ion_mass)*1000000<frag_mass_tol]
    data=check_and_apply_mass_data(temp,data,ion_mass,peak,column)
    return data
 
def check_and_apply_mass_data(temp,data,ion_mass,peak,column):

    if temp.empty:
        return data
    elif len(temp)==1:
      
        j=temp.index[0]
        
    else:
        j=np.argmax(temp.intensity)
    exp_mz=data.get_value(j,"mz")
    ann=data.get_value(j,column)
    mass_dev_da=exp_mz-ion_mass
    mass_dev_ppm=(mass_dev_da/ion_mass)*1000000

    if ann:
        #print ann,len(ann)

        data.set_value(j,column,ann+"|"+peak+"/"+str(mass_dev_da)+"/"+str(mass_dev_ppm))
    else:
        data.set_value(j,column,peak+"/"+str(mass_dev_da)+"/"+str(mass_dev_ppm))
    
    return data


def check_and_apply_pre_mass_dalton(data,ion_mass,peak,frag_mass_tol):
    temp=data[abs(data.mz-ion_mass)<frag_mass_tol]
    data=check_and_apply_pre_data(temp,data,ion_mass,peak)
    return data


def check_and_apply_pre_mass_ppm(data,ion_mass,peak,frag_mass_tol):

    temp=data[(abs(data.mz-ion_mass)/ion_mass)*1000000<frag_mass_tol]
    data=check_and_apply_pre_data(temp,data,ion_mass,peak)
    return data

def check_and_apply_pre_data(temp,data,ion_mass,peak):
   
    if temp.empty:
        return data
    i=temp.index
    for j in i:
        exp_mz=data.get_value(j,"mz")
        mass_dev_dalton=exp_mz-ion_mass
        mass_dev_ppm=(mass_dev_dalton/ion_mass)*1000000
        ann=data.get_value(j,"Precursors")
        if not ann:
            data.set_value(j,"Precursors",peak+"/"+str(mass_dev_dalton)+"/"+str(mass_dev_ppm))
        else:
            data.set_value(j,"Precursors",ann+"|"+peak+"/"+str(mass_dev_dalton)+"/"+str(mass_dev_ppm))

    return data



