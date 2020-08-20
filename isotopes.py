global PROTON_MASS
PROTON_MASS=1.00727647



def get_isotope_peaks(pep_exp_mr):
    isotope_peaks=3
    pep_exp_mr=float(pep_exp_mr)
    if pep_exp_mr<1000.0:
        isotope_peaks=3
    elif 1000.0<=pep_exp_mr<2000.0:
        isotope_peaks=4
    elif 2000.0<=pep_exp_mr<3000.0:
        isotope_peaks=5
    elif 3000.0<=pep_exp_mr<4000.0:
        isotope_peaks=6
    elif 4000.0<=pep_exp_mr<5000.0:
        isotope_peaks=7
    elif 5000.0<=pep_exp_mr<6000.0:
        isotope_peaks=8
    elif 6000.0<=pep_exp_mr<7000.0:
        isotope_peaks=9
    elif pep_exp_mr>=7000.0:
        isotope_peaks=10
    else:
        isotope_peaks=10
        
    return isotope_peaks    
 



def annotate_fragment_isotopes_none(data,charge,base,neutral_mass,column,check_and_apply_mass,frag_mass_tol):
    return data

def annotate_fragment_isotopes_all(data,charge,base,neutral_mass,column,check_and_apply_mass,frag_mass_tol):

    isotopes=get_isotope_peaks(neutral_mass)
    for isotope in xrange(1, isotopes+1):
        c13_1mz=(neutral_mass+isotope*PROTON_MASS+PROTON_MASS*charge)/charge
        c13_1peak=base+str(isotope)+"C13/"+str(c13_1mz)
        data=check_and_apply_mass(data,c13_1mz,c13_1peak,column,frag_mass_tol)
    return data
