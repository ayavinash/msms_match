import os,sys
from Tkinter import Tk,Label,Button,Checkbutton,Entry,StringVar,DISABLED,NORMAL,END,W,E,Frame,BooleanVar,IntVar,Radiobutton
import ScrolledText
from ttk import Combobox
from tkFileDialog import askopenfilename
#import spectrum_plotter_oop_functions as library
from tkMessageBox import showerror,showinfo
import pandas as pd
from pyteomics import mass,parser
from StringIO import StringIO
import numpy as np
import pandas as pd
from spectrum_plotter_oop_functions import *
from writer import write_match_to_file
from tkColorChooser import askcolor


class Spectrum_plotter_GUI:
	def __init__(self,master):
		self.master=master
		master.title("/////////////////Match Generator//////////////////////",)
		self.peak_list_frame=Frame(master,width=100,height=50,background="blue")
		self.peak_list_frame.peak_list=ScrolledText.ScrolledText(self.peak_list_frame,width=40,height=50)
		self.peak_list_frame.label=Label(self.peak_list_frame,text="Enter Peak List Below")
		self.peak_list_frame.clear_button=Button(self.peak_list_frame,command=self.clear_peak_list,width=20,text="Clear")
		self.peak_list_frame.label.grid(row=0)
		self.peak_list_frame.peak_list.grid(row=1)
		self.peak_list_frame.clear_button.grid(row=2,padx=10,pady=10)

		#self.peak_list_frame.grid(row=0)


		self.parameter_frame=Frame(master)
		self.parameter_frame.pep_seq=StringVar()
		self.parameter_frame.pep_charge=StringVar()
		self.parameter_frame.pep_mz=StringVar()
		self.parameter_frame.mass_tol=StringVar()
		self.parameter_frame.mass_tol_unit=StringVar()
		self.parameter_frame.max_water_loss=StringVar()
		self.parameter_frame.max_ammonia_loss=StringVar()

		self.parameter_frame.count_basic_residues=BooleanVar()
		self.parameter_frame.annotate_c13_peaks=BooleanVar()
		self.parameter_frame.match_file_name=StringVar()

		
		self.parameter_frame.peptide_label=Label(self.parameter_frame,text="Peptide Sequence >")
		self.parameter_frame.peptide_entry=Entry(self.parameter_frame,textvariable=self.parameter_frame.pep_seq,width=80)
		self.parameter_frame.peptide_label.grid(row=0,column=1)
		self.parameter_frame.peptide_entry.grid(row=0,column=2,sticky="W")

		self.parameter_frame.pep_charge_label=Label(self.parameter_frame,text="Peptide Charge: ")
		self.parameter_frame.pep_charge_box=Combobox(self.parameter_frame,textvariable=self.parameter_frame.pep_charge,values=range(1,101))
		self.parameter_frame.pep_charge_label.grid(row=1,column=1)
		self.parameter_frame.pep_charge_box.grid(row=1,column=2,sticky="W")

		self.parameter_frame.masstol_label=Label(self.parameter_frame,text="Mass Tolerance: ")
		self.parameter_frame.masstol_entry=Entry(self.parameter_frame,textvariable=self.parameter_frame.mass_tol,width=10)
		self.parameter_frame.masstol_unit_label=Label(self.parameter_frame,text="Mass Tolerance Unit: ")
		self.parameter_frame.mass_tol_unit_box=Combobox(self.parameter_frame,textvariable=self.parameter_frame.mass_tol_unit,values=["Dalton","ppm"])
		self.parameter_frame.max_water_loss_box=Combobox(self.parameter_frame,textvariable=self.parameter_frame.max_water_loss,values=["All"]+range(0,101))
		self.parameter_frame.max_ammonia_loss_box=Combobox(self.parameter_frame,textvariable=self.parameter_frame.max_ammonia_loss,values=["All"]+range(0,101))
		self.parameter_frame.max_water_loss_label=Label(self.parameter_frame,text="Max water loss")
		self.parameter_frame.max_ammonia_loss_label=Label(self.parameter_frame,text="Max ammonia loss")
		self.parameter_frame.count_basic_residues_button=Checkbutton(self.parameter_frame,text="Count Basic Residues (RKH)",variable=self.parameter_frame.count_basic_residues)
		self.parameter_frame.annotate_c13_peaks_button=Checkbutton(self.parameter_frame,text="Annotate C13 isotopes",variable=self.parameter_frame.annotate_c13_peaks)
		self.parameter_frame.reset_all_button=Button(self.parameter_frame,command=self.set_defaults,text="Reset All",background="grey")
		self.parameter_frame.match_file_label=Label(self.parameter_frame,text="Match file name")
		self.parameter_frame.match_file_name_entry=Entry(self.parameter_frame,textvariable=self.parameter_frame.match_file_name)


		
		self.parameter_frame.masstol_label.grid(row=2,column=1)
		self.parameter_frame.masstol_unit_label.grid(row=3,column=1)
		self.parameter_frame.masstol_entry.grid(row=2,column=2,sticky="W")
		self.parameter_frame.mass_tol_unit_box.grid(row=3,column=2,sticky="W")
		self.parameter_frame.max_water_loss_label.grid(row=4,column=1,sticky="w")
		self.parameter_frame.max_water_loss_box.grid(row=4,column=2,sticky="w")
		self.parameter_frame.max_ammonia_loss_label.grid(row=5,column=1,sticky="W")
		self.parameter_frame.max_ammonia_loss_box.grid(row=5,column=2,sticky="W")
		self.parameter_frame.count_basic_residues_button.grid(row=6,column=1,sticky="W")
		self.parameter_frame.annotate_c13_peaks_button.grid(row=7,column=1,sticky="W")
		self.parameter_frame.reset_all_button.grid(row=8,column=0,sticky="W")
		self.parameter_frame.match_file_label.grid(row=9,column=0,sticky="w")
		self.parameter_frame.match_file_name_entry.grid(row=9,column=1,sticky="w")





		self.theo_frame=Frame(master,borderwidth=10)
		self.theo_frame.theo=BooleanVar()
		self.theo_frame.a=BooleanVar()
		self.theo_frame.a_plus1=BooleanVar()
		self.theo_frame.b=BooleanVar()
		self.theo_frame.c=BooleanVar()
		self.theo_frame.c_plus1=BooleanVar()
		self.theo_frame.c_plus2=BooleanVar()
		self.theo_frame.c_minus1=BooleanVar()

		self.theo_frame.x=BooleanVar()
		self.theo_frame.x_plus1=BooleanVar()
		self.theo_frame.y=BooleanVar()
		self.theo_frame.z=BooleanVar()
		self.theo_frame.z=BooleanVar()
		self.theo_frame.z_plus1=BooleanVar()
		self.theo_frame.z_plus2=BooleanVar()
		self.theo_frame.z_plus3=BooleanVar()

		self.theo_frame.max_charge=StringVar()
		self.theo_frame.max_loss=StringVar()
		self.theo_frame.color=StringVar()





		self.theo_frame.theo_button=Checkbutton(self.theo_frame,text="Ion Series",variable=self.theo_frame.theo)
		self.theo_frame.a_button=Checkbutton(self.theo_frame,text="a",variable=self.theo_frame.a)
		self.theo_frame.a_plus1_button=Checkbutton(self.theo_frame,text="a+1",variable=self.theo_frame.a_plus1)
		self.theo_frame.b_button=Checkbutton(self.theo_frame,text="b",variable=self.theo_frame.b)
		self.theo_frame.c_button=Checkbutton(self.theo_frame,text="c",variable=self.theo_frame.c)
		self.theo_frame.c_plus1_button=Checkbutton(self.theo_frame,text="c+1",variable=self.theo_frame.c_plus1)
		self.theo_frame.c_plus2_button=Checkbutton(self.theo_frame,text="c+2",variable=self.theo_frame.c_plus2)
		self.theo_frame.c_minus1_button=Checkbutton(self.theo_frame,text="c-1",variable=self.theo_frame.c_minus1)

		self.theo_frame.x_button=Checkbutton(self.theo_frame,text="x",variable=self.theo_frame.x)
		self.theo_frame.x_plus1_button=Checkbutton(self.theo_frame,text="x+1",variable=self.theo_frame.x_plus1)
		self.theo_frame.y_button=Checkbutton(self.theo_frame,text="y",variable=self.theo_frame.y)
		self.theo_frame.z_button=Checkbutton(self.theo_frame,text="z",variable=self.theo_frame.z)
		self.theo_frame.z_plus1_button=Checkbutton(self.theo_frame,text="z+1",variable=self.theo_frame.z_plus1)
		self.theo_frame.z_plus2_button=Checkbutton(self.theo_frame,text="z+2",variable=self.theo_frame.z_plus2)
		self.theo_frame.z_plus3_button=Checkbutton(self.theo_frame,text="z+3",variable=self.theo_frame.z_plus3)
		

		self.theo_frame.max_charge_label=Label(self.theo_frame,text="Max Charge")
		self.theo_frame.max_loss_label=Label(self.theo_frame,text="Max Neutral Loss")
		#self.theo_frame=theo_frame



		self.theo_frame.max_charge_box=Combobox(self.theo_frame,textvariable=self.theo_frame.max_charge,values=["All"]+range(1,101))
		self.theo_frame.max_loss_box=Combobox(self.theo_frame,textvariable=self.theo_frame.max_loss,values=["All"]+range(0,101))
		self.theo_frame.color_button=Button(self.theo_frame,text="Select color",command=self.get_theo_color)


		self.theo_frame.theo_button.grid(row=0)
		self.theo_frame.a_button.grid(row=1,column=0)
		self.theo_frame.b_button.grid(row=2,column=0)
		self.theo_frame.c_button.grid(row=3,column=0)
		self.theo_frame.c_plus1_button.grid(row=4,column=0)
		self.theo_frame.c_plus2_button.grid(row=5,column=0)
		self.theo_frame.c_minus1_button.grid(row=6,column=0)
		self.theo_frame.a_plus1_button.grid(row=7,column=0)

		self.theo_frame.x_button.grid(row=1,column=1)
		self.theo_frame.y_button.grid(row=2,column=1)
		self.theo_frame.z_button.grid(row=3,column=1)
		self.theo_frame.z_plus1_button.grid(row=4,column=1)
		self.theo_frame.z_plus2_button.grid(row=5,column=1)
		self.theo_frame.z_plus3_button.grid(row=6,column=1)
		self.theo_frame.x_plus1_button.grid(row=7,column=1)

		self.theo_frame.max_charge_label.grid(row=8,column=0)
		self.theo_frame.max_charge_box.grid(row=8,column=1)
		self.theo_frame.max_loss_label.grid(row=9,column=0)
		self.theo_frame.max_loss_box.grid(row=9,column=1)
		self.theo_frame.color_button.grid(row=11,column=0)





###Internal fragments

		self.int_frame=Frame(master,borderwidth=10)
		self.int_frame.int=BooleanVar()
		self.int_frame.a=BooleanVar()
		self.int_frame.b=BooleanVar()
		self.int_frame.max_charge=StringVar()
		self.int_frame.max_loss=StringVar()
		self.int_frame.max_len=StringVar()
		self.int_frame.nterm_is_p=BooleanVar()
		self.int_frame.cterm_is_deqwh=BooleanVar()	
		self.int_frame.color=StringVar()


		self.int_frame.int_button=Checkbutton(self.int_frame,text="Internal Ions",variable=self.int_frame.int)
		self.int_frame.a_button=Checkbutton(self.int_frame,text="a type",variable=self.int_frame.a)
		self.int_frame.b_button=Checkbutton(self.int_frame,text="b type",variable=self.int_frame.b)

		self.int_frame.max_charge_label=Label(self.int_frame,text="Max Charge")
		self.int_frame.max_loss_label=Label(self.int_frame,text="Max Neutral Loss")
		self.int_frame.max_len_label=Label(self.int_frame,text="Max Length")
		
		self.int_frame.p_button=Checkbutton(self.int_frame,text="N-term in P",variable=self.int_frame.nterm_is_p)
		self.int_frame.deqwh_button=Checkbutton(self.int_frame,text="C-term in (DEQWH)",variable=self.int_frame.cterm_is_deqwh)


		self.int_frame.max_charge_box=Combobox(self.int_frame,textvariable=self.int_frame.max_charge,values=["All"]+range(1,101))
		self.int_frame.max_loss_box=Combobox(self.int_frame,textvariable=self.int_frame.max_loss,values=["All"]+range(0,101))
		self.int_frame.max_len_box=Combobox(self.int_frame,textvariable=self.int_frame.max_len,values=["All"]+range(0,101))
		self.int_frame.color_button=Button(self.int_frame,text="Select Color",command=self.get_int_color)



		self.int_frame.int_button.grid(row=0)
		self.int_frame.a_button.grid(row=1,column=0)
		self.int_frame.b_button.grid(row=1,column=1)

		self.int_frame.max_charge_label.grid(row=2,column=0)
		self.int_frame.max_charge_box.grid(row=2,column=1)
		self.int_frame.max_loss_label.grid(row=3,column=0)
		self.int_frame.max_loss_box.grid(row=3,column=1)
		self.int_frame.max_len_label.grid(row=4,column=0)
		self.int_frame.max_len_box.grid(row=4,column=1)
		self.int_frame.p_button.grid(row=5,column=0)
		self.int_frame.deqwh_button.grid(row=5,column=1)
		self.int_frame.color_button.grid(row=6,column=0)



###Precursors

		self.pre_frame=Frame(master,borderwidth=10)
		self.pre_frame.pre=BooleanVar()
		self.pre_frame.remove_peaks=BooleanVar()
		self.pre_frame.min_charge=StringVar()
		self.pre_frame.max_loss=StringVar()
		self.pre_frame.color=StringVar()
		


		self.pre_frame.pre_button=Checkbutton(self.pre_frame,text="Precursors",variable=self.pre_frame.pre)
		self.pre_frame.remove_button=Checkbutton(self.pre_frame,text="Remove Precursors",variable=self.pre_frame.remove_peaks)

		self.pre_frame.min_charge_label=Label(self.pre_frame,text="Min Charge")
		self.pre_frame.max_loss_label=Label(self.pre_frame,text="Max Neutral Loss")
		
		self.pre_frame.min_charge_box=Combobox(self.pre_frame,textvariable=self.pre_frame.min_charge,values=["All"]+range(1,101))
		self.pre_frame.max_loss_box=Combobox(self.pre_frame,textvariable=self.pre_frame.max_loss,values=["All"]+range(0,101))


		self.pre_frame.pre_button.grid(row=0)
		self.pre_frame.color_button=Button(self.pre_frame,text="Select color",command=self.get_pre_color)

		self.pre_frame.min_charge_label.grid(row=1,column=0)
		self.pre_frame.min_charge_box.grid(row=1,column=1)
		self.pre_frame.max_loss_label.grid(row=2,column=0)
		self.pre_frame.max_loss_box.grid(row=2,column=1)
		self.pre_frame.remove_button.grid(row=3,column=0)
		self.pre_frame.color_button.grid(row=4,column=0)

		

##Immonium ions

		self.immo_frame=Frame(master,borderwidth=10)
		self.immo_frame.immo=BooleanVar()
		self.immo_frame.immo_button=Checkbutton(self.immo_frame,text="Immonium Ions",variable=self.immo_frame.immo)
		self.immo_frame.immo_button.pack()
		self.immo_frame.color=StringVar()
		self.immo_frame.color_button=Button(self.immo_frame,text="Select color",command=self.get_immo_color)
		self.immo_frame.color_button.pack()



## MAscot
		self.mascot_frame=Frame(master,borderwidth=10)
		self.mascot_frame.mascot=BooleanVar()
		self.mascot_frame.mascot_button=Checkbutton(self.mascot_frame,text="Mascot",variable=self.mascot_frame.mascot)
		self.mascot_frame.mascot_button.pack()
		self.mascot_frame.color=StringVar()
		self.mascot_frame.color_button=Button(self.mascot_frame,text="Select color",command=self.get_mascot_color)
		self.mascot_frame.color_button.pack()

		self.isotope_frame=Frame(master,borderwidth=10)
		self.isotope_frame.minimum_isotopes_label=Label(self.isotope_frame,text="Minimum isotope peaks")
		self.isotope_frame.minimum_isotopes=StringVar()
		self.isotope_frame.minimum_isotopes_box=Combobox(self.isotope_frame,textvariable=self.isotope_frame.minimum_isotopes,values=range(1,101))
		self.isotope_frame.strictly_monoisotopic=BooleanVar()
		self.isotope_frame.strictly_monoisotopic_button=Checkbutton(self.isotope_frame,text="Monoistopic peak should be present",variable=self.isotope_frame.strictly_monoisotopic)
		self.isotope_frame.minimum_peaks_to_consider_no_0C13_label=Label(self.isotope_frame,text="Minimum isotopes to allow missed monoisotopic")
		self.isotope_frame.minimum_peaks_to_consider_no_0C13=StringVar()
		self.isotope_frame.minimum_peaks_to_consider_no_0C13_box=Combobox(self.isotope_frame,textvariable=self.isotope_frame.minimum_peaks_to_consider_no_0C13,values=range(1,101))

		self.isotope_frame.minimum_isotopes_label.grid(row=0,column=0)
		self.isotope_frame.minimum_isotopes_box.grid(row=0,column=1)
		self.isotope_frame.strictly_monoisotopic_button.grid(row=1,column=0)
		self.isotope_frame.minimum_peaks_to_consider_no_0C13_label.grid(row=2,column=0)
		self.isotope_frame.minimum_peaks_to_consider_no_0C13_box.grid(row=2,column=1)




		self.dvw_frame=Frame(master,borderwidth=10)
		self.dvw_frame.dvw=BooleanVar()
		self.dvw_frame.da=BooleanVar()
		self.dvw_frame.db=BooleanVar()
		self.dvw_frame.v=BooleanVar()
		self.dvw_frame.wa=BooleanVar()
		self.dvw_frame.wb=BooleanVar()

		self.dvw_frame.dvw_button=Checkbutton(self.dvw_frame,text="Side chain fragments",variable=self.dvw_frame.dvw)
		self.dvw_frame.da_button=Checkbutton(self.dvw_frame,text="da",variable=self.dvw_frame.da)
		self.dvw_frame.db_button=Checkbutton(self.dvw_frame,text="db",variable=self.dvw_frame.db)
		self.dvw_frame.v_button=Checkbutton(self.dvw_frame,text="v",variable=self.dvw_frame.v)
		self.dvw_frame.wa_button=Checkbutton(self.dvw_frame,text="wa",variable=self.dvw_frame.wa)
		self.dvw_frame.wb_button=Checkbutton(self.dvw_frame,text="wb",variable=self.dvw_frame.wb)
		
		self.dvw_frame.dvw_button.grid(row=0,column=0)
		self.dvw_frame.da_button.grid(row=1,column=0)
		self.dvw_frame.db_button.grid(row=1,column=1)
		self.dvw_frame.v_button.grid(row=2,column=0)
		self.dvw_frame.wa_button.grid(row=3,column=0)
		self.dvw_frame.wb_button.grid(row=3,column=1)







		self.view_frame=Frame(master,borderwidth=10)
		self.view_frame.within_label=Label(self.view_frame,text="Do not Label Peaks within")
		self.view_frame.within=StringVar()
		self.view_frame.within_entry=Entry(self.view_frame,textvariable=self.view_frame.within)

		self.view_frame.noise=StringVar()
		self.view_frame.noise_label=Label(self.view_frame,text="Noise level")
		self.view_frame.noise_entry=Entry(self.view_frame,textvariable=self.view_frame.noise)

		self.view_frame.match_type=IntVar()
		self.view_frame.match_type_loss_button=Radiobutton(self.view_frame,text="Select peak with least neutral loss",variable=self.view_frame.match_type,value=1)
		self.view_frame.match_type_mass_button=Radiobutton(self.view_frame,text="Select peak with least mass deviation",variable=self.view_frame.match_type,value=2)

		self.view_frame.apply_to_all=BooleanVar()
		self.view_frame.apply_to_all_button=Checkbutton(self.view_frame,text="Apply above rule to all groups",variable=self.view_frame.apply_to_all)



		self.view_frame.show_sequence=BooleanVar()
		self.view_frame.show_mass=BooleanVar()
		self.view_frame.show_mass_dev=BooleanVar()
		self.view_frame.show_mass_dev_ppm=BooleanVar()
		
		self.view_frame.show_sequence_button=Checkbutton(self.view_frame,text="Show sequence",variable=self.view_frame.show_sequence)
		self.view_frame.show_mass_button=Checkbutton(self.view_frame,text="Show m/z",variable=self.view_frame.show_mass)
		self.view_frame.show_mass_dev_button=Checkbutton(self.view_frame,text="Show mass dev",variable=self.view_frame.show_mass_dev)
		self.view_frame.show_mass_dev_ppm_button=Checkbutton(self.view_frame,text="Show mass dev ppm",variable=self.view_frame.show_mass_dev_ppm)





		self.view_frame.noise_label.grid(row=0,column=0)
		self.view_frame.noise_entry.grid(row=0,column=1)
		self.view_frame.within_label.grid(row=1,column=0)
		self.view_frame.within_entry.grid(row=1,column=1)
		self.view_frame.match_type_loss_button.grid(row=2,column=0)
		self.view_frame.match_type_mass_button.grid(row=3,column=0)
		self.view_frame.apply_to_all_button.grid(row=4,column=0)
		self.view_frame.show_sequence_button.grid(row=5,column=0)
		self.view_frame.show_mass_button.grid(row=5,column=1)
		self.view_frame.show_mass_dev_button.grid(row=6,column=0)
		self.view_frame.show_mass_dev_ppm_button.grid(row=6,column=1)


		BOLD=("Courier","24","bold")
		self.update_button=Button(master,command=self.run,text="Plot",background="red",font=BOLD)

		


####Load match from file


		

		self.file_frame=Frame(master)
		self.file_frame.match_file_path=StringVar()
		select_file_label=Label(self.file_frame,text="Select Match file: ")
		select_file_entry=Entry(self.file_frame,textvariable=self.file_frame.match_file_path,width=50)
		select_file_button=Button(self.file_frame,text="Browse",command=self.get_match_file)
		select_file_label.grid(row=0,column=0)
		select_file_entry.grid(row=0,column=2)
		select_file_button.grid(row=0,column=4)





		self.peak_list_frame.pack(side="left")
		self.parameter_frame.pack(side="top")
		self.theo_frame.pack(side="left")
		self.int_frame.pack(side="left")
		self.immo_frame.pack(side="left")
		self.mascot_frame.pack(side="left")
		self.dvw_frame.pack(side="left")
		self.isotope_frame.pack(side="bottom")
		self.pre_frame.pack(side="bottom")
		self.view_frame.pack(side="bottom")
		self.update_button.pack(side="bottom")
		self.file_frame.pack()

		self.set_defaults()



	def clear_peak_list(self):
		print "Hello"
		self.peak_list_frame.peak_list.delete(1.0,END)



	def update_plot(self):
		print "plotting"



	def get_match_file(self):
	    
	    match_file_path_input=askopenfilename(title="Select file",initialdir="../tmp/")
	    print match_file_path_input
	    self.file_frame.match_file_path.set(match_file_path_input)
	    match_file_path_input=self.file_frame.match_file_path.get()
	    if not match_file_path_input:
	        fileMessage=showerror("Input File Error","Select a match file to plot !")
	    

        

	def set_defaults(self):
		self.parameter_frame.pep_seq.set("H-GTAAAAAAAAAAAAAKVPAK-OH")
		self.parameter_frame.pep_charge.set("3")
		self.parameter_frame.mass_tol.set("0.6")
		self.parameter_frame.mass_tol_unit.set("da")
		self.parameter_frame.max_water_loss.set("1")
		self.parameter_frame.max_ammonia_loss.set("1")
		self.parameter_frame.count_basic_residues.set(1)
		self.parameter_frame.annotate_c13_peaks.set(0)
		self.parameter_frame.match_file_name.set("input_spectra.txt")


		self.theo_frame.theo.set(1)
		self.theo_frame.a.set(1)
		self.theo_frame.a_plus1.set(0)

		self.theo_frame.b.set(1)
		self.theo_frame.c.set(0)
		self.theo_frame.c_plus1.set(0)
		self.theo_frame.c_plus2.set(0)
		self.theo_frame.c_minus1.set(0)

		self.theo_frame.x.set(0)
		self.theo_frame.x_plus1.set(0)

		self.theo_frame.y.set(1)
		self.theo_frame.z.set(0)
		self.theo_frame.z_plus1.set(0)
		self.theo_frame.z_plus2.set(0)
		self.theo_frame.z_plus3.set(0)

		self.theo_frame.max_charge.set("All")
		self.theo_frame.max_loss.set("1")
		#self.theo_frame.theo_color.set((0,0,0))


		self.int_frame.int.set(1)
		self.int_frame.a.set(1)
		self.int_frame.b.set(1)
		self.int_frame.max_charge.set("All")
		self.int_frame.max_loss.set("1")
		self.int_frame.max_len.set("All")
		self.int_frame.nterm_is_p.set(0)
		self.int_frame.cterm_is_deqwh.set(0)

		self.mascot_frame.mascot.set(0)
		self.immo_frame.immo.set(0)

		self.pre_frame.pre.set(1)
		self.pre_frame.min_charge.set("All")
		self.pre_frame.max_loss.set(1)
		self.pre_frame.remove_peaks.set(0)

		self.view_frame.noise.set("DNL")
		self.view_frame.within.set("0.5")
		self.view_frame.match_type.set(1)
		self.view_frame.apply_to_all.set(0)

		self.view_frame.show_sequence.set(0)
		self.view_frame.show_mass.set(0)
		self.view_frame.show_mass_dev.set(0)
		self.view_frame.show_mass_dev_ppm.set(0)

		self.theo_frame.color.set("#000000")
		self.int_frame.color.set("#006400")
		self.pre_frame.color.set("#f4a460")
		self.mascot_frame.color.set("#ff0000")
		self.immo_frame.color.set("#ff00ff")


		self.isotope_frame.minimum_isotopes.set("3")
		self.isotope_frame.strictly_monoisotopic.set(1)
		self.isotope_frame.minimum_peaks_to_consider_no_0C13.set("5")

		self.dvw_frame.dvw.set(0)
		self.dvw_frame.da.set(0)
		self.dvw_frame.db.set(0)
		self.dvw_frame.v.set(0)
		self.dvw_frame.wa.set(0)
		self.dvw_frame.wb.set(0)





	def get_header_from_file(self,match_file_path):
	    my_file=open(match_file_path,"rb")
	    title=my_file.next().rstrip("\r\n").lstrip("#").split("Title:")[1]
	    source=my_file.next().rstrip("\r\n").lstrip("#").split("Source:")[1]
	    peptide=my_file.next().rstrip("\r\n").lstrip("#").split("Peptide:")[1]
	    modx_seq=my_file.next().rstrip("\r\n").lstrip("#").split("Modx_seq:")[1]
	    charge=my_file.next().rstrip("\r\n").lstrip("#").split("Charge:")[1]
	    mz=my_file.next().rstrip("\r\n").lstrip("#").split("ExpMz:")[1]
	    score=my_file.next().rstrip("\r\n").lstrip("#").split("Score:")[1]
	    interference=my_file.next().rstrip("\r\n").lstrip("#").split("Interference:")[1]
	    noise_type=my_file.next().rstrip("\r\n").lstrip("#").split("Noise_type:")[1]
	    noise_level=my_file.next().rstrip("\r\n").lstrip("#").split("Noise_Level:")[1]
	    my_file.close()
	    return (title,mz,peptide,charge,score,modx_seq,interference,noise_level,noise_type)


	def get_header_from_program(self):
		if not self.parameter_frame.pep_seq.get():
			showerror("Input Error","Provide a peptide sequence!")
		if not self.parameter_frame.pep_charge.get():
			showerror("Input Error","Provide charge of precursor!")

		charge=self.parameter_frame.pep_charge.get()
		modx_seq=self.parameter_frame.pep_seq.get()
		title="Peak list"
		mz="some mz"
		score="Some Score"
		interference="Not known"

		noise_level= "Not known"
		peptide=modx_seq[2:-3]

		return title,mz,peptide,charge,score,modx_seq,interference,noise_level




	def create_data_frame_from_peak_list(self,all_settings):
		print all_settings["peak_list"]
		#sys.exit()

		data=StringIO(all_settings["peak_list"])
		#print data
		#if " " in data:
		#	data=pd.read_table(data,names=["mz","intensity"],sep=" ")
		#if "\t" in data:
		data=pd.read_table(data,names=["mz","intensity"],sep="\t")


		#print data
		#print data.intensity
		#sys.exit()


		#data["Noise_level"]=np.nan
		data["Mascot"]=""
		data["Theoretical"]=""
		data["Precursors"]=""
		data["Internal_fragments"]=""
		data["Immonium_ions"]=""
		data["Side_chain_losses"]=""

		if all_settings["noise_type"]=="DNL":
			data,Noise_level=annotate_signal_noise_paper(data)
			all_settings["noise_level"]=Noise_level
		
		original_data=data.copy()
		all_settings["original_data"]=original_data
		max_mz,min_mz=data.mz.iloc[-1],data.mz.iloc[0]
		all_settings["max_mz"]=max_mz
		all_settings["min_mz"]=min_mz
		all_settings["data"]=data
		
		return all_settings


	def get_current_settings(self):
		all_settings={}
		noise=self.view_frame.noise.get()
		
		if noise!="DNL":
			try:
				float(noise)
				all_settings["noise_level"]=float(noise)
				all_settings["noise_type"]="Absolute"

			except:
				showerror("ValueError","Provide a meaningful noise level!")
		else:
			all_settings["noise_type"]="DNL"



		within=self.view_frame.within.get()
		if within!="5.0":
			try:
				float(within)
			except:
				showerror("Value Error","Provide a value >=0 for label display")
		all_settings["within"]=within
		
		mass_tol=self.parameter_frame.mass_tol.get()

		try:
			float(mass_tol)
		except:
			showerror("Value Error","Provide a mass tolerance!")
		all_settings["mass_tol"]=float(mass_tol)
		all_settings["mass_tol_unit"]=self.parameter_frame.mass_tol_unit.get()



		match_file_path=self.file_frame.match_file_path.get()

		peak_list=self.peak_list_frame.peak_list.get(1.0,END)
		print "peak_list_type",type(peak_list)

		if not match_file_path and not peak_list:
		    showerror("Input Error","Provide peak list or Select a match file !")
		if match_file_path:
			all_settings["match_file_path"]=match_file_path
			all_settings["Mascot"]=self.mascot_frame.mascot.get()
		else:
			all_settings["match_file_path"]=0
			self.mascot_frame.mascot.set(0)
			all_settings["Mascot"]=0
		all_settings["peak_list"]=peak_list
		all_settings["pep_seq"]=self.parameter_frame.pep_seq.get()
		all_settings["peptide"]=get_peptide_from_pep_seq(all_settings["pep_seq"])
		all_settings["pep_charge"]=int(self.parameter_frame.pep_charge.get())
		
		if all_settings["pep_charge"]>2:
			max_frag_charge=all_settings["pep_charge"]-1
		else:
			max_frag_charge=all_settings["pep_charge"]
		all_settings["max_frag_charge"]=max_frag_charge


		all_settings["mass_tol"]=self.parameter_frame.mass_tol.get()
		all_settings["mass_tol_unit"]=self.parameter_frame.mass_tol_unit.get()
		all_settings["Theoretical"]=self.theo_frame.theo.get()
		all_settings["Internal_fragments"]=self.int_frame.int.get()
		all_settings["Precursors"]=self.pre_frame.pre.get()
		all_settings["Immonium_ions"]=self.immo_frame.immo.get()
		#all_settings["Noise_level"]=noise
		all_settings["Theo_max_loss"]=self.theo_frame.max_loss.get()
		all_settings["Theo_max_charge"]=self.theo_frame.max_charge.get()
		all_settings["Int_max_loss"]=self.int_frame.max_loss.get()
		all_settings["Int_max_charge"]=self.int_frame.max_charge.get()
		all_settings["Int_max_len"]=self.int_frame.max_len.get()
		all_settings["Int_atype"]=self.int_frame.a.get()
		all_settings["Int_btype"]=self.int_frame.b.get()
		all_settings["Int_nterm_P"]=self.int_frame.nterm_is_p.get()
		all_settings["Int_cterm_DEQWH"]=self.int_frame.cterm_is_deqwh.get()
		all_settings["Pre_max_loss"]=self.pre_frame.max_loss.get()
		all_settings["Pre_remove_peaks"]=self.pre_frame.remove_peaks.get()
		all_settings["Pre_min_charge"]=self.pre_frame.min_charge.get()
		all_settings["max_global_water_loss"]=self.parameter_frame.max_water_loss.get()
		all_settings["max_global_ammonia_loss"]=self.parameter_frame.max_ammonia_loss.get()
		all_settings["count_basic_residues"]=self.parameter_frame.count_basic_residues.get()
		all_settings["match_type"]=self.view_frame.match_type.get()
		all_settings["apply_to_all"]=self.view_frame.apply_to_all.get()
		all_settings["annotate_c13_peaks"]=self.parameter_frame.annotate_c13_peaks.get()
		all_settings["match_file_name"]=self.parameter_frame.match_file_name.get()


		all_settings["show_mass"]=self.view_frame.show_mass.get()
		all_settings["show_sequence"]=self.view_frame.show_sequence.get()
		all_settings["show_mass_dev"]=self.view_frame.show_mass_dev.get()
		all_settings["show_mass_dev_ppm"]=self.view_frame.show_mass_dev_ppm.get()
		all_settings["theo_color"]=self.theo_frame.color.get()
		all_settings["int_color"]=self.int_frame.color.get()
		all_settings["mascot_color"]=self.mascot_frame.color.get()
		all_settings["pre_color"]=self.pre_frame.color.get()
		all_settings["immo_color"]=self.immo_frame.color.get()

		all_settings["minimum_isotopes"]=self.isotope_frame.minimum_isotopes.get()
		all_settings["strictly_monoisotopic"]=self.isotope_frame.strictly_monoisotopic.get()
		all_settings["minimum_peaks_to_consider_no_0C13"]=self.isotope_frame.minimum_peaks_to_consider_no_0C13.get()

		all_settings["Side_chain_losses"]=self.dvw_frame.dvw.get()

		if not all_settings["Side_chain_losses"]:
			all_settings["da"]=self.dvw_frame.da.set(0)
			all_settings["db"]=self.dvw_frame.db.set(0)
			all_settings["v"]=self.dvw_frame.v.set(0)
			all_settings["wa"]=self.dvw_frame.wa.set(0)
			all_settings["wb"]=self.dvw_frame.wb.set(0)
		else:
			all_settings["da"]=self.dvw_frame.da.get()
			all_settings["db"]=self.dvw_frame.db.get()
			all_settings["v"]=self.dvw_frame.v.get()
			all_settings["wa"]=self.dvw_frame.wa.get()
			all_settings["wb"]=self.dvw_frame.wb.get()

		dvw_ions_list=["da","db","v","wa","wb"]
		dvw_ions=[]
		for dvw_ion in dvw_ions_list:
			if all_settings[dvw_ion]:
				dvw_ions.append(dvw_ion)
		all_settings["dvw_ions"]=dvw_ions				

		nterm_series,cterm_series,int_series=set(),set(),set()
		if self.theo_frame.a.get():
			nterm_series.add("a")
		if self.theo_frame.b.get():
			nterm_series.add("b")
		if self.theo_frame.c.get():
			nterm_series.add("c")
		if self.theo_frame.c_plus1.get():
			nterm_series.add("c+1")
		if self.theo_frame.c_plus2.get():
			nterm_series.add("c+2")
		if self.theo_frame.c_minus1.get():
			nterm_series.add("c-1")
		if self.theo_frame.a_plus1.get():
			nterm_series.add("a+1")


		if self.theo_frame.x.get():
			cterm_series.add("x")
		if self.theo_frame.y.get():
			cterm_series.add("y")
		if self.theo_frame.z.get():
			cterm_series.add("z")

		if self.theo_frame.z_plus1.get():
			cterm_series.add("z+1")
		if self.theo_frame.z_plus2.get():
			cterm_series.add("z+2")
		if self.theo_frame.z_plus3.get():
			cterm_series.add("z+3")
		if self.theo_frame.x_plus1.get():
			cterm_series.add("x+1")
			
		
		if self.int_frame.a.get():
			int_series.add("a")
		if self.int_frame.b.get():
			int_series.add("b")

		all_settings["abc"]=tuple(nterm_series)
		all_settings["xyz"]=tuple(cterm_series)
		all_settings["ab_int"]=tuple(int_series)

		return all_settings


	def run(self):

		all_settings=self.get_current_settings()
		print all_settings.keys()
		print all_settings["strictly_monoisotopic"]
		#print all_settings
		#print all_settings["data"]
		initiate_globals(all_settings)
		
		if all_settings["peak_list"]!="\n":
			all_settings=self.create_data_frame_from_peak_list(all_settings)
			
			if all_settings["Theoretical"]:
				all_settings=annotate_theoretical(all_settings)
			if all_settings["Immonium_ions"]:
				all_settings=annotate_immonium_ions(all_settings)
			if all_settings["Internal_fragments"]:
				all_settings=annotate_internal_fragments(all_settings)
			if all_settings["Precursors"]:
				all_settings=annotate_precorsor(all_settings)
			if all_settings["Side_chain_losses"]:
				all_sttings=annotate_side_chain_losses(all_settings)


			#print all_settings["data"]

			title,mz,peptide,charge,score,modx_seq,interference,noise_level=self.get_header_from_program()
			all_settings["title"]=title
			all_settings["mz"]=mz
			all_settings["score"]=score
			all_settings["interference"]=interference

			write_match_to_file(all_settings)

		elif all_settings["match_file_path"]:
				#try:

			title,mz,peptide,charge,score,modx_seq,interference,noise_level,noise_type=self.get_header_from_file(all_settings["match_file_path"])

			#data=pd.read_table(all_settings["match_file_path"],comment="#",names=["mz","intensity","Noise_level","Mascot","Theoretical","Precursors","Immonium_ions","Internal_fragments"],dtype={"mz":np.float64,"intensity":np.float64,
			#	"Noise_level":np.float64,"Mascot":np.str,"Theoretical":np.str,"Precursors":np.str,"Immonium_ions":np.str,"Internal_fragments":np.str})
			data=pd.read_table(all_settings["match_file_path"],comment="#")
			all_settings["data"]=data
			all_settings["original_data"]=data
			all_settings["pep_seq"]=modx_seq
			self.parameter_frame.pep_seq.set(modx_seq)

			all_settings["pep_charge"]=int(charge)
			self.parameter_frame.pep_charge.set(int(charge))
			
			all_settings["title"]=title
			all_settings["mz"]=mz
			all_settings["score"]=score
			all_settings["interference"]=interference
			all_settings["noise_type"]=noise_type
			if "noise_level" not in all_settings:
				all_settings["noise_level"]=noise_level

			all_settings["min_mz"]=min(data.mz)
			all_settings["max_mz"]=max(data.mz)
			print data.columns
			if "Theoretical" not in data.columns:
				if all_settings["Theoretical"]:
					data["Theoretical"]=""
					all_settings=annotate_theoretical(all_settings)
					
			if "Immonium_ions" not in data.columns:
				if all_settings["Immonium_ions"]:
					data["Immonium_ions"]=""
					all_settings=annotate_immonium_ions(all_settings)
			if "Internal_fragments" not in data.columns:
				if all_settings["Internal_fragments"]:
					data["Internal_fragments"]=""
					all_settings=annotate_internal_fragments(all_settings)
			if "Precursors" not in data.columns:
				if all_settings["Precursors"]:
					data["Precursors"]=""
					all_settings=annotate_precorsor(all_settings)
			if "Side_chain_losses" not in data.columns:
				if all_settings["Side_chain_losses"]:
					data["Side_chain_losses"]=""
					all_settings=annotate_side_chain_losses(all_settings)

					
		else:
			showerror("Input Error","Provide peak list or a match file to procede!")

		all_settings=filter_by_annotation(all_settings)



		if all_settings["apply_to_all"]:
			all_settings=assign_annotation_to_all(all_settings)
		else:
			if all_settings["match_type"]==1:
				all_settings=assign_annotation_by_loss(all_settings)
			else:
				all_settings=assign_annotation_by_deviation(all_settings)


		spectrum_plotter(all_settings)


		self.theo_frame.color.set("#000000")
		self.int_frame.color.set("#006400")
		self.pre_frame.color.set("#f4a460")
		self.mascot_frame.color.set("#ff0000")
		self.immo_frame.color.set("#ff00ff")




	def get_theo_color(self):
	    new_color=askcolor()
	    print new_color
	    if new_color[0]:
	    	self.theo_frame.color.set(new_color[1])
	    else:
    		self.theo_frame.color.set("#000000")

	def get_int_color(self):
	    new_color=askcolor()
	    print new_color
	    if new_color[0]:
	    	self.int_frame.color.set(new_color[1])
	    else:
    		self.int_frame.color.set("#006400")

	def get_pre_color(self):
	    new_color=askcolor()
	    print new_color
	    if new_color[0]:
	    	self.pre_frame.color.set(new_color[1])
	    else:
    		self.pre_frame.color.set("#f4a460")
		
	def get_mascot_color(self):
	    new_color=askcolor()
	    print new_color
	    if new_color[0]:
	    	self.mascot_frame.color.set(new_color[1])
	    else:
	    	self.mascot_frame.color.set("#ff0000")
	    	
	def get_immo_color(self):
	    new_color=askcolor()
	    print new_color
	    if new_color[0]:
	    	self.immo_frame.color.set(new_color[1])
	    else:
    		self.immo_frame.color.set("#ff00ff")


if __name__=="__main__":
	root=Tk()
	Spectrum_plotter_gui=Spectrum_plotter_GUI(root)
	root.mainloop()
	print "got here"


