# FH_group_self_ion_imp_W_TGS_github_V1
Data and code archive for the Self-ion implanted tungsten TGS paper - Thermal diffusivity degradation and point defect density in self-ion implanted tungsten 
Hi there !

This is a quick readme to the complete data set and processing scripts for all the data used in the study - Thermal diffusivity degradation and point defect density 
in self-ion implanted tungsten"

Please feel free to contact the authors regarding any issues encountered -
Felix Hofmann (feliox.hofmann@eng.ox.ac.uk)
Abdallah Reza (mohamed.reza@eng.ox.ac.uk)

----------------------------------------------------------------------------------------------------------------------------------------------------------

Figure 2b

figure2b_code.m uses the SRIM_20MEV_data_new.mat file, which is the SRIM output data, for simulations detailed in section 2.1. It then plots the damage and
implanted ion density profiles for the 1 dpa sample. Superimposed is the TGS probing depth.



Figure 3 

The Raw TGS data is given. 

map_fitting_2d_archive code will do the fitting of the data and produce a set of plots and a processed data file, simillar to the ones in 
the processed data subfolder. This is for each sample. 


summary_map_helsinki.m code, will then use these processed data files, i.e. the thermal diffusivity maps, and give the thermal diffusivity vs. dose plot 
and data set - helsinki_summary_data_4_4.mat

summary_plots_with_refs_archive.m then uses this data set, with data from literature, to plot the final FIGURE 3. 




Figure 4

analysis_TC_vs_dose_archive_fig_4.m uses the summary data (thermal diffusivity vs dose data) from figure 3, and calculates the defect densities 
using the KT model detailed in the paper section 3.2 and then gives FIGURE 4



Figure 5

analysis_TC_vs_dose_archive_figure_5.m uses the thermal diffusivity vs dose data set, simillar to in Figure 4, and calculates the TGS + KT defect density prediction.

It then includes the MD +TEM defect densities, from the data file felix_calc_MDTEM.mat, to give the comparison plot, that is FIGURE 5. The felix_calc_MDTEM.mat file
is what was obtained from the calculation method detailed in the paper section 3.3 and implemented in the excel sheet given (MD_TEM_comparison_calculation_Felix).


Figure 6

analysis_TC_vs_dose_archive_figure_6.m uses the MD + TEM defect density calculation data and the kinetic theory model to then give the back calculated expected 
thermal diffusivity reduction. It then plots this with the original TGS measured thermal diffusivity vs dose data, to give FIGURE 6.
