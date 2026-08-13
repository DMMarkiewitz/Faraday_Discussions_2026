Collected here are the critical files to reproduce the results displayed in "Thermodynamic effects of solid electrolyte interphase formation from solvation and ionic association in water-in-salt electrolytes" by Daniel M Markiewitz, Michael McEldrew, Conor ME Phelan, Qianlu Zheng, Jasper Singh, Robert S Weatherup, Rosa M Espinosa-Marzal, Martin Z Bazant, Zachary AH Goodwin (https://arxiv.org/abs/2602.23875).

---

The files are set up to be run directly by the user, with the meshes and data files required in the working directory. The molecular dynamics (MD) simulation data needed for the simulation analysis scripts and pre-built meshes for the theory scripts can be regenerated individually, but for accessibility, they are available at https://doi.org/10.5281/zenodo.20450228. These larger files should be placed in the Runners folder for the scripts to run and create the key figures and results shown in the manuscript.

The files in the LAMMPS folder are those needed to run the MD simulations to generate the data used. This data is also available at https://doi.org/10.5281/zenodo.20450228. The protocol behind these MD simulations has been previously outlined in "Theory of the double layer in water-in-salt electrolytes" by Michael McEldrew, Zachary AH Goodwin, Alexei A Kornyshev, Martin Z Bazant (https://arxiv.org/abs/1808.06118).

In the Active folder, there are the scripts to create the meshes, solve for the theory predictions, and merge plots to create the figures seen in the published manuscript. In order for these scripts to execute, the data stored at https://doi.org/10.5281/zenodo.20450228 must be in the working directory. To prevent overwriting the sample data, we provide the initial set of meshes/data sets needed for the code to run and reproduce the key figures at https://doi.org/10.5281/zenodo.20450228; these files have an S_ at the start, so the default save names do not override the manuscript data sets. Without these manuscript data sets, the scripts do have a hierarchy based on previously generated data. 



For the theory scripts, without the manuscript data, the following scripts should be run first:

- FD_2026_EDL_WISE_STI_Theory (Generates meshes of the key variables' dependence on the non-dimensionalized electrostatic potential and its gradient. This treatment supports the electrical double layer solver as these meshes speed up computation of key quantities such as charge density, non-dimensional dielectric constant, and association probabilities.)

- Theory_Cluster_Distribution_Plot_Bulk (Generates Bulk Cluster Distribution, i.e., concentrations of all clusters up to some cutoff. This script calculates the concentrations of ionic clusters from the association constant, based on the number of anions and cations present. From these concentrations, the distribution is plotted to highlight the clusters' net charge bias and the distribution of large multi-ionic aggregates.)

- Activity_Bulk (Calculates the bulk activity of species as a function of molality. From the functionalities, volume fractions, and the association constant, one can calculate the change in the chemical potential of the solution's species and plot how they change with concentration compared against a reference state. The reference state, as the name suggests, can be chosen; here we chose the 0.5m solution as the reference.)

- Association_Probabilities_Bulk (Calculates the bulk association probability as a function of molality.)



Following executing these scripts for the concentrations of interest, one can run the next set of theory scripts:

- Implicit_Poisson_Solver_Sticky_WiSE_EDL_BVP_Panel (Generates Panels for the theoretical EDL predictions. These panels show the species volume fractions, key cluster volume fractions, association probabilities, and the product of the ionic association probabilities as a function of distance from the electrode. To achieve these predictions, the non-dimensionalized version of the modified Poisson equation from the manuscript is numerically solved with a BVP solver using the previously generated meshes to speed up the calculations.)
	
- Implicit_Poisson_Solver_Sti_WiSE_EDL_BVP_Act_Calc (Generates Plots of the activity near an electrode and meshes. The activity is calculated using the equations shown in the manuscript from the output of the electrical double layer solver. Additionally, the data for these curves are saved so predictions at various concentrations can be merged into a single figure as was done with the activity of the species within the electrical double layer plot in the manuscript.)
	
- Theory_Cluster_Distribution_Plot_Positive_Surface (Generates Cluster Distribution Plots for the positive electrode. This code is similar to the bulk script; here one works with the electrical double layer quantities.)
	
- Theory_Cluster_Distribution_Plot_Negative_Surface (Generates Cluster Distribution Plots for the negative electrode. This code is analogous to the positive electrode script.)



Lastly, with the output of these scripts, one can run the last theory script:

- Activity_Merge_Plotter (Merges the EDL activity plots into a composite figure.)



Executing these scripts will reproduce the key theory figures and findings displayed in our manuscript.



For the MD simulations analysis, the scripts can be executed independently. The main data set needed comes from the readDumpFile function, which processes the MD simulations .xyz file to condense the data into the essential elements for our analysis. We included the necessary ones for the 0.2 C/m^2 analysis, but the other surface charges can be generated by uncommenting the readDumpFile function in either MD processing script.

For the MD simulations analysis:

- WiSE_EDL_MD_Processing_Panels (Generates the panels and the association constants. These panels show species volume fractions, key cluster volume fractions, association probabilities, and the product of the ionic association probabilities as a function of distance from the electrode. The association constants are calculated from the association probabilities.)

- WiSE_EDL_MD_Cluster_Dist (Generates the cluster distribution plots by calculating the concentration of individual clusters in the electrical double layer and showing the corresponding distribution at a specified location in the electrical double layer.)

---

For any questions, please reach out to us.
@ dmm385@mit.edu
@ zac.goodwin@materials.ox.ac.uk / zac.goodwin13@gmail.com
@ bazant@mit.edu