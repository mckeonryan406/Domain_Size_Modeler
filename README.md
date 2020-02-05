## INVERSE DOMAIN SEARCH FOR MDD ANALYSIS

Created By: Ryan McKeon 
Date: 10 March 2017 
Contact: [ryan.e.mckeon@dartmouth.edu](mailto:ryan.e.mckeon@dartmouth.edu)



#### PURPOSE

***General Description:*** This program built in Matlab uses an inverse approach to fit observed data where the inputs controlling the observed data are not expressly known and must be inferred iteratively by the model.  Specifically, this set of tools is directed at geochronological applications and the program was written in conjunction with research during my postdoc at Caltech, however the general idea of using observed results and known behavior of the controls on those results can be applied to a wide range of different problems.  This model iteratively produces inputs, forward models what results those inputs would create, compares that to the observed data, and then uses a Bayesian Information Criterion to determine the best combination of model fit and simplicity in determining the most likely inputs that produced the results observed.    

***Geochronological Nitty-gritty:*** This code generates a diffusion domain structure for Multiple Diffusion Domain analysis from step heating data (either 4He/3He or Ar/Ar). At present, this code uses a SPHEREICAL diffusion domain only and the analytical solution for which is from Farley et al., 2010 G-cubed. This code uses the Bayesian Information Criterion score (BIC) to balance the improvement of a model fit vs. the number of domains used to generate said fit. I.E. is it better to use fewer domains and introduce less complexity to the model? In the end you can decide, but the results here may suprise you.

#### UPDATES

**27 August 2018**: Added physical domain size output and fixed bugs --> Using newly published kinetics for helium diffusion from hematite from Farley 2018 GCA, the code now estimates the physical size of the modeled domains.  Bug fixes addessed issues with final plots and domain size calculation.



#### QUICK START GUIDE

1. Download and unzip *Domain_inverse_search_v2.3_CURRENT.zip*. The unzipped folder contains all of the necessary scripts that are associated with the inverse model.  
2. Use one of the example files from our Michigan Hematite samples published in Farley and McKeon 2015 - Geology that are included with the code.
3. Fire up Domain_Model_MASTER.m and run it as a script!! Hit the green play botton and watch the Command Window for instructions. This is hard wired to start with MI43d2.EXAMPLE.txt -- this is a good place to start.
4. When prompted for kinetics info, enter "1" (no quotes) to regress steps from the Arrhenius plot, then enter "9 11 13 15" (no quotes) to regress these steps to derive the diffusion kinetics.
5. Watch the model run, out of the box it comes set up to attempt to fit the observed data with 2, 3, 4, 6, and 9 domains respectively. This will take 10 minutes or so to run. At the end, you will have a 4 panel plot window appear summarizing the results of the model run.
6. Fire up Final_Results_Plotter.m and set the value of "ndomains" to 6. This will export a pdf with an Arrhenius Plot and Delta (Lovera's r/ro) plot for the 6 domain model run. It will also export a .csv file with all of the modeled and observed data used and produced by the model.

#### HOW THE CODE WORKS

The Domain Modeler reads in observed step heating results (input file defined above and described below) and asks the user to choose how to give the model diffusion kinetics, either by manually entering them in the Command Window when prompted, or by regressing step heating data prior to starting the inverse search. The kinetics data sets the reference domain size to an arbitrary (but mathmatically simple) value of 1 and assigns it the user specified (through regression or manual typing) Activation Energy (Ea in kcal/mol) and diffusity (Ln(Do/a2). The code then iteratively changes the size parameter of the domains it creates to try to fit the observed data. Chaning the size changes the diffusivity of each model domain relative to the reference from the kinetics, the Ea is held constant for ALL domains. The code also adjusts the gas content of each of the modeled domains, splitting the 3He yield between varying sizes of domains to fit the observed step heating results. This serach uses the r/ro plot of Lovera et al. (1991) JGR -- Here referred to as delta --
to optimize the modeled data with the observed results. The search is optimized by fitting both the step gas release (F) and the cumulative gas release (Fcum) of the model with the observed results (after Farley and Flowers, 2012 EPSL).

This search is accomplished in two steps, the first is purely random, where the domains are assigned a random size and proportion of the total gas with no attempt to hone the fit. The section step takes the best fit from the random attempts and hones the size and gas fraction of each domain in a systematic fashion. This code is set up to iterate through attempting to fit the data using a different number of domains and for each specific number of domains it tries multiple times from scratch to fit the data in order to avoid hammering away for a long time on a local minimum of mistfit or BIC score.

This code uses the Bayesian Information Criterion (BIC) (calcualted in master_domain.m) to determine the best fitting model. As indicated above this metric favors simplier models (i.e. fewer domains) and at the conclusion of the code, the output is given for the domain structure that produces the lowest BIC score even if it is visually not as good a fit as other models with more domains. You can plot and export the results for model runs with more domains and a higher BIC (and likely a lower misfit score) by using the Final_Results_Plotter.m script which calls upon Master_Results_Table 3D array that stores all of the output from the model. The format for INPUTS and OUTPUTS are explained next.

#### Inputs and Outputs:

Inputs: obs_data_file = observed data array - Defined Above! 
- col 1 = step # 
- col 2 = Temp (deg C) 
- col 3 = time (hr) 
- col 4 = Fcum observed

Outpus: generates an Excel file called "modeled_data_OUT.xls", which stores the results from the model iteration with the lowest BIC score for the whole modeling session. I.E. if you try a wide range of differnt numbers of domains, only one single combination of goodness of fit with allowed complexity will be saved. This file has the following structure:

- row 1 = misfit score, BIC (bayesian information criterion), ndomains, Ea, lnDo/a2, slope 
- row 2 = Gas fraction of each domain modeled 
- row 3 = Relative domain size for each domain modeled 
- row 4 = log(Do/a2) for each domain --> for input to QtQT 
- row 5 = closure temperature for each domain (assuming a 10 degC/Ma cooling rate) 
- row 6 - n:
>> - col 1 = step number 
>> - col 2 = Temp (deg C) 
>> - col 3 = time (hr) 
>> - col 4 = Fcum modeled 
>> - col 5 = 10000/K 
>> - col 6 = lnDa2 modeled 
>> - col 7 = delta modeled --> ln(r/ro) in Lovera-speak 
>> - col 8 = delta OBSERVED --> included for plotting comparison 
>> - col 9 = lnDa2 OBSERVED --> included for plotting 
>> - col 10 = Fcum OBSERVED --> included for plotting

Also generates plots -- The running Delta plot shows the progress of each iteration as the fit is improved. At the conclusion of the code, a pdf with plots of Delta (model vs. observed), arrhenius data (observed with modeled and reference domains of 100 nm, 10 um, and 1 mm radius) and a comparison of best BIC scores for different numbers of domains are saved as "domain_model_plots.pdf".

#### IMPORTANT 

Running this code as a script (by hitting the green play botton above rather than as a function from the command window) saves the in-run variables to the Workspace. This is recommended for this code because it allows you to quickly plot up the final results for any specific number of domains from the previous model run using Final_Results_Plotter.m instead of just getting the iteration with the lowest BIC score as this code is hardwired to do.
