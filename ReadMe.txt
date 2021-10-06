********************************************************************
* Explanations about the data and code for the cross-sensory study *
********************************************************************
In order to run the code you first need to do several little things:

1) Download the zip file that includes all data and code from github. Place it wherever you like and extract it.
   Make sure you have one main directory named 'Cross-sensory-study-main' that includes several directories and code files.
2) Add the 'functions' folder to matlab paths (through 'set path').
3) Make sure you have 'psignifit' on your computer, or download it from the following link (and add to matlab paths):

https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/neuronale-informationsverarbeitung/research/software/psignifit/

psignifit version 2.5.6, Wichmann and Hill (2001a)
********************************************************************

CODE
----

There are two code files:

- MAIN: run this code to perform the following steps:
	1) Analyzing raw data and generating psychometric curves
	2) Priors Model - fitting logistic regression model to extract beta coefficients
	3) Expanded priors model - fitting expanded model to extract betas when the 3 priors are separated
	4) Plot figures 
see further details inside 'MAIN' on each step.

- SIMULATION: run this code to perform the following steps:
	1) Simulating participants' choices (how exactly? details inside 'SIMULATION')
	2) Fitting the expanded priors model (on simulated data)
	3) Plot scatters to see original betas vs. recovered betas
see further details and explanation about how the simulation works inside 'SIMULATION'.
--------------------------------------------------------------------

DIRECTORIES AND THEIR CONTENT
-----------------------------

- functions - includes matlab functions that are needed to run the code. Please add this folder to matlab paths.
- Plots - after running 'MAIN', this folder will contain the following plots:
	1) Summarized PSE difference (figure 3B)
	2) M3 results summary (figure 4) - for every subplot in figure 4 a different plot is produced here
	3) AIC Bayesian model comparisons (figure 5)
	4) BIC Bayesian model comparisons (supplementary figure 1)
- ves-ves (and the other three conditions) - this folder includes several things for this condition:
	1) Excel file named 'All input_ves-ves' that includes information regarding the pariticpants that 
	performed this condition. It also includes information regarding the sessions they performed (like
	(whether it was excluded etc.).
	2) All - after running 'MAIN', this folder will contain a 'Results' file for this condition.
	'Results' file includes all the data of the participants that performed this condition:
		- AnaData: data after fitting psychometric functions. Includes thresholds, biases, amount of trials etc.
		- dataRec: for each participant there is one data record that includes:
			- trialCount: number of trial (from 1 to the amount of trials performed)
			- stimType: stimulus type of each trial (1-ves, 2-vis)
			- coherence: visual coherence level on each trial (always 100%)
			- dir: heading direction for each trial
			- response: participant's response on each trial (1-left, 2-right)
			- isTestTrial: a flag that indicates whether a trial is a test trial (1) or a prior trial (0)
		- sess_* and subj_*: information regarding the participants and their performed sessions. 
		For example: 
			- sess_AGE includes the age of each participant
			- subj_ID includes the ID of the participant
	3) Priors model - After running 'MAIN', this folder will contain two results files of the model fitting:
		1) ModelResults_avg3priors+RMS - model M3 (in which 3 prior stimuli are averaged) results for each participant. 
		   This includes all four beta coefficients and AIC and BIC for each model.
		2) ModelResults_3priors+RMS - model M3E (in which 3 prior stimuli are separated) results for each participant. 
		   This includes all six beta coefficients and AIC and BIC for each model.
	   Prior model includes also the sub-folder 'Simulation': after running the code 'SIMULATION' for this 
	   condition (by default simulation is running for vis-vis), this folder will contain:
		1) Simulated data in regular order and in reverse order (when betas for prior_stimulus1 and prior_stimulus3
		   are switched.
		2) Recovered beta coefficients for when simulated in regular order and for when simulated in reverse order.
		3) Plots (supplementary figures 2 and 3) - for each subplot in these figures a different plot is produced here.
		
	4) Psychometric curves - after running 'MAIN', this folder will contain all the participants' 
	   psychometric curves in 3 formats:
		1) matlab figure
		2) PDF file
		3) tiff file
		the name of the psychometric file is identical to the name of the raw data file.
	5) Raw - this file contains all raw data for the participants in this condition.
********************************************************************
Good luck!

If you have any questions please contact me at: shalom.shir@gmail.com
or my supervisor Dr. Adam Zaidel at: ajzaidel@gmail.com

Shir.