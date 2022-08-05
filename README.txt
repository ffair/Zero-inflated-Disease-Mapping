Supplementary information / reproducible research files for the manuscript 
Title: "Spatial Correlated Incidence Modeling with Zero-Inflation"

Authors: Feifei Wang, Haofeng Li, Han Wang, Yang Li
In case of questions or comments please contact yang.li@ruc.edu.cn

Although parallel computing has been used, the simulation and real data analysis are still a time-consuming process.
Therefore, R workspace is provided here for the purpose of checking.

The ZIP model used for comparison needs the R package 'pscl' to complete the model fitting process, and it can be installed
via 'install.packages("pscl")'
The DM and ZBSJ models used for comparison needs the R package 'INLA', and it can be installed via 'install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)'
The details of 'INLA' package can be found in 'https://www.r-inla.org/home'

The code was written/evaluated in R with the following software versions:
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

loaded via a namespace (and not attached):
[1] compiler_4.2.1 tools_4.2.1   


Simulations were run in parallel on a Linux server with software versions:
R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.7 LTS

Matrix products: default
BLAS:   /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.18.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets 
[7] methods   base     

other attached packages:
[1] MASS_7.3-54       pscl_1.5.5        doParallel_1.0.16
[4] iterators_1.0.13  foreach_1.5.1    

loaded via a namespace (and not attached):
[1] compiler_4.0.5   tools_4.0.5      tinytex_0.30    
[4] codetools_0.2-18 xfun_0.27       




This folder contains the following data and files that can be used to reproduce all analysis, tables and figures of the manuscript.
It contains three subfolders containing the following files:

./code/:
	./01 Simulation with Artificial Geographical Map/:
		./01 First Situation/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model, simulate with
			ZDM, ZIP and POR model and generate true values in simulation.

			simulation.R:
			An R script that perform simulation with artificial geographical map in the first situation (section 3.2). The
			resulting tables are stored in list variable 'tab', with each dimension regarding a type of simulation setting. 
			The resulting figures willed be automatic stored in the folder ./results/.
		
		./Second Situation/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model, simulate with
			ZDM, ZIP and POR model and generate true values in simulation.

			simulation.R:
			An R script that perform simulation with artificial geographical map in the second situation (section 3.2). The
			resulting tables are stored in list variable 'tab', with each dimension regarding a type of simulation setting. 
			The resulting figures willed be automatic stored in the folder ./results/.

	./02 Simulation with Virginia Geographical Map/:
		./01 alpha=0.8_beta=-6/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model,
			generate true values in simulation.

			simulation.R:
			An R script that perform simulation with Virginia geographical map (section 3.3) when alpha=0.8
			and beta = -6. The resulting tables are stored in list variable 'tab'. The resulting figures willed be 
			automatic stored in the folder ./results/.

		./02 alpha=0.05_beta=-6/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model,
			generate true values in simulation.

			simulation.R:
			An R script that perform simulation with Virginia geographical map (section 3.3) when alpha=0.05
			and beta = -6. The resulting tables are stored in list variable 'tab'. The resulting figures willed be 
			automatic stored in the folder ./results/.

		./03 alpha=0.8_beta=-8.2/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model,
			generate true values in simulation.

			simulation.R:
			An R script that perform simulation with Virginia geographical map (section 3.3) when alpha=0.8
			and beta = -8.2. The resulting tables are stored in list variable 'tab'. The resulting figures willed be 
			automatic stored in the folder ./results/.

		./04 alpha=0.05_beta=-8.2/:
			function.R:
			An R script that contains functions to generate spatial random effect, MCMC for ZDM model,
			generate true values in simulation.

			simulation.R:
			An R script that perform simulation with Virginia geographical map (section 3.3) when alpha=0.05
			and beta = -8.2. The resulting tables are stored in list variable 'tab'. The resulting figures willed be 
			automatic stored in the folder ./results/.

	./03 Application on the Lyme Disease in Virginia/:
		./Ecoregion 0/:
			function.R:
				An R script that contains fundamental functions for MCMC process of ZDM model.

			Ecoregion 0.R:
				An R script that contains code for all the 5 models to fit on the data of Ecoregion 0 of Virginia (section 4).

		./Ecoregion 1/:
			function.R:
				An R script that contains fundamental functions for MCMC process of ZDM model.

			Ecoregion 1.R:
				An R script that contains code for all the 5 models to fit on the data of Ecoregion 1 of Virginia (section 4).

./data/:
A folder contains Lyme disease data in Virginia (lyme_data.csv), adjacency matrix of artificial geographical map (W.csv),
adjacency matrix of Ecoregion 0 (W0.csv) and adjacency matrix of Ecoregion 1 (W1.csv).

./results/:
	./01 Simulation with Artificial Geographical Map/:
		./01 First Situation/:
		A subfolder that contains R workspace for simulation with artificial geographical map in the first situation, and will
		contain resulting figures of the simulation setting.

			summary_simulation_artificial_first_situation.RData:
			R workspace for simulation with artificial geographical map in the first situation.
	
		./02 Second Situation/:
		A subfolder that contains R workspace for simulation with artificial geographical map in the second situation, and will
		contain resulting figures of the simulation setting.

			summary_simulation_artificial_first_situation.RData:
			R workspace for simulation with artificial geographical map in the second situation.

	./02 Simulation with Virginia Geographical Map/:
	A subfolder that contains R workspace for simulation with Virginia geographical map, and will contain resulting 
	figures of the simulation setting.

		alpha=0.8_beta=-6.RData:
		R workspace for simulation with Virginia geographical map under alpha=0.05 and beta = -6.

		alpha=0.05_beta=-6.RData:
		R workspace for simulation with Virginia geographical map under alpha=0.05 and beta = -6.

		alpha=0.8_beta=-8.2.RData:
		R workspace for simulation with Virginia geographical map under alpha=0.05 and beta = -6.

		alpha=0.05_beta=-8.2.RData:
		R workspace for simulation with Virginia geographical map under alpha=0.05 and beta = -6.

	./03 Application on the Lyme Disease in Virginia/:
	A subfolder contains resulting 	figures of the simulation setting.