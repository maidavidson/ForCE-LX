# ForCE-LX
Public version of ForCE-LX

This folder contains the code for ForCE-LX.	

This is one-line model for the shoreline displacement due to cross-shore/longshore sediment transport, sea level rise and sediment sources and sinks. The model is designed to accommodate complex coastlines including geological features and anthropogenic structures.	

The code is provided in this instance to support a publication in review. Although the repository is publicly available without licence, the model is not yet ready for general use as it lacks the support guidance required to run the model on new sites effectively. The model runs provided only support the publication in review. As such the authors offer no support for this model/code and accept no liability for its application and use currently.	 
This said further support is planned for the future and this site will be updated accordingly when the paper is published, and the support material is ready.

The RunForce.m program initiates the model run. RunForce.m allows the user to open setup files that are located in any of the following locations:	

Projects/Perrnporth/inputData/setupPerranporth.m	

Projects/Analytical/inputData/setupConstantWaveTest_LS.m	

Projects/Analytical/inputData/setupVariableWaveTest_LSXS.m	

Having loaded the setup file, which contains all the user input variables including the names of the data files used to initialise the model, RunForCE.m then runs the model code ForCE.m. 	
Note that the setup file can be opened and edited by the user to change the model parameters.	

Output data files are then deposited in the outputData folders at the same directory level as the inputData directory. 	
