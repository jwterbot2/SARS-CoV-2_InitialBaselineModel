#Loads data into the workspace based on a list of models contained within the file named by $Model_ListFile

#Define management constants, input/output stems, and number of replications each model was run.
##TAR ARCHIVE INFILE STILL NEEDS DEFINED
{
  DEBUG = TRUE;
  GARBAGE_COLLECT = FALSE;
  LOAD_FROM_WORKSPACE = FALSE;
  WORKSPACE_LOADED = FALSE;
  
  Infile_Prefix = "sc2-full2_";
  Infile_Suffix = ".output.csv";
  Tar_File = "NO_TAR";
  Outfile_Stem = "sc2-full2_allOSGOutput_1000samp_1000Window_";
  Rep_Num = 5;
  
  Workspace_Name = "sc2-full2_allOSGOutput_in-R-Storage.RData";
  Model_List_File = "figureModels2.txt";
}

Parameter_Space = getParameterSpace();
Request_Space = getRequestSpace();

Model_List <- read.csv(Model_List_File, header=FALSE, na.strings="na")
Model_Cores = Model_List[,1];
Model_Names = lapply(lapply(1:length(Model_List[[2]]), function(i) lapply(Model_List[,c(2:length(Model_List))], "[[", i)), function(i) paste(i, collapse=","))
Model_Cores = setNames(Model_Cores, Model_Names);
if (!WORKSPACE_LOADED) 
{
  if (LOAD_FROM_WORKSPACE) 
  {
    WORKSPACE_LOADED = loadFromWorkspace(Workspace_Name);
  }
  else
  {
    # Open missing files document and store to data frame ($Missing_Files).
    Missing_Files = if (Tar_File=="NO_TAR") NULL else read.csv(paste(Infile_Prefix, "Missing-Outfiles", Infile_Suffix, sep=""), header=FALSE);
    
    Data_Storage = list() 
    Data_Storage = loadFromFiles(Model_Cores, Parameter_Space[[2]], Request_Space, Rep_Num, Infile_Prefix, Infile_Suffix, Tar_File, Missing_Files)
  }
}

