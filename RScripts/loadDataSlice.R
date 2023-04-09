#Loads data into the workspace based on what is requested in functionBooklet.R's @getParameterSpace and @getRequestSpace

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
}


Parameter_Space = getParameterSpace();
Request_Space = getRequestSpace();

#Load work space; change this to the "big boy" when testing the full dataset.
if (!WORKSPACE_LOADED) 
{
  if (LOAD_FROM_WORKSPACE) 
  {
    WORKSPACE_LOADED = loadFromWorkspace(Workspace_Name);
  }
  else
  {
    File_Cores=getFileCores(Request_Space[[1]], Parameter_Space[[1]]);
    
    # Open missing files document and store to data frame ($Missing_Files).
    Missing_Files = if (Tar_File=="NO_TAR") NULL else read.csv(paste(Infile_Prefix, "Missing-Outfiles", Infile_Suffix, sep=""), header=FALSE);
    
    Data_Storage = loadFromFiles(File_Cores, Parameter_Space[[2]], Request_Space, Rep_Num, Infile_Prefix, Infile_Suffix, Tar_File, Missing_Files)
  }
}