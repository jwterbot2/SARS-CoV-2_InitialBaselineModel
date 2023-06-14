#import required packages
{
  library(ggplot2);
  library(gridExtra);
}

#function booklet

#@thankyou Documentation
#Requires nothing, returns nothing; prints "You're Welcome!").
thankyou <- function() print("You're Welcome!");

#@getFileCores Documentation
{
  #$models :  A list of Boolean lists. Each Boolean list describes the parameters sets used in that model.
  #$param_set : A list of parameter sets as Booleans; each set is one or more lists of parameter values.
  #             For example, (Required Parameters (4), Recombination (1), Progeny Skew (2), and 
  #             Distribution of Fitness Effects (1).
  #The length of each list in $models should be the same as the length of $param_set. 
  #<-Returns:
  #$file_cores :  A list of "file cores" which are the inner stem of the infile names composed of period
  #               separated numbers which encode the parameter values used. The names of the elements in
  #               this list describe the parameters' values in detail.
}
getFileCores <- function(models, param_set)
{
  getSubsetCodes <- function(param_subset, subset_mult)
  {
    code_list = c();
    if(length(param_subset)>1)
    {
      partial_codes = getSubsetCodes(param_subset[-1], subset_mult[-1]);
      if(partial_codes[1]<0)
      {
        return(setNames(c(-1),paste("na",names(partial_codes)[1],sep=",")));
      }
      for(i in 1:length(param_subset[[1]]))
      {
        if(param_subset[[1]][[i]])
        {
          for (j in 1:length(partial_codes))
          {
            code_list = append(code_list, ((i-1)*subset_mult[[1]] + partial_codes[j]));
            names(code_list)[length(code_list)] = paste(names(param_subset[[1]][i]), ",", names(partial_codes)[j], sep="");
          }
        }
      }
      return(code_list);
    }
    else
    {
      for(i in 1:length(param_subset[[1]]))
      {
        if(param_subset[[1]][[i]])
        {
          code_list = append(code_list, (i-1)*subset_mult[[1]]);
          names(code_list)[length(code_list)] = names(param_subset[[1]][i]);
        }
      }
      if(length(code_list)!=0)
      {
        return(code_list);
      }
      return(c("na"=-1));
    }
  }
  getParamSetCodes <- function(model, param_codes)
  {
    if(is.null(names(model)))
    {
      modelName = "-";
    }
    else
    {
      modelName = names(model)[1];
      model = model[[1]];
    }
    param_set_codes = list();
    for (i in length(model):1)
    {
      cur_codes = list();
      if (!model[i])
      {
        n_subparams = gregexpr(",",names(param_codes[[i]])[1])[[1]];
        if(n_subparams[1]<0)
        {
          n_subparams = 0;
        }
        else
        {
          n_subparams = length(n_subparams);
        }
        n_subparams = n_subparams+1;
        cur_codes = setNames(-1, paste(replicate(n_subparams, "na"), collapse=","));
      }
      else
      {
        cur_codes = param_codes[[i]];
      }
      
      if (length(param_set_codes) == 0)
      {
        param_set_codes = cur_codes;
        names(param_set_codes) = names(cur_codes);
      }
      else
      {
        new_codes = list();
        for (j in 1:length(cur_codes))
        {
          for(k in 1:length(param_set_codes))
          {
            new_codes = append(new_codes, paste(cur_codes[j],".", param_set_codes[k], sep=""));
            names(new_codes)[length(new_codes)] = paste(names(cur_codes)[j], names(param_set_codes)[k],sep=",");
          }
        }
        param_set_codes = new_codes;
      }
    }
    param_set_codes = setNames(paste(modelName, ".", param_set_codes, sep=""), names(param_set_codes));
    return(param_set_codes);
  }
  
  param_codes = list();
  for (subset in param_set)
  {
    subset_lens=list();
    subset_mults=list(1);
    for(i in 1:length(subset))
    {
      if(i!=1)
      {
        subset_mults = append(subset_mults, (subset_mults[[i-1]]*subset_lens[[i-1]]));
      }
      subset_lens = append(subset_lens, length(subset[[i]]));
    }
    param_codes = append(param_codes, list(sort(getSubsetCodes(subset, subset_mults))));
  }
  
  file_cores = list();
  for (i in 1:length(models))
  {
    file_cores = append(file_cores, getParamSetCodes(models[i], param_codes));
  }
  
  return(file_cores);
}

#@makeStorageBranch Documentation
{
  #$file_core : A string detailing the center stem of the infile and which parameter values were used for
  #             that simulation (found in name).
  #$param_names : A list of the names of the parameters used in the models.
  #$sample_sizes :  A list of strings of which sample sizes to use. The names of each element should 
  #                 correspond to any appendix added to the infile.
  #$window_sizes :  A list of strings of which window sizes to use. Currently all windows for a $file_core
  #                 and window size combination should be in the same infile. Future implementations may
  #                 change this and add window size information to file name, but this function will need
  #                 revising to allow this.
  #$n_reps :  Number of replicates run for each $file_core.
  #$stats_to_get :  A list of strings of which statistics to get from the infile.
  #$infile_prefix : A string that all infiles should begin with.
  #$infile_suffix : A string that all infiles should end with.
  #$infile_archive : Option not fully implemented; this would be a compressed archive containing the infiles. Currently this should be "NO_TAR"
  #$missing_files : Option not fully implemented; a file or list of which files in the set are not contained in the folder or archive.
  #Creates an element of the data storage, $new_branch, populates it with contextual information, and data.
  #<-Returns:
  #$new_branch :  A multidimensional list for a set of model parameter that can be placed into a list with other branches. Contained in each dimension is:
  #					[1] : Model core (name) and...
  #					[2] : List of parameters and...
  #					[3] : Replicates containing...
  #					[4] : Sample Sizes containing...
  #					[5] : File Name and...
  #					[6] : Window Sizes containing...
  #					[7] : Statistics
}
makeStorageBranch <- function(file_core, param_names, sample_sizes, window_sizes, n_reps, stats_to_get,
                              infile_prefix, infile_suffix, infile_archive=NULL, missing_files=NULL)
{
  #We will be using these values multiple times.
  n_sample_sizes = length(sample_sizes);
  n_window_sizes = length(window_sizes);
  n_stats_to_get = length(stats_to_get);
  
  #Creates the skeleton for the data storage branch named $file_core. Initializes with the values and names for the parameters from $file_core.
  #Currently does not get the parameter values from the file_core if they are not provided in the file_core's name.
  if (!is.null(names(file_core)))
  {
    new_branch = setNames(list(list(setNames(strsplit(names(file_core)[1], split=",")[[1]], param_names), replicate(n=n_reps, expr=list(setNames(replicate(n=n_sample_sizes, expr=list(list("", setNames(replicate(n=n_window_sizes, expr=list(setNames(replicate(n=n_stats_to_get, expr=list()), stats_to_get))), names(window_sizes))))), names(sample_sizes)))))), file_core);
  }
  else
  {
    new_branch = setNames(list(list(param_names, replicate(n=n_reps, expr=list(setNames(replicate(n=n_sample_sizes, expr=list(list("", setNames(replicate(n=n_window_sizes, expr=list(setNames(replicate(n=n_stats_to_get, expr=list()), stats_to_get))), names(window_sizes))))), names(sample_sizes)))))), file_core);
  }
  
  
  n_sample_sizes = length(sample_sizes);
  n_window_sizes = length(window_sizes);
  n_stats_to_get = length(stats_to_get);
  
  #For each rep of this model core
  for (rep in 1:n_reps)
  {
    #For each sample size being analyzed
    for (i_sample_size in 1:n_sample_sizes)
    {
      #The name of the infile for this combination containing the stats.
      cur_filename = paste(infile_prefix, file_core, ".", (rep-1), names(sample_sizes)[i_sample_size], infile_suffix, sep="");
      if (DEBUG) print(cur_filename);
      #Check if the file is missing versus a list of missing infiles and their reasons (must have 2 columns);
      file_missing = filterFile(cur_filename, missing_files, infile_archive);
      if(file_missing[1] >= 0)
      {
        #If the file was missing, trim this branch at the 5th dimension with the data frame with two elements, "na" and the reason the file is missing.
        new_branch[[1]][[2]][[rep]][[i_sample_size]] = file_missing
      }
      else
      {
        #If the file is not missing, set the corresponding element as the $cur_filename to store the filename.
        new_branch[[1]][[2]][[rep]][[i_sample_size]][[1]] = cur_filename;
        if (infile_archive != "NO_TAR") untar(infile_archive, cur_filename, verbose = TRUE);
        cur_data = read.csv(cur_filename);
        if (infile_archive != "NO_TAR") file.remove(cur_filename);
        
        for(i_window_size in 1:n_window_sizes)
        {
          for(i_stat_to_get in 1:n_stats_to_get)
          {
            #Check if stat 9, "nUniqHapl" is checked only alongside WindowSize[1]==T or error.
            if (i_stat_to_get==9 & i_window_size > 1)
            {
              stat_data_frame = data.frame("na", reason="Stat not available for windows.");
            }
            else if (i_stat_to_get >= 15)
            {
              if(i_window_size > 1)
              {
                stat_data_frame = data.frame("na", reason="Stat not available for windows.");
              }
              else
              {
                stat_data_frame =  cur_data[(cur_data$size==window_sizes[i_window_size] & cur_data$stat==stats_to_get[i_stat_to_get]), c((ncol(cur_data)-4), (ncol(cur_data)-3), ncol(cur_data))];
              }
            }
            else
            {
              stat_data_frame =  cur_data[(cur_data$size==window_sizes[i_window_size] & cur_data$stat==stats_to_get[i_stat_to_get]), c((ncol(cur_data)-4), (ncol(cur_data)-3), ncol(cur_data))];
            }
            
            new_branch[[1]][[2]][[rep]][[i_sample_size]][[2]][[i_window_size]][[i_stat_to_get]] =  stat_data_frame;
          }
        }
      }
    }
  }
  
  return(new_branch);
}

#@filterFile Documentation
{
  #$file_to_get : A string naming the file to check against the missing file list.
  #$files_missing : A data frame with the first column serving as a list of missing files and the second 
  #                 column the cause of the missing file.
  #$infile_archive : A filename corresponding to where infiles can be found if using a compressed archive.
  #Checks the file against all of the missing files;
  #<-Returns:
  #$is_missing :  If the file is missing, this is a data frame with two elements: "na" and the reason the 
  #               file is missing. If the file is not missing, it returns -1.
}
filterFile <- function(file_to_get, missing_files=NULL, infile_archive=NULL)
{
  #If the file is not found in the list then the default (-1) is returned to communicate the file exists.
  is_missing = -1;
  if ((is.null(infile_archive) | infile_archive == "NO_TAR") & !file.exists(file_to_get))
  {
    is_missing = data.frame(filename="na", reason="File not in directory already.");
  }
  else if (is.null(missing_files) & infile_archive != "NO_TAR")
  {
    ##IMPLEMENT THIS OPTION TO CHECK TAR FILE FOR FILE
    is_missing = data.frame(filename="na", reason="Cannot use tar file without missing file document");
  }
  else
  {
    #For each row in the $missing_files data frame...
    for (row in rownames(missing_files))
    {
      #Check if the file in that row is equal to $file_to_get
      if (file_to_get == missing_files[row,][1]) 
      {
        #If so, set $is_missing to the two element data frame described above and break from loop.
        is_missing = data.frame(filename="na", reason=missing_files[row,2]);
        break;
      }
    }
  }
  return(is_missing);
}

#@loadFromWorkspace Documentation
{
  #$filename : The filename with the workspace data
  #Simple, redundant function that loads a workspace file.
  #<-Returns:
  #TRUE
}
loadFromWorkspace <- function(filename)
{
  load(filename);
  return(TRUE);
}

#@loadFromFiles Documentation
{
  #$file_cores : A list of file cores (aka model cores) that also are the middle of the infile names.
  #$param_names : A list of the names of the parameters used in the models.
  #$request_space : A list of lists; a list of model parameters being examined, a list of which sample sizes to get, a list of which window sizes to get, and a list of which stats to get.
  #$n_reps : Number of replicates each file core (aka model core) should have.
  #$infile_prefix : A string that all infiles should begin with.
  #$infile_suffix : A string that all infiles should end with.
  #$infile_archive : Option not fully implemented; this would be a compressed archive containing the infiles. Currently this should be "NO_TAR"
  #$missing_files : Option not fully implemented; a file or list of which files in the set are not contained in the folder or archive.
  #Loads requested data for each file core (aka model core) in $file_cores.
  #<-Returns:
  #Data_Storage : A list of storage branches (multidimensional lists) that contain the requested data
}
loadFromFiles <- function(file_cores, param_names, request_space, n_reps, infile_prefix, infile_suffix, infile_archive, missing_files)
{
  Data_Storage = list();
  
  #For loop is done using iterator in order to preserve and pass along name information of the model cores which contains the parameter level information
  for ( i in 1:length(file_cores))
  {
    Cur_Branch = makeStorageBranch(file_cores[i], param_names, request_space[[2]], request_space[[3]], n_reps, request_space[[4]], infile_prefix, infile_suffix, infile_archive, missing_files);
    Data_Storage = append(Data_Storage, Cur_Branch);
    
  }
  return(Data_Storage);
}

#@getParameterSpace Documentation
{
  #Define parameters in each parameter set and then assembles the master parameter set, $Params.
  #Each parameter set consists of Boolean lists named for a parameter. The Boolean values describes
  #if that parameter value (the name of that element) is being used. Changes to the requested models
  #need to be made by altering this function.
  #<-Returns:
  #Params : List of lists of which parameter levels to get
  #Param_Names : List of parameter names.
}
getParameterSpace <- function()
{
  ReqdParams = setNames(list(setNames(list(setNames(list(1,1,1), list("Lo", "Md", "Hi")),
                                           setNames(list(1,1,1), list("Lo", "Md", "Hi")),
                                           setNames(list(1,1,1), list("Lo", "Md", "Hi")),
                                           setNames(list(1,1,1,1), list("Lo", "Md", "Hi", "VLo"))),
                                      list("mu", "K", "init", "runtime"))), "ReqdParams");
  
  RecombParams = setNames(list(setNames(list(setNames(list(1,1,1), list("Lo", "Md", "Hi"))),
                                        list("r"))), "RecombParams");
  
  
  ProgSkewParams = setNames(list(setNames(list(setNames(list(1,1,1), list("Lo", "Md", "Hi")),
                                               setNames(list(1,1,1), list("Lo", "Md", "Hi"))),
                                          list("xi", "burstN"))), "ProgSkewParams");
  
  DFEParams = setNames(list(setNames(list(setNames(list(1,1,1), list("Lo", "Md", "Hi"))),
                                     list("dfeDist"))), "DFEParams");
  
  Params = c(ReqdParams, RecombParams, ProgSkewParams, DFEParams);
  Param_Names = list();
  for ( Param_Set in Params)
  {
    Param_Names = append(Param_Names, names(Param_Set));
  }
  
  return(list(Params, Param_Names));
}

#@getRequestSpace Documentation
{
  #Define which models(name = in file core), sample sizes (name = file appendix), window sizes, and statistics to get.
  #Note the Boolean list from @as.logical() dictates which are being analyzed. Changes to the requested data sets
  #need to be made by altering this function
  #<-Returns:
  #List of lists with the models to get, the sample sizes to get, the window sizes to get, and the statistics to get
}
getRequestSpace <- function()
{
  #Note that Models_To_Get is a list of Boolean lists, each Boolean list controls which parameter sets are inluded in the model.
  #The length of each Boolean list should be the same as the length of $Params.
  Models_To_Get = setNames(list(c(1,0,0,0), c(1,1,0,0), c(1,1,1,0), c(1,1,1,1)), 
                           list("8","9","11","15"))[
                             as.logical(list(1,1,1,1))];
  
  SampleSizes_To_Get = setNames(list("100", "1000"), 
                                list("_partial_100", "_partial_1000"))[
                                  as.logical(list(1,1))];
  
  WindowSizes_To_Get = setNames(list("30000", "100", "1000", "2000"), 
                                list("Full", "100", "1000", "2000"))[
                                  as.logical(list(1,0,1,0))];
  
  Stats_To_Get = list("fixBurnMuts", "fixNewMuts", "fixSampMuts",  
                      "thetaPi", "thetaW", "thetaH", 
                      "fayWuH", "tajimasD",
                      "nUniqHapl", "haplDiv", 
                      "gH1", "gH12", "gH123", "gH2/H1",
                      "totalVars", "fixFiltMuts", 
                      "errFiltMuts", "filtSNPs", "filtUniqHapl")[
                        as.logical(list(0,0,0,
                                        1,0,0,
                                        0,0,
                                        0,0,
                                        0,0,0,0,
                                        1,0,
                                        0,1,0))];
  return(list(Models_To_Get, SampleSizes_To_Get, WindowSizes_To_Get, Stats_To_Get));
}

#@drawBranchFigures Documentation
{
  #$data_branch : A multidimensional list with data for a model.
  #$samplesizes_to_get : List of sample sizes to draw figures for
  #$windowsizes_to_get : List of which window sizes to draw figures for
  #$stats_to_get : List of statistics to draw figures of
  #$figure_title : Optional string of the figure title; otherwise defaults to the name of $data_branch
  #$filter_models : Boolean value for if a filtering criteria should be applied
  #Partially written function to draw figures of all combinations of sample size, window size, and statistic requested.
  #Usage not recommended.
}
drawBranchFigures <- function(data_branch, samplesizes_to_get, windowsizes_to_get, stats_to_get, figure_title=NULL, filter_models=FALSE)
{
  #figures = list();
  for(i_samplesize in 1:length(samplesizes_to_get))
  {
    for(i_windowsize in 1:length(windowsizes_to_get))
    {
      stats_data = NULL;
      for(i_stats_to_get in 1:length(stats_to_get))
      {
        stats_data = getStatForBranch(data_branch, samplesizes_to_get[i_samplesize], windowsizes_to_get[i_windowsize], stats_to_get[[i_stats_to_get]]);
        
        if (!is.null(stats_data))
        {
          if(filter_models)
          {
            title = if (is.null(figure_title)) names(data_branch)[[1]] else paste(strsplit(names(data_branch)[[1]], "\\.")[[1]][[1]], figure_title,sep=",");
            subtitle=paste("Window Size:", windowsizes_to_get[[i_windowsize]], "| Sample Size:", samplesizes_to_get[[i_samplesize]], sep=" ");
            # figures = append(figures, plot(stats_data[[1]], stats_data[[2]], main=names(data_branch)[[1]], sub=subtitle,  xlab="Genome Position", ylab=stats_to_get[i_stats_to_get]));
            plot(stats_data[[2]][stats_data$filter==0.02], stats_data[[3]][stats_data$filter==0.02], main=title, sub=subtitle,  xlab="Genome Position", ylab=stats_to_get[i_stats_to_get]);
            #  #FILTERING LOGIC TO REWRITE TO CONSIDER FILTERED SNPS
            #  filter_max = 5;
            # # print(stats_data[[3]]);
            mean_stat = mean(as.numeric(stats_data[[3]][stats_data$filter==0.02]));
            print(mean_stat);
            # # if(mean_stat > filter_max) {draw_figure=FALSE};
            #  mean_stat = mean(as.numeric(stats_data[[3]][stats_data$filter==0.02]));
            # # print(mean_stat);
            #  if(mean_stat < filter_max & mean_stat > 0) 
            #  {
            #    #print("POSSIBLE!")\
            #    print(mean_stat);
            #  }
          }
          else
          {
            title = if (is.null(figure_title)) names(data_branch)[[1]] else paste(strsplit(names(data_branch)[[1]], "\\.")[[1]][[1]], figure_title,sep=",");
            subtitle=paste("Window Size:", windowsizes_to_get[[i_windowsize]], "| Sample Size:", samplesizes_to_get[[i_samplesize]], sep=" ");
            # figures = append(figures, plot(stats_data[[1]], stats_data[[2]], main=names(data_branch)[[1]], sub=subtitle,  xlab="Genome Position", ylab=stats_to_get[i_stats_to_get]));
            plot(stats_data[[2]], stats_data[[3]], main=title, sub=subtitle,  xlab="Genome Position", ylab=stats_to_get[i_stats_to_get]);
          }
        }
        else if (DEBUG)
        {
          print(paste("Model core: ", names(data_branch)[[1]], " does not have any data for this combination of sample size, window size, and statistic"));
        }
      }
    }
  }
  return(FALSE);
}

#@getStatForBranch Documentation
{
  #$data_branch : A multidimensional list with data for a model
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$stat_to_get : A single stat for which to get the data
  #Gets the requested data from the $data_branch of the requested sample size, window size, and statistic.
  #<-Returns:
  #$stats_data : The requested data as a matrix
}
getStatForBranch <- function(data_branch, samplesize_to_get, windowsize_to_get, stat_to_get)
{
  stats_data = NULL;
  for(i_rep in 1:length(data_branch[[1]][[2]]))
  {
    cur_filename = data_branch[[1]][[2]][[i_rep]][[names(samplesize_to_get)]][[1]]
    if (cur_filename != "na")
    {
      z = data_branch[[1]][[2]][[i_rep]][[names(samplesize_to_get)]][[2]][[names(windowsize_to_get)]][[stat_to_get]];
      stats_data = if (is.null(stats_data)) z else rbind(stats_data, z);
    }
    else if (DEBUG)
    {
      writeLines(paste("At rep:", (i_rep-1), "\nInfile was missing\nReason given is:\n", data_branch[[1]][[2]][[i_rep]][[names(samplesize_to_get)]][[2]]));
    }
  }
  if(!is.null(stats_data)) stats_data$value = as.numeric(stats_data$value);
  return(stats_data);
}

#@filterByStat Documentation
{
  #$data_branch : A multidimensional list with data for a model
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$stat_to_get : A single stat for which to get the data
  #$variant_filter_value : Value of minimum frequency for variants to be used (defaults to no filter, 0)
  #$min : Minimum mean value for the requested statistic (defaults to 0)
  #$max : Maximum mean value for the requested statistic (defaults to 1)
  #$useMinMax : List of two Boolean values detailing if $min and $max should be used (defaults to TRUE for both)
  #$inclusive : List of two Boolean values detailing if the $min and $max are inclusive or not (defaults to exclusive $min and $max)
  #Checks if a models mean value for a requested statistic falls within a requested threshold as defined by $min, $max, $useMinMax, and $inclusive.
  #<-Returns:
  #A Boolean value of if the model passes the threshold criteria, or FALSE if the model has no data for the requested statistic, or TRUE if neither $min nor $max are used.
}
filterByStat <- function(data_branch, samplesize_to_get, windowsize_to_get, stat_to_get, variant_filter_value=0, min=0, max=1, useMinMax= c(TRUE,TRUE), inclusive=c(FALSE,FALSE))
{
  if (useMinMax[1] || useMinMax[2])
  {
    stat_data = getStatForBranch(data_branch, samplesize_to_get, windowsize_to_get, stat_to_get);
    if(is.null(stat_data))
    {
      if (DEBUG) writeLines("No data available for this combination of sample size, window size, filtering frequency, and statistic; returning 'FALSE'");
      return(FALSE);
    }
    stat_mean_value = mean(as.numeric(subset(stat_data, filter==variant_filter_value)$value));
    return((useMinMax[1] && (stat_mean_value > min || (inclusive[1] && stat_mean_value == min))) && (useMinMax[2] && (stat_mean_value < max || (inclusive[2] && stat_mean_value == max)))); 
  }
  if (DEBUG) writeLines("No filtering criteria selected, returning 'TRUE'.");
  return(TRUE);
}

#@filterByStat Documentation
{
  #$data_branch : A multidimensional list with data for a model
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$filter_value : Value of minimum frequency for variants to be used for the SNP count graph
  #$snp_y_max : Max value for the y axis of the SNP count graph (defaults to 10)
  #$tpi_y_max : Max value for the y axis of the theta pi graph (defautls to 0.001)
  #Generates a figure with two panels. The left panel is a jitter plot with the number of SNPs in each replicate of the requested model using the filtered data as described by $filter_value.
  #The right panel is of a jitter plot of the values of theta pi using the full data set; the x-axis of this panel corresponds to the position along the genome
  #<-Returns:
  #$combinedFig : The two panel figure
 }
getBranchDiversityFigure_v1 <- function(data_branch, samplesize_to_get, windowsize_to_get, filter_value, snp_y_max=10, tpi_y_max=0.001)
{
  filtered_snps = subset(getStatForBranch(data_branch, samplesize_to_get, setNames("30000", "Full"), "filtSNPs"), filter==filter_value);
  theta_pi = getStatForBranch(data_branch, samplesize_to_get, windowsize_to_get, "thetaPi")
  
  title_Size = 12
  axis_Title_Size = title_Size - 1;
  axis_TickLabel_Size = title_Size - 3;
  point_Alpha = 0.75;
  point_Size = 0.8;
  line_Size = 0.8;
  
  filtered_snps_fig = ggplot(filtered_snps, aes(x = start, y = value)) + geom_jitter(width = 0.1, height = 0) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_text(size=axis_Title_Size), axis.text.y=element_text(size=axis_TickLabel_Size), axis.title.x=element_blank(), title = element_text(size = title_Size)) +
    theme(plot.margin=unit(c(0.3,0.3,1.3,0.4), "cm")) + 
    scale_x_continuous(breaks = NULL, limits = c(0.85, 1.15)) + scale_y_continuous(breaks = seq(0, snp_y_max, by = 2), limits = c(0,snp_y_max)) + labs(y="Number of Segregating SNPs");
  theta_pi_fig = ggplot(theta_pi, aes(x = start, y = value)) + geom_jitter(width = 50, height = 0) + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=axis_TickLabel_Size), axis.title=element_text(size=axis_Title_Size), title = element_text(size = title_Size)) +
    ylim(0, tpi_y_max) + labs(y="Theta-Pi", x="Start Site of Window") +
    theme(plot.margin=unit(c(0.3, 0.3, 0.3, 0.4), "cm"));
  
  combinedFig = grid.arrange(grobs = list(filtered_snps_fig, theta_pi_fig), ncol=2, widths = c(1,4));
  
  return(combinedFig);
}

#@getRequiredPlot_v1 Documentation
{
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$stat_to_get : Which statistic (from array in @getRequestSpace) is being used for this figure
  #$mu_lvl : The (zero-) index of the mu value to use based on array in @getParameterSpace
  #$K_lvl : The (zero-) index of the K value to use based on array in @getParameterSpace
  #$init_lvl : The (zero-) index of the init value to use based on array in @getParameterSpace
  #Generates plots for all models of a given set of required parameters (but all runtime values), partly hardcoded with values specific to the original study regarding the vertical lines that seperate the runtime categories
  #<-Returns:
  #$figure : The requested figure
}
getRequiredPlot_v1 <- function(samplesize_to_get, windowsize_to_get, stat_to_get, mu_lvl=0, K_lvl=0, init_lvl=0, filt_val = 0.02, ymax = NULL, plot_data=NULL)
{
  if(is.null(plot_data))
  {
    plot_data = getRequiredPlotDataFrame(samplesize_to_get, windowsize_to_get, stat_to_get, mu_lvl, K_lvl, init_lvl, filt_val);
  }
  if(is.null(ymax))
  {
    ymax = max(plot_data$Max) + 10
  }
  figure = ggplot(data=plot_data) + geom_segment(aes(x = x, xend = x, y = Min, yend = Max)) + geom_point(aes(x=x, y=Mean), size = 0.5) + 
    geom_vline(xintercept=c(117, 239, 361), color = "grey50", linetype = "dashed") + geom_hline(yintercept=5, color = "blue") + 
    theme_bw() + labs(y = "SNPs") +  ylim(0,ymax) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid = element_blank());
  return(figure)
}

#@getRequiredPlot_v1_wColor Documentation
{
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$stat_to_get : Which statistic (from array in @getRequestSpace) is being used for this figure
  #$mu_lvl : The (zero-) index of the mu value to use based on array in @getParameterSpace
  #$K_lvl : The (zero-) index of the K value to use based on array in @getParameterSpace
  #$init_lvl : The (zero-) index of the init value to use based on array in @getParameterSpace
  #Generates same plot as @getRequiredPlot_v1 but with line segments colored based on model complexity
}
getRequiredPlot_v1_wColor <- function(samplesize_to_get, windowsize_to_get, stat_to_get, mu_lvl=0, K_lvl=0, init_lvl=0, filt_val = 0.02, ymax = NULL, plot_data=NULL)
{
  if(is.null(plot_data))
  {
    plot_data = getRequiredPlotDataFrame(samplesize_to_get, windowsize_to_get, stat_to_get, mu_lvl, K_lvl, init_lvl, filt_val);
  }
  if(is.null(ymax))
  {
    ymax = max(plot_data$Max) + 10
  }
  #Add column with hard-coded values of model complexity (i.e. only works with the specific model set used in the original study)
  plot_data["compl"] = c(rep("Bottleneck",4), rep("Recombination",12), rep("ProgSkew", 108), rep("DFE", 324))
  ggplot(data=plot_data) + geom_segment(aes(x = x, xend = x, y = Min, yend = Max, color = compl)) + geom_point(aes(x=x, y=Mean), size = 0.5) + 
    geom_vline(xintercept=c(117, 239, 361), color = "grey50", linetype = "dashed") + geom_hline(yintercept=5, color = "blue") + 
    theme_bw() + labs(y = "SNPs") +  ylim(0,ymax) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid = element_blank());
}

#@getRequiredPlotDataFrame Documentation
{
  #$samplesize_to_get : A single sample size for which to get the data
  #$windowsize_to_get : A single window size for which to get the data
  #$stat_to_get : Which statistic (from array in @getRequestSpace) is being used for this figure
  #$mu_lvl : The (zero-) index of the mu value to use based on array in @getParameterSpace
  #$K_lvl : The (zero-) index of the K value to use based on array in @getParameterSpace
  #$init_lvl : The (zero-) index of the init value to use based on array in @getParameterSpace
  #<-Returns:
  #$plot_data : The min/max/mean data for the requested required parameter space across all run-times for the requested sample size, window size, and statistic.
}
getRequiredPlotDataFrame <- function(samplesize_to_get, windowsize_to_get, stat_to_get, mu_lvl, K_lvl, init_lvl, filt_val = 0.02)
{
  getData <- function(model_prefix, required_array, model_suffix, x_offsets, cur_x, filter_value)
  {
    tmp_data = data.frame("Name" = character(), "x" = numeric(), "Min" = numeric(), "Max" = numeric(), "Mean" = numeric());
    for (i in 1:length(required_array))
    {
      cur_model = paste(model_prefix, ".", required_array[i], ".", model_suffix, sep = "");
      cur_data = getStatForBranch(Data_Storage[cur_model], samplesize_to_get, windowsize_to_get, stat_to_get)
      if(!is.null(cur_data))
      {
        if(stat_to_get == "filtSNPs")
        {
          cur_data = subset(cur_data, filter==filter_value);
        }
        tmp_data[nrow(tmp_data)+1,] = list(cur_model, x_offsets[i] + cur_x, min(cur_data$value), max(cur_data$value), mean(cur_data$value));
      }
      else
      {
        print("No Data");
        #tmp_data[nrow(tmp_data)+1,] = list(cur_model, x_offsets[i] + cur_x, 0, 0, 0); #This line will add a bump at the x-axis for models missing all data
        tmp_data[nrow(tmp_data)+1,] = list(cur_model, x_offsets[i] + cur_x, NA, NA, NA); #This line will add a row of NAs for models missing all data
      }
    }
    return(tmp_data);
  }
  plot_data = data.frame("Name" = character(), "x" = numeric(), "Min" = numeric(), "Max" = numeric(), "Mean" = numeric());
  cur_x = 1;
  x_off = c(0, 122, 244, 366);
  req_base = (1*mu_lvl)+(3*K_lvl)+(9*init_lvl);
  req_ary = c((req_base+81), req_base, (req_base+27), (req_base+54));
  
  #Bottleneck Models
  plot_data = rbind(plot_data, getData("8", req_ary, "-1.-1.-1", x_off, cur_x, filt_val));
  cur_x = cur_x + 1;
  
  #Recombination Models
  for(recomb_lvl in 0:2)
  {
    plot_data = rbind(plot_data, getData("9", req_ary, paste(recomb_lvl,".-1.-1", sep=""), x_off, cur_x, filt_val));
    cur_x = cur_x + 1;
  }
  
  #Progeny Skew Models
  for(pskew_lvl in 0:8)
  {
    for(recomb_lvl in 0:2)
    {
      plot_data = rbind(plot_data, getData("11", req_ary, paste(recomb_lvl,".",pskew_lvl,".-1", sep=""), x_off, cur_x, filt_val));
      cur_x = cur_x + 1;
    }
  }
  
  #DFE Models
  for(dfe_lvl in 0:2)
  {
    for(pskew_lvl in 0:8)
    {
      for(recomb_lvl in 0:2)
      {
        plot_data = rbind(plot_data, getData("15", req_ary, paste(recomb_lvl,".",pskew_lvl,".",dfe_lvl, sep=""), x_off, cur_x, filt_val));
        cur_x = cur_x + 1;
      }
    }
  }
  
  return(plot_data);
}