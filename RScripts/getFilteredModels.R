#Requires data to have been loaded into the workspace already. Checks each data branch in $Data_Storage against a set of thresholds/filters
#Models that clear the threshold are written to an output file.

#Define management constants, namely the filtering criteria being used.
{
  Filter_Sample_Size = setNames("100", "_partial_100")
  Filter_Window_Size = setNames("30000", "Full")
  Stat_To_Filter = "filtSNPs";
  Filter_Min = 0;
  Filter_Max = 5;
  Filter_UseMinMax = c(TRUE, TRUE);
  Filter_Inclusive = c(FALSE, TRUE);
  Variant_Filter_Value = 0.02;
  Filtered_Model_File = "filteredModels.out.txt";
}

for (i in c(1:length(Data_Storage)))
{
  if (filterByStat(Data_Storage[i], Filter_Sample_Size, Filter_Window_Size, Stat_To_Filter, variant_filter_value=Variant_Filter_Value, min=Filter_Min, max=Filter_Max, useMinMax = Filter_UseMinMax, inclusive=Filter_Inclusive))
  {
    cat(names(Data_Storage[i]), paste(Data_Storage[[i]][[1]], collapse=","), file = Filtered_Model_File, append = TRUE, sep = ",");
    cat("\n", file = Filtered_Model_File, append = TRUE);
  }
}