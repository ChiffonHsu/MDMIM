> library(xcms)
> 
  > # create PeakGroupsParam with defaults
  > pg_params <- PeakGroupsParam()
  > 
    > # run retention time alignment using PeakGroups method
    > xdata3_aligned <- adjustRtime(xdata3, param = pg_params)
    Performing retention time correction using 33 peak groups.
    [=============================================================================================] 1023/1023 (100%) in  2s
    Applying retention time adjustment to the identified chromatographic peaks ... OK
    > 
      > # check RT differences after alignment
      > summary(rtime(xdata3_aligned, adjusted = TRUE) - rtime(xdata3_aligned, adjusted = FALSE))
    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    -76.07127  -1.69245   0.01310   0.02719   1.37364  44.73131 
    > 
      > # plot adjusted RT deviation
      > massprocesser::plot_adjusted_rt(xdata3_aligned, group_for_figure = "mix")