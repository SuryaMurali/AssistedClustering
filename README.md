# AssistedClustering
#-------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------#
#  This program is to perform Automated K-means clustering combined with Ward's clustering, #
#  focusing on simultaneous variable selection and outlier detection.                       #
#                                                                                           #
# <Simultaneous Variable Selection and outlier detection for Automated K-means clustering>  # 
#                                                                                           #
#  Step 0 : initialization - Data Transformation, Sampling                                  #
#       0-1 : Detect outliers(global outlier) before clustering( Not Used ! )               #
#  Step 1 : Perform automated K-means clustering for each variable                          # 
#  Step 2 : Adjusted Rand Index for each pair of variables                                  #
#  Step 3 : Find highest top 5 adjusted Rand Index and calculate ratio of SSB/SST           #
#  Step 4 : Select two variables which have the highest ratio of SSB/SST                    #
#  Step 5 : Run K-means using the selected variables                                        #
#       5-1 : Detect outliers for each cluster (local outlier)                              #
#       5-2 : Detect global outliers                                                        #
#  Step 6 : K-means for unselected variables and  Compute Adjusted Rand Index               #
#           between Selected Vars and Unselected Vars                                       #
#  Step 7 : If max rand index > critical value, add selected variable                       #
#              and go to Step 5                                                             #
#           Else Stop                                                                       #
#  Step 8 : Last Kmeans using selected variables                                            #
#       8-1 : Detect outliers for each cluster (local outlier)                              #
#       8-2 : Detect global outliers                                                        #
                                                          #
#                                                                                           #  
#  Needed package for this program                                                          #
#    install.packages("sampling")                                                           #
#    install.packages("mclust")                                                             #
#                                                                                           #
#  How to Run :                                                                             #

# Very Simple!
#    result <- SVOKmeans(data=<your_dataset>)                                               #
#    Then follow the interactive prompt to specify necessary parameters                     #
