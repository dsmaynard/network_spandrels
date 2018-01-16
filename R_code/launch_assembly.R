
rm(list=ls())


source("assembly_functions.R")

# optional arguments: ntraits, bin.prob
outcom<-assemble_community(50,scenario="traits_binary_immigration")

# optional arguments: ntraits
outcom<-assemble_community(50,scenario="traits_continuous_immigration")

# optional arguments: ntraits, bin.prob, mutation.rate
outcom<-assemble_community(50,scenario="traits_binary_radiation")

# optional arguments: ntraits, mutation.rate
outcom<-assemble_community(50,scenario="traits_continuous_radiation")

# optional arguments: NONE
outcom<-assemble_community(100,scenario="notraits_immigration_search")

# optional arguments: mutation.rate
outcom<-assemble_community(50,scenario="notraits_radiation_search") 

# optional arguments: NONE
outcom<-assemble_community(50,scenario="notraits_immigration_unif")



#####################
# PLOTTING RESULTS
#####################

# Plotting observed vs. expected seigenvalue distributions. 
# The null density line is semicricle law if no-traits; Marchenko-Pastur if traits
compare_spectra(outcom)

## plotting trends
p1<-ggplot(outcom$history,aes(x=t,y=n))+geom_line() 
p2<-ggplot(outcom$history,aes(x=t,y=tot_abundance))+geom_line()
p3<-ggplot(outcom$history,aes(x=t,y=mean_abundance))+geom_line()
p4<-ggplot(outcom$history,aes(x=t,y=mean_interaction))+geom_line()
p5<-ggplot(outcom$history,aes(x=t,y=var_interaction))+geom_line()   
  

#plot
grid.arrange(p1,p2,p3,p4,p5,ncol=1)


###########################################################
# calculating binary interaction matrix and network metrics
##############################################################

# convert the interaction matrix to binary, setting the top 10% equal to one
Abin<-binarize_A(outcom,cutval=0.9)

# calculate network metrics
motifs_func(Abin)
nestedness_func(Abin)
modularity_func(Abin)

