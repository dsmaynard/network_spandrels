###################################################################################
# Maynard DS, Serván CA, Allesina S. "Network spandrels reflect ecological assembly".
# Ecology Letters, 2018, doi: 10.1111/ele.12912
# contact: dmaynard@uchicago.edu or sallesina@uchicago.edu
#####################################################################################


rm(list=ls())

## for integrating dynamics
library(deSolve)

# for calculating network metrics
library(igraph)

# various functions for manipulating data
library(tidyverse)

# for the Marchenko-Pastur distribution
library(RMTstat) 

# for plotting multiple grobs in ggplot
library(gridExtra)

##############################################################################################
# The main assembly function. The "scenario" argument calls one of the 7 immigration scenarios.
# see the file launch_assembly.R for examples
##############################################################################################


assemble_community <- function(target_size=50, scenario="notraits_immigration_search", a_ii=-1, r=1, show_prog=T, full_history=T, ...){
   
	### find the correct function for the specified scenario
	add_species<-match.fun(scenario)


	# initialize the first community comprised of a single species
	A <- matrix(c(a_ii), 1, 1)
	x <- -r/A[1,1]
	P <- NULL

	# initialize the one-species community
	current_comm<-list(x=x,A=A,P=NULL)
	
	history_results<-NULL
	comm_size<-1
	timestep<-1
		
	while(comm_size<target_size){
		
		# check if printing this iteration
		if (show_prog & ((timestep%%20==0 & timestep<100) | (timestep %% 100 == 0))) {
			print(paste("Iteration ",timestep,". Found ",comm_size," of ",target_size," species.",sep=""))
		}

		# get a new A (and P)
		new_AP<-add_species(x=current_comm$x,r=r,A=current_comm$A,P=current_comm$P,...)
	
		# check stability, feasbility, and integrate until found. return the new community
		current_comm<-integrate_and_prune(x=current_comm$x,r=r,A=new_AP$A,P=new_AP$P,...)
		

		# get the new community size
		comm_size<-nrow(current_comm$A)
		
		# do you want to return the assmebly trends
		if(full_history){
			history_results <- rbind(history_results, get_stats(status=current_comm,timestep=timestep,r=r,...))
		}
		
		# update the time
		timestep<-timestep+1		
	}
	
	# clean up the variables and add meta-data to the return list
	if(!is.null(current_comm$P)){
		colnames(current_comm$P)<-paste("S",1:comm_size,sep="")
		current_comm$ntraits<-nrow(current_comm$P)
	}
	rownames(current_comm$A)<-colnames(current_comm$A)<-paste("S",1:comm_size,sep="")
	names(current_comm$x)<-paste("S",1:comm_size,sep="")
	current_comm$scenario<-scenario
	current_comm$r_invader<-r

	# return only the final results or the assembly history as well
	if(full_history){
		return(append(current_comm,list(history=history_results)))
	}
	else{
		return(current_comm)
	}
}


## Calculate statistics for assembly trends
get_stats <- function(status,timestep,r,...){
	
	n <- length(status$x)
	x <- status$x
	r <- rep(r,n)
	A <- as.matrix(status$A)

	## stats on coefficients
	mu <- mean(A[upper.tri(A)])
	sigma2 <- mean(A[upper.tri(A)]^2) - mu^2
	d <- mean(diag(A))
	## stats on abundance
	mean_abundance <- mean(x)
	tot_abundance <- sum(x)

	return(data.frame(
		t = timestep,
		n =n,
		mean_interaction = mu,
		var_interaction = sigma2,
		mean_abundance = mean_abundance,
		tot_abundance = tot_abundance
	))
}

#######################################################################################################
# Functions for generating a new species trait vector, either randomly or by mutating a parent's vector
########################################################################################################

## generate a new species trait vector in the continuous traits case
new_continuous<-function(P,parent,old_tmp,ntraits=200,mutation.rate=0.03,...){

	# first check to see if we are initializing a new species
	if(is.null(P)){
		return(P=as.matrix(rnorm(ntraits)))
	}
	
	P<-as.matrix(P)
	k<-nrow(P)

	if(missing(parent)){

		if(missing(old_tmp)){
			tmp <- rnorm(k)

		}
		else{
			# if hill climbing, select a trait and resample
			tmp <- old_tmp
			nflip <- sample.int(length(tmp), 1)   
			tmp[nflip]<-rnorm(1)
		}
	}
	else{

		#add a new variable with small variation to get new
		tmp <- P[,parent]
		tmp <-tmp+rnorm(k,mean=0,sd=mutation.rate)

	}
	return(tmp)
}


## generate a new species trait vector in the binary traits case
new_binary<-function(P,parent,old_tmp,ntraits=200,bin.prob=0.5,mutation.rate=0.03,...){
 
	# parent is present only in radiation
	# old_tmp is the vector that is iteratively flipped inder radiation
	
	# first check to see if we are initializing a new species
	if(is.null(P)){
		nsuccess<-round(ntraits*bin.prob)
		tmp <- sample(c(rep(1,nsuccess),rep(0,ntraits-nsuccess)),size=ntraits,replace=F)
		return(P=as.matrix(tmp))
	}
	
	#otherwise 
	P<-as.matrix(P)
	k<-nrow(P)

	if(missing(parent)){
		if(missing(old_tmp)){
			nsuccess<-round(k*bin.prob)
			tmp <- sample(c(rep(1,nsuccess),rep(0,k-nsuccess)),size=k,replace=F)
		}
		else{
			## pick a random one and flip it to climb the hill
			tmp<-old_tmp
			ones<-seq(1:k)[tmp==1]
			zeros<-seq(1:k)[tmp==0]
			flip.ones<-sample(ones,1)
			flip.zeros<-sample(zeros,1)
			tmp[c(flip.ones,flip.zeros)]<-1-tmp[c(flip.ones,flip.zeros)]
		}
	}
	else{
		## find mutation.rate*k/2 that are =1 and the same number that are =0, and flip them
		tmp <- P[,parent]
		nflip<-max(1,min(sum(tmp),round(k*mutation.rate/2)))
		ones<-seq(1:k)[tmp==1]
		zeros<-seq(1:k)[tmp==0]
		flip.ones<-sample(ones,nflip)
		flip.zeros<-sample(zeros,nflip)
		tmp[c(flip.ones,flip.zeros)]<-1-tmp[c(flip.ones,flip.zeros)]
	}
	return(tmp)
}




#################################################################
# FUNCTIONS FOR FINDING NEW INVADERS BASED ON THE CHOSEN SCENARIO
#################################################################


### immigration, binary traits
traits_binary_immigration <- function(x, r, A, P,THRESH=1e-6,...){
	
	# initialize the trait vector if only one species so far	
	if(is.null(P)){
		P<-as.matrix(new_binary(P=NULL,...))
	}
	
	n <- length(x)
	k <- nrow(P)

	while(TRUE){
		
		# get a new trait vector
		tmp <- new_binary(P=P,...)
		
		# calculate the interaction coefficents from this new vector
		avec<-get_A_from_P_binary(P=P,newvec=tmp)

		# calculate invader's growth rate
		growth_invader <- r + sum(avec * x)
		
		# check invasion criterion
		if(growth_invader>0){
			break
		}
		else{
			while(growth_invader<0){
				# flip two indices and see if it makes things better
				tmp1<-new_binary(P=P,old_tmp=tmp,...)
				avec1<-get_A_from_P_binary(P=P,newvec=tmp1)

				# recheck invasion criterion
				growth_invader1 <- r + sum(avec1 * x)
				if(growth_invader1 > growth_invader){
					tmp<-tmp1
					growth_invader<-growth_invader1
					avec<-avec1
				}
			}
			break
		}
	}

	# add the new vector to the A matrix
	Anew <- rbind(A, avec)
	Anew <- cbind(Anew, c(avec, A[1,1]))

	# add the new vector to the P matrix	
	Pnew <- as.matrix(cbind(P, tmp))
	# Anew <- get_A_cov(Pnew)
	
	return(list(A = Anew, P = Pnew))
}


# immigration, continuous traits
traits_continuous_immigration <- function(x, r, A, P,THRESH=1e-6,...){

	# initialize the trait vector if only one species so far	
	if(is.null(P)){
		P<-as.matrix(new_continuous(P=NULL,...))
	}
	
	n <- length(x)
	k <- nrow(P)

	while(TRUE){
		
		# get a new trait vector
		tmp <- new_continuous(P=P,...)
		
		# calculate the interaction coefficents from this new vector
		avec<-get_A_from_P_continuous(P=P,newvec=tmp)

		# calculate invader's growth rate
		growth_invader <- r + sum(avec * x)
		
		# check invasion criterion
		if(growth_invader>0){
			break
		}
		else{
			while(growth_invader<0){
				# flip two indices and see if it makes things better
				tmp1<-new_continuous(P=P,old_tmp=tmp,...)
				avec1<-get_A_from_P_continuous(P=P,newvec=tmp1)

				# recheck invasion criterion
				growth_invader1 <- r + sum(avec1 * x)
				if(growth_invader1 > growth_invader){
					tmp<-tmp1
					growth_invader<-growth_invader1
					avec<-avec1
				}
			}
			break
		}
	}

	# add the new vector to the A matrix
	Anew <- rbind(A, avec)
	Anew <- cbind(Anew, c(avec, A[1,1]))

	# add the new vector to the P matrix	
	Pnew <- as.matrix(cbind(P, tmp))
	# Anew <- get_A_cov(Pnew)
	
	return(list(A = Anew, P = Pnew))
}

### radiation, binary traits
traits_binary_radiation <- function(x, r, A, P, THRESH=1e-6, ...){

	# initialize the trait vector if only one species so far
	if(is.null(P)){
		P<-as.matrix(new_binary(P=NULL,...))
	}
	
	n <- length(x)
	k <- nrow(P)
	
	while(TRUE){
		
		# select a parent to mutate
		parent <- sample.int(n, 1) 
		
		# perturb the parent's trait vector
		tmp <- new_binary(P=P,parent=parent,...)
		
		# get the new interaction coefficients
		avec<-get_A_from_P_binary(P=P,newvec=tmp)

		# check the invasion criterion
		growth_invader <- r + sum(avec * x)
		if(growth_invader>THRESH){
			break
		}	
	}
	
	# add the new vector to the A matrix
	Anew <- rbind(A, avec)
	Anew <- cbind(Anew, c(avec, A[1,1]))

	# add the new vector to the P matrix	
	Pnew <- as.matrix(cbind(P, tmp))   
	# Anew <- get_A_cov(Pnew)
	# 
	return(list(A = Anew, P = Pnew))
}


## radiation, continuous traits
traits_continuous_radiation <- function(x, r, A, P, THRESH=1e-6, ...){

	# initialize the trait matrix if only one species so far
	if(is.null(P)){
		P<-as.matrix(new_continuous(P=NULL,...))
	}
	
	n <- length(x)
	k <- nrow(P)
	
	while(TRUE){
		
		# select a parent to mutate
		parent <- sample.int(n, 1) 
		
		# perturb the parent's trait vector
		tmp <- new_continuous(P=P,parent=parent,...)
		
		# get the new interaction coefficients
		avec<-get_A_from_P_continuous(P=P,newvec=tmp)

		# check the invasion criterion
		growth_invader <- r + sum(avec * x)
		if(growth_invader>THRESH){
			break
		}	
	}
	
	# add the new vector to the A matrix
	Anew <- rbind(A, avec)
	Anew <- cbind(Anew, c(avec, A[1,1]))

	# add the new vector to the P matrix	
	Pnew <- as.matrix(cbind(P, tmp))   
	# Anew <- get_A_cov(Pnew)
	# 
	return(list(A = Anew, P = Pnew))
}


## immigration, no-traits case
notraits_immigration_search<- function(x, r, A,THRESH=1e-6, ...){

	n <- length(x)

	# get a random set of coefficents
	avec <- -runif(n)
	
	# get the invader's growth rate
	growth_invader <- r + sum(avec * x)
	
	while(growth_invader < THRESH){
		
		# resample a single coefficient
		i <- sample.int(n, 1)
		tmp1 <- avec
		tmp1[i] <- -runif(1)

		#recheck invasion criterion
		growth_rate1 <- r + sum(tmp1 * x)
		if (growth_rate1 > growth_invader){
			growth_invader <- growth_rate1
			avec<-tmp1
		}
	}

	# add the new vector to the A matrix
	Anew <- rbind(A, avec)
	Anew <- cbind(Anew, c(avec, A[1,1]))
	
	return(list(A=Anew,P=NULL))
}

## radiation, no-traits case
notraits_radiation_search <- function(x, r, A, mutation.rate=0.03, THRESH=1e-6, ...){

	n <- length(x)

	## parent
	while(TRUE){
		
		## sample a parent
		parent <- sample.int(n, 1) 
		tmp <- A[parent,]
		
		## mutate its existing interaction coefficients my scaline them with a random value
		tmp <- runif(n, 1 - mutation.rate, 1 + mutation.rate) * tmp
		tmp[tmp > 0 ] <- -0.0001
		tmp[tmp < -1] <- -1 
		
		# recheck invasion criterion
		growth_invader <- r + sum(tmp * x)
		if (growth_invader > THRESH) {
			break
		}
	}
	
	# add the new vector to the A matrix
	Anew <- rbind(A, tmp)
	Anew <- cbind(Anew, c(tmp, A[1,1]))
	
	return(list(A = Anew,P=NULL))
}



## immigration, no traits, via uniform sampling of the space of feasible solutions
notraits_immigration_unif <- function(x, r, A,...){

	# see Supplemental information for details
	
	n <- length(x)
	
	r<-rep(r,n)
	
	D <- matrix(0, n , n)
	diag(D) <- r/x

	while(TRUE){
		w <- c(0, sort(runif(n)), 1)
		y <- (w - lag(w))[2:(n+2)]
		y <- y
		v1 <- y[1:n]
		avec <- -(D %*% v1)
		if(all(avec > - 1)) break
	}
	
	# add the new vector to the A matrix
	Anew <- rbind(A, as.numeric(avec))
	Anew <- cbind(Anew, c(avec, A[1,1]))
	
	return(list(A=Anew,P=NULL))
}



#############################################################################################
# calculating the interaction coefficients (a_ij) for a new invader based on its trait values
#############################################################################################

get_A_from_P_continuous<-function(P,newvec){
	
	P<-as.matrix(P)	
	newvec<-newvec/sqrt(sum(newvec^2))
	P<-t(t(P)/sqrt(colSums(P^2))) 
	ret.val<-as.numeric(-t(P)%*%newvec)
	ret.val<-(ret.val-1)/2		
	
	return(ret.val)
}

		
get_A_from_P_binary<-function(P,newvec){

	P<-as.matrix(P)
	newvec<-newvec/sqrt(sum(newvec^2))
	P<-t(t(P)/sqrt(colSums(P^2))) 
	ret.val<-as.numeric(-t(P)%*%newvec)
	
	return(ret.val)
}

#####################################
# STABILITY, INTEGRATING AND PRUNING
#####################################

## get the rate of change at a particular point in time
LV <- function(time, state, params){

	with(as.list(params), {
		state[state < THRESH] <- 0
		dxdt <- state * (r + A %*% state)
		return(list(dxdt))
	})
}

## integrate to get long-term dynamics
integrate_LV <- function(x, r, A, THRESH=1e-6,length_of_integration=1e10,...){

	pars <- list(r = r, A = A,THRESH=THRESH)
	times <- seq(0, length_of_integration, length.out = 2)
	out <- as.data.frame(ode(x, times, LV, pars))  
	final_dens <- out[2, -1]
	final_dens[final_dens < THRESH] <- 0
	return(final_dens)
}

## prune the species that are now extinct
prune_sys <- function(x, A, P,THRESH=1e-6,...){
	extant <- x > THRESH
	
	if(!is.null(P)){
		P <- P[, extant]  
	}
	A <- A[extant, extant]
	x <- x[extant]
	return(list(x=x,A=A,P=P))
}

## check for stability
is_A_stable <- function(A){
	if (max(eigen(A, only.values = TRUE, symmetric = TRUE)$values) < 0) return(TRUE)
	return(FALSE)
}


# integrate the dynamics and prune, continuing until a stable and feasible community is obtained.
integrate_and_prune<-function(x, r, A, P,INITIAL=1e-5,THRESH=1e-6,...){
	
	# add the species in low abundance
	x<-c(x,INITIAL)
	
	## Dynamics
	while(TRUE){

		# get the growthrate vector
		rvec<-rep(r,length(x))
		
		# check whether invader is embedded in the community
		if (is_A_stable(A)){
			
			# If matrix is stable, check feasibility
			x_star <- -rvec %*% solve(A)
			
			if (sum(x_star > THRESH) == length(x_star)){
				##feasible and stable: set to equilibrium
				return(list(x=x_star,A=A,P=P))
			}

		}
		## 3) if not, 
		## 3a) integrate numerically
		x <- integrate_LV(x,rvec,A,...)

		## 3b) prune and reassign
		pruned_comm<- prune_sys(x,A,P)
		x<-pruned_comm$x
		A<-pruned_comm$A
		P<-pruned_comm$P
	}
	
}

###############################################################################
# Functions for calculating (and plotting) the observed and expected spectrum
#############################################################################

# observed spectrum of the interaction matrix
observed_spectrum_A <- function(A){
	return(data.frame(eValues=eigen(A, symmetric = TRUE)$values))
}

# observed spectrum of the scaled covariance matrix
observed_spectrum_P <- function(P){
	
	P<-scale(P)
	
	if(ncol(P)>nrow(P)){
		P<-scale(t(P))
	}
			
	svd(P)$d
	eV<-data_frame(eValues=eigen(P%*%t(P)/nrow(P))$values) %>% filter(eValues>1e-6)
	
	return(-eV)
}

# Expected distribution of the eigenvalues of A, given by the semicircle law
expected_A <- function(A){
	
	mu <- mean(A[upper.tri(A)])
	sigma2 <- mean(A[upper.tri(A)]^2) - mu^2
	d <- mean(diag(A))
	n <- nrow(A)
	x1 <- seq(-2, 2, by = 0.01)
	y1 <- sqrt(4 - x1^2) / (pi)
	x1 <- x1 * sqrt(n * sigma2) + d
	x1 <- -(x1 - mu)
	
	return(data_frame(x = -x1, density = y1))
}

# Expected distribution of the eigenvalues of the covariance mmtrix, given by Marchenko Pastur
expected_P <- function(P,min=0){
	
	if(ncol(P)>nrow(P)){
		P<-scale(t(P))
	}
	else{
		P<-scale(P)
	}
	
	evp<-observed_spectrum_P(P)
	
   	x=seq(0,max(abs(evp))+1,.001)
   	y=dmp(x,ndf=nrow(P),pdim=ncol(P))
	
   	outdat<-data_frame(x=-x,density=y) %>% filter(density>0)
   	
	return(outdat)
}


# plot the observed and expected spectra
compare_spectra<-function(outcom){
	
	if(is.null(outcom$P)){
		# if using no-traits scenarios
		observed<-observed_spectrum_A(outcom$A)
		expected<-expected_A(outcom$A)
	}
	else{
		# doing trait niche overlap model
		observed<-observed_spectrum_P(outcom$P)
		expected<-expected_P(outcom$P)
	}

	
	p1<-ggplot(observed, aes(x = eValues))+ 
		geom_histogram(aes(y = ..density..), binwidth = 0.2, colour = alpha("black",.2),alpha=0.7 ) +
		theme_bw() + 
		ylab("Density")+
		xlab(expression(paste(lambda, ", Interaction Matrix",sep=""))) + 
		geom_line(data = expected, aes(x = x, y = density), size = .5) +
		theme(legend.position="none")

	show(p1)
}



###############################################
# modularity, nestendess, motifs, plotting, etc
###############################################


# return a binary interaction matrix
binarize_A<-function(outcom,cutval=0.9){
	
	A<-outcom$A
	cutval<-quantile(A[upper.tri(A)], probs=cutval)
	Abin<-matrix(1,nrow(A),ncol(A))
	Abin[A<=cutval]<-0
	diag(Abin)<-0
	
	return(Abin)
}

# nestedness as maximum singular value
nestedness_func<-function(R){
	return(as.numeric(max(svd(R)$d)))
}

# total number of motifs of size=size
motifs_func<-function(R,size=3){
	return(count_motifs(graph_from_adjacency_matrix(R,mode="undirected"), size = size))
}

# modularity as short random walks
modularity_func<-function(R,steps=5){
	R<-graph_from_adjacency_matrix(R,mode="undirected")
	return(modularity(R,cluster_walktrap(R,steps=steps)$membership))
}


	
	
	

