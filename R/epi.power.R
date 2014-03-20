epi.power <- function(truebetas,model,distn,full=rep(TRUE,9),reduced,N,alpha=.05,pa=NULL,pb=NULL,ve=NA,plotit=FALSE){

	#####
	# Check valid input
	#####
	
	# validate full model
	if(!(length(truebetas)==9)){
		stop("'truebetas' must be a vector of the 9 parameter values for the saturated two-locus GLM model")
	}
	if(!is.numeric(truebetas)){
		stop("The parameter values in 'truebetas' must be numeric")
	}
	if(any(is.na(truebetas))){
		stop("All parameter values for the saturated model must be specified in 'truebetas'")
	} 
	
	# validate specification of full, reduced model for testing
	if(!is.logical(full) || !is.logical(reduced) || !(length(full)==9) || !(length(reduced)==9) ){
		stop("The values 'full' and 'reduced' must be logical vectors of length 9 indicating whether each of the 9 parameters for the saturated model are included in the full and reduced model, respectively")
	}
	if(any(reduced & !full)){
		stop("Parameter(s) #",paste(which(reduced & !full),sep=",")," only present in the reduced model. The reduced model must be a subset of the full model.")
	}
	if(all(full==reduced)){
		stop("The vectors 'full' and 'reduced' must differ to indicate the desired hypothesis test")
	}
	
	# validate sample size(s)
	if(missing(N)){
		stop("Must specify sample size N")
	}
	if(!all(round(N)==N)){
		N <- round(N)
		warning(paste("N must be an integer. Rounding to N=",N,sep=""))
	}
	if(!is.numeric(N) || any(N<1) ){
		stop("N must be positive and numeric")
	}
	
	### Validate allele freqs
	if(is.null(pa)){
		warning("Allele frequency for first locus not given. Defaulting to pa=.5")
		pa <- .5
	}
	if(is.null(pb)){
		warning("Allele frequency for second locus not given. Defaulting to pb=.5")
		pb <- .5
	}
	if(!is.numeric(pa) || !is.numeric(pb)){
		stop("Allele frequencies pa and pb must be numeric")
	}
	if(any(pa <= 0) || any(pa >= 1) || any(pb <= 0) || any(pb >= 1)){
		stop("Allele frequencies pa and pb must be between 0 and 1")
	}
	if((length(pa)==3 && (abs(sum(pa)-1)>.Machine$double.eps ^ 0.5)) || (length(pb)==3 && (abs(sum(pb)-1)>.Machine$double.eps ^ 0.5))){
		stop("Vector of genotype frequencies must sum to 1")
	}else if((length(pa)>1 && !(length(pa)==3)) || (length(pb)>1 && !(length(pb)==3))){
		stop("pa and pb must be either single numeric values giving allele frequencies for the two loci, or length 3 numeric vectors giving genotype frequencies")
	}
	
	### convert all allele freqs to genotype freq vectors
	
	if(length(pa) == 1){
		pa.V <- c((1-pa)^2, 2*pa*(1-pa), pa^2)
	}else if(length(pa)==3){
		pa.V <- pa
	}
	if(length(pb) == 1){
		pb.V <- c((1-pb)^2, 2*pb*(1-pb), pb^2)
	}else if(length(pb)==3){
		pb.V <- pb
	}
	
	# verify other values
	
	if(length(alpha)>1 || !is.numeric(alpha) || alpha >= 1 || alpha <= 0){
		stop("alpha must be a single numeric value between 0 and 1 giving the threshold for significance")
	}
	
	if( !is.logical(plotit) || is.na(plotit) || length(plotit)>1 ){
		warning("plotit must be either TRUE or FALSE. Defaulting to TRUE.")
		plotit <- TRUE
	}

	if(!is.na(ve) && !(distn=="normal")){
		warning("ve is not used for 'logit' or 'probit' distribution model")
		ve <- NA
	}
	
	#####
	# Compute power
	#####
	
	if(distn=="normal"){
	
		# use "pwr" to compute power from Cohen f^2 effect size
#		require("pwr")
	
		# get base design matrix
		x <- rep(c(0,1,2),3)
		z <- c(rep(0,3),rep(1,3),rep(2,3))
		design.M <- epi.design(x=x,z=z,model=model,pa=pa.V,pb=pb.V)
		ytrue <- as.matrix(design.M) %*% as.matrix(truebetas,ncol=1)
	
		# get vector of cell (genotype) probabilities
		p.V <- pa.V[x+1] * pb.V[z+1]
	
		# get variance pieces of true model
		mod <- epi.info(betas=truebetas,pa=pa.V,pb=pb.V,model=model,distn=distn,ve=ve)
		vary <- mod$Var_y
		vare <- mod$error_variance
	
		# get multiple R^2 between full, true models
		if(!all(full)){
			fit1 <- lm(ytrue~.,data=cbind(ytrue,design.M[,full, drop=FALSE]),weights=p.V)
			RRft <- summary(fit1)$r.squared
		}else{
			RRft <- 1
		}
		
		# multiple R^2 for full model
		RR_f <- (RRft * (vary-vare))/vary
		
		# get multiple R^2 between true, reduced models
		fit2 <- lm(ytrue~.,data=cbind(ytrue,design.M[,reduced, drop=FALSE]),weights=p.V)
		RRrt <- summary(fit2)$r.squared
		
		# multiple R^2 for reduced model
		RR_r <- (RRrt * (vary-vare))/vary
		
		# compute effect size f^2
		ff <- (RR_f-RR_r)/(1-RR_f)
		
		# get power for effect size
		u <- sum(full)-sum(reduced) 	# numerator df
		v <- N-sum(full) 				# denominator df
		pow <- pwr.f2.test(u=u,v=v,f2=ff,sig.level=alpha)$power
	
	}else if(distn=="logit"){
	
		# use asypow package for power from asymptotic likelihood ratio
		# require() here and not Import due to asypow's license
		if(!suppressWarnings(require("asypow"))){
			stop("Computing asymptotic power of likelihood ratio test for logit models requires package 'asypow'; please install asypow and retry. Note, use of asypow is subject to the ACM license (http://www.acm.org/publications/policies/softwarecrnotice).")
		}
		
		# get base design matrix
		x <- rep(c(0,1,2),3)
		z <- c(rep(0,3),rep(1,3),rep(2,3))
		design.M <- epi.design(x=x,z=z,model=model,pa=pa.V,pb=pb.V)
		
		# get vector of cell (genotype) probabilities
		p.V <- pa.V[x+1] * pb.V[z+1]
		
		# if not full model not saturated, get parameters
		if(any(!full)){
			ytrue <- as.matrix(design.M) %*% as.matrix(truebetas,ncol=1)
			fit1 <- lm(ytrue~.,data=cbind(ytrue,design.M[,full, drop=FALSE]),weights=p.V)
			B.V <- as.vector(coef(fit1))
		}else{
			B.V <- truebetas
		}
		
		# compute power with asypow
		inf.M <- asypow::info.mvlogistic(coef=B.V,design=as.matrix(design.M[,full]),rss=p.V)
		
		# setup constraint matrix 
		#(sets parameters not in reduced model to 0; see asypow documentation)
		con.M <- matrix(NA,sum(full)-sum(reduced),3)
		con.M[,1] <- 1
		con.M[,2] <- which(full&(!reduced))
		con.M[,3] <- 0
		
		# non-centrality parameter
		nc <- asypow::asypow.noncent(theta.ha=B.V,info.mat=inf.M,constraints=con.M)
		
		# loop sample sizes and compute with asypow
		pow <- rep(NA,length(N))
		for(i in 1:length(N)){
			pow[i] <- asypow::asypow.power(nc,N[i],alpha)
		}
	
	}else if(distn=="probit"){
	
		# use asypow package for power from asymptotic likelihood ratio
		# require() here and not Import due to asypow's license
		if(!suppressWarnings(require("asypow"))){
			stop("Computing asymptotic power of likelihood ratio test for probit models requires package 'asypow'; please install asypow and retry. Note, use of asypow is subject to the ACM license (http://www.acm.org/publications/policies/softwarecrnotice).")
		}
		
		# get base design matrix
		x <- rep(c(0,1,2),3)
		z <- c(rep(0,3),rep(1,3),rep(2,3))
		design.M <- epi.design(x=x,z=z,model=model,pa=pa.V,pb=pb.V)
		
		# get vector of cell (genotype) probabilities
		p.V <- pa.V[x+1] * pb.V[z+1]
		
		# if not full model not saturated, get parameters
		if(any(!full)){
			ytrue <- as.matrix(design.M) %*% as.matrix(truebetas,ncol=1)
			fit1 <- lm(ytrue~.,data=cbind(ytrue,design.M[,full, drop=FALSE]),weights=p.V)
			B.V <- as.vector(coef(fit1))
		}else{
			B.V <- truebetas
		}
		
		# loop design cells to get expected information
		# computing manually since asypow doesn't include probit
		inf.M <- matrix(0,ncol(design.M[,full]),ncol(design.M[,full]))
		for(i in 1:nrow(design.M)){
			vars <- matrix(as.numeric(design.M[i,full]),ncol=1)
			b <- matrix(B.V,ncol=1)
			xb <- t(vars) %*% b
			# require index [1,1] to use as scaler
			EI <- (dnorm(xb)*dnorm(xb)/(pnorm(xb)*(1-pnorm(xb))))[1,1] * vars %*% t(vars)
			inf.M <- inf.M + p.V[i]*EI
		}
		
		### compute power with asypow
		
		# setup constraint matrix 
		#(sets parameters not in reduced model to 0; see asypow documentation)
		con.M <- matrix(NA,sum(full)-sum(reduced),3)
		con.M[,1] <- 1
		con.M[,2] <- which(full&(!reduced))
		con.M[,3] <- 0
		
		# non-centrality parameter
		nc <- asypow::asypow.noncent(theta.ha=B.V,info.mat=inf.M,constraints=con.M)
		
		# loop sample sizes
		pow <- rep(NA,length(N))
		for(i in 1:length(N)){
			pow[i] <- asypow::asypow.power(nc,N[i],alpha)
		}
	
	}else{
		stop("Distribution must be 'normal', 'logit', or 'probit'.")
	}

	####
	# Format output and finish
	####
	
	out <- data.frame(N=N,power=pow)
	
	# optional plot
	if(plotit){
		plot(out, type = "l", main = "Power Analysis")
	}
	
	return(out)
}
