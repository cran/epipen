epi.data <-
function(betas,model,distn,x=NULL,z=NULL,n=NULL,ncase=NULL,ncontrol=NULL,pa=NULL,pb=NULL,ve=NA,verbose=FALSE){

	#####
	# Validate input 
	# if needed, x,z generated from rbinom(n,2,p)
	#####

	B.V <- betas

	if(!(length(B.V)==9)){
		stop("'betas' must be a vector of the 9 parameter values for the selected GLM model")
	}
	if(!is.numeric(B.V)){
		stop("The parameter values in 'betas' must be numeric")
	}
	if(any(is.na(B.V))){
		stop("All parameter values must be specified in 'betas'")
	}
	
	if(!is.null(x)){
		if(!is.numeric(x)){
			stop("Input data must be numeric")
		}
		if(any(is.na(match(x,c(0,1,2,NA))))){
			stop("Invalid values in x. All values must be 0,1,2, or NA.")
		}
		if(!is.null(pa)){
			warning("Data given for first locus. Ignoring allele frequency pa")
			pa <- NULL
		}
		if(!is.null(n) || !is.null(ncase) || !is.null(ncontrol)){
			warning("Sample size is ignored when data is provided for either locus")
		}
		n <- length(x)
		if(!is.null(z)){
			if(!is.numeric(z)){
				stop("Input data must be numeric")
			}
			if(any(is.na(match(z,c(0,1,2,NA))))){
				stop("Invalid values in z. All values must be 0,1,2, or NA.")
			}
			if(length(z)!=n){
				stop("Data for SNPs x and z must be of the same dimension")
			}
			if(!is.null(pb)){
				warning("Data given for second locus. Ignoring allele frequency pb")
				pb <- NULL
			}
		}else{
			if(is.null(pb)){
				stop("Need to provide either data or allele frequency for the second locus")
			}
			if(!is.numeric(pb)){
				stop("pb must be numeric")
			}
			if(length(pb)==1){
				if(pb <=0 || pb >= 1){
					stop("allele frequency pb must be between 0 and 1")
				}
				z <- rbinom(n,2,pb)
			}else if(length(pb)==3){
				if(abs(sum(pb)-1)>.Machine$double.eps ^ 0.5){
					stop("genotype frequencies pb must sum to 1")
				}
				zprob <- runif(n)
				z <- rep(NA,n)
				z[zprob < pb[1]] <- 0
				z[zprob >= pb[1] & zprob < pb[1]+pb[2]] <- 1
				z[zprob >= pb[1]+pb[2]] <-2
			}else{
				stop("pb must be either the allele frequency or the three genotype frequencies for the second locus")
			}
		}
	}else{
		if(is.null(pa)){
			stop("Need to provide either data or allele frequency for the first locus")
		}
		if(!is.numeric(pa)){
			stop("pa must be numeric")
		}
		if(length(pa)==1){
			if(pa <=0 || pa >= 1){
				stop("allele frequency pa must be between 0 and 1")
			}
		}else if(length(pa)==3){
			if(abs(sum(pa)-1)>.Machine$double.eps ^ 0.5){
				stop("genotype frequencies pa must sum to 1")
			}
		}else{
			stop("pa must be either the allele frequency or the three genotype frequencies for the first locus")
		}
		if(!is.null(z)){
			if(!is.numeric(z)){
				stop("Input data must be numeric")
			}
			if(any(is.na(match(z,c(0,1,2,NA))))){
				stop("Invalid values in z. All values must be 0,1,2, or NA.")
			}
			if(!is.null(pb)){
				warning("Data given for second locus. Ignoring allele frequency pb")
				pb <- NULL
			}
			if(!is.null(n) || !is.null(ncase) || !is.null(ncontrol)){
				warning("Sample size is ignored when data is provided for either locus")
			}
			n <- length(z)
			if(length(pa)==1){
				x <- rbinom(n,2,pa)
			}else{
				xprob <- runif(n)
				x <- rep(NA,n)
				x[xprob < pa[1]] <- 0
				x[xprob >= pa[1] & xprob < pa[1]+pa[2]] <- 1
				x[xprob >= pa[1]+pa[2]] <-2
			}
		}else{
			if(is.null(pb)){
				stop("Need to provide either data or allele frequency for the second locus")
			}
			if(length(pb)>1 || !is.numeric(pb)){
				stop("pb must be a single numeric value giving allele frequency between 0 and 1")
			}

			if(is.null(n) && is.null(ncase) && is.null(ncontrol) ){
				stop("Must provide sample size to generate SNP data")
			}else{
				if(!is.null(n)){
					if(!is.null(ncase) || !is.null(ncontrol)){
						stop("Must specify either \"n\" or number of cases and controls, not both")
					}
				}else{
					if(!(dist=="logistic")){
						stop("Cases and controls can only be specified for the logistic model")
					}
					if(is.null(ncase) || is.null(ncontrol)){
						stop("Number of both cases and controls must be specified")
					}
					n <- ncase + ncontrol
				}
			}
			if(length(pa)==1){
				x <- rbinom(n,2,pa)
			}else{
				xprob <- runif(n)
				x <- rep(NA,n)
				x[xprob < pa[1]] <- 0
				x[xprob >= pa[1] & xprob < pa[1]+pa[2]] <- 1
				x[xprob >= pa[1]+pa[2]] <-2
			}
			if(length(pb)==1){
				z <- rbinom(n,2,pb)
			}else{
				zprob <- runif(n)
				z <- rep(NA,n)
				z[zprob < pb[1]] <- 0
				z[zprob >= pb[1] & zprob < pb[1]+pb[2]] <- 1
				z[zprob >= pb[1]+pb[2]] <-2
			}
		}
	}

	if(!(distn=="logit") && !(distn=="normal") && !(distn=="probit")){
		stop("GLM distribution must be 'normal', 'logit', or 'probit'.")
	}

	####
	# setup
	####

	# if using provided data, use empirical allele freqs
	
	if(is.null(pa)){
		genos <- c(0,1,2)
		pa <- table(factor(x,genos))/sum(!is.na(x))
	}
	if(is.null(pb)){
		genos <- c(0,1,2)
		pb <- table(factor(z,genos))/sum(!is.na(z))
	}
	
	# get model info
	mod <- epi.convert(modelin=model,params=betas,modelout=c(model,"varex"),
				pa=pa,pb=pb,distn=distn,ve=ve,verbose=TRUE)

	####
	# get design matrix
	####

	dat.M <- epi.design(x=x,z=z,model=model,pa=pa,pb=pb)

	####
	# compute linear model
	####

	ytrue <- as.matrix(dat.M) %*% matrix(betas,ncol=1)

	####
	# add error, dichotomize if needed
	####

	if(distn=="normal"){
		ve <- mod$Additional_info$error_variance
		y <- ytrue+rnorm(n,sd=sqrt(ve))
	}else if(distn=="logit"){
		ve <- mod$Additional_info$error_variance # expect(pi^2)/3
		ystar <- ytrue+rlogis(n)
		yprob <- exp(ystar)/(1+exp(ystar))
		y <- rbinom(n,1,yprob)
	}else if(distn=="probit"){
		ve <- mod$Additional_info$error_variance # expect 1
		ystar <- ytrue+rnorm(n,sd=sqrt(ve))
		yprob <- pnorm(ytrue)
		y <- as.numeric(ystar>0)
	}


	#####
	# loop for case/control if needed
	#####

	if((distn=="logit" || distn=="probit") && !is.null(ncase) && !is.null(ncontrol)){

		# loop to get enough cases and controls
		while(sum(y)<ncase || sum(!y)<ncontrol){
			
			# generate more genotype data
			if(length(pa)==1){
				xnew <- rbinom(n,2,pa)
			}else{
				xnewprob <- runif(n)
				xnew <- rep(NA,n)
				xnew[xnewprob < pa[1]] <- 0
				xnew[xnewprob >= pa[1] & xnewprob < pa[1]+pa[2]] <- 1
				xnew[xnewprob >= pa[1]+pa[2]] <-2
			}
			if(length(pb)==1){
				znew <- rbinom(n,2,pb)
			}else{
				znewprob <- runif(n)
				znew <- rep(NA,n)
				znew[znewprob < pb[1]] <- 0
				znew[znewprob >= pb[1] & znewprob < pb[1]+pb[2]] <- 1
				znew[znewprob >= pb[1]+pb[2]] <-2
			}
			newdat.M <- epi.design(x=xnew,z=znew,model=model,pa=pa,pb=pb)
			
			# generate more phenotype data
			if(distn=="logit"){
				newytrue <- cbind(rep(1,n),as.matrix(newdat.M)) %*% matrix(betas,ncol=1)
				newystar <- ytrue+rlogis(n)
				newyprob <- exp(newystar)/(1+exp(newystar))
				newy <- rbinom(n,1,newyprob)
			}else if(distn=="probit"){
				newytrue <- cbind(rep(1,n),as.matrix(newdat.M)) %*% matrix(betas,ncol=1)
				newystar <- ytrue+rnorm(n,sd=sqrt(ve))
				newyprob <- pnorm(newytrue)
				newy <- as.numeric(newystar>0)
			}

			# append to existing
			y <- c(y,newy)
			yprob <- c(yprob,newyprob)
			ystar <- c(ystar,newystar)
			ytrue <- c(ytrue,newytrue)
			x <- c(x,xnew)
			z <- c(z,znew)
		}

		# reduce to target total sample size if needed
		if(length(y)>(ncase+ncontrol)){
			case <- which(y==1)
			control <- which(y==0)
			pickcase <- sample(case,ncase,replace=F)
			pickcontrol <- sample(control,ncontrol,replace=F)
			saved <- sort(c(pickcase,pickcontrol))
			
			y <- y[saved]
			yprob <- yprob[saved]
			ystar <- ystar[saved]
			ytrue <- ytrue[saved]
			x <- x[saved]
			z <- z[saved]
		}
	}

	####
	# Format, output
	####

	if(verbose){
		if(distn=="normal"){
			fulldat <- data.frame(y=y,XB=ytrue,x=x,z=z)
		}else{ #distn=logit or probit
			fulldat <- data.frame(y=y,ystar=ystar,XB=ytrue,yprob=yprob,x=x,z=z)
		}
		out <- list(Full_Data=fulldat,Model_Info=mod)
		return(out)
	}else{
		outdat <- data.frame(y=y,x=x,z=z)
		return(outdat)
	}
}

