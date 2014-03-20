epi.info <-
function(betas,pa=NULL,pb=NULL,model=NULL,distn,ve=NA){

	#####
	# Validate input
	# Note: extra validate assumed to be done by calling function
	#####

	### Validate allele freqs
	
#	if(is.null(pa)){
#		warning("Allele frequency for first locus not given. Defaulting to pa=.5")
#		pa <- .5
#	}
#	if(is.null(pb)){
#		warning("Allele frequency for second locus not given. Defaulting to pb=.5")
#		pb <- .5
#	}
#	if(!is.numeric(pa) || !is.numeric(pb)){
#		stop("Allele frequencies pa and pb must be numeric")
#	}
#	if(any(pa <= 0) || any(pa >= 1) || any(pb <= 0) || any(pb >= 1)){
#		stop("Allele frequencies pa and pb must be between 0 and 1")
#	}
#	if((length(pa)==3 && (abs(sum(pa)-1)>.Machine$double.eps ^ 0.5)) || (length(pb)==3 && (abs(sum(pb)-1)>.Machine$double.eps ^ 0.5))){
#		stop("Vector of genotype frequencies must sum to 1")
#	}else if((length(pa)>1 && !(length(pa)==3)) || (length(pb)>1 && !(length(pb)==3))){
#		stop("pa and pb must be either single numeric values giving allele frequencies for the two loci, or length 3 numeric vectors giving genotype frequencies")
#	}

	# validate betas
	if(!(length(betas)==9)){
		stop("betas must be a vector of the 9 parameter values for the selected linear model")
	}
	if(!is.numeric(betas)){
		stop("The parameter values in 'betas' must be numeric")
	}
	if(any(is.na(betas))){
		stop("All parameter values must be specified in 'betas'")
	} 

	# validate model spec
#	if(!(distn=="logit") && !(distn=="normal") && !(distn=="probit")){
#		stop("GLM family must be 'normal', 'logit', or 'probit'.")
#	}
	
	if(is.null(model)){
		warning("Model not specified. Defaulting to 'NOIA_S'.")
		model <- "NOIA_S"
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
	
	
	#####
	# Use betas from NOIA_S, convert from input if needed
	#####
	
	if(model=="NOIA_S"){
	
		B.V <- betas
		names(B.V) <- c("mu","alpha_a","delta_a",
			"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
			"alphadelta_ab","deltadelta_ab")
			
	}else{
	
		B.V <- epi.convert(modelin=model,params=betas,modelout="NOIA_S",
					pa=pa,pb=pb,distn=distn,ve=NA,verbose=FALSE)$NOIA_S
					
	}
	
	#####
	# Conditional Expected Values
	#####

	### NOIA_S design matrix
	# denominators for recessive effects
	denomA <- pa.V[1]+pa.V[3]-(pa.V[1]-pa.V[3])^2
	denomB <- pb.V[1]+pb.V[3]-(pb.V[1]-pb.V[3])^2
	
	# univariate design matrices
	designA.M <- matrix(c(
		1,-pa.V[2]-2*pa.V[3],-2*pa.V[2]*pa.V[3]/denomA,
		1,1-pa.V[2]-2*pa.V[3],4*pa.V[1]*pa.V[3]/denomA,
		1,2-pa.V[2]-2*pa.V[3],-2*pa.V[1]*pa.V[2]/denomA),3,3,byrow=TRUE)
	
	designB.M <- matrix(c(
		1,-pb.V[2]-2*pb.V[3],-2*pb.V[2]*pb.V[3]/denomB,
		1,1-pb.V[2]-2*pb.V[3],4*pb.V[1]*pb.V[3]/denomB,
		1,2-pb.V[2]-2*pb.V[3],-2*pb.V[1]*pb.V[2]/denomB),3,3,byrow=TRUE)
	
	# combine two loci in design
	# rearrange columns for nicer presentation
	design.M <- kronecker(designB.M,designA.M)
	design.M <- cbind(design.M[,1:4],design.M[,7],design.M[,5:6],design.M[,8:9])
	
	
	# get expected values
	XB <- design.M %*% as.matrix(B.V,ncol=1)
	Eu <- matrix(XB,3,3,byrow=TRUE)
	
	if(distn=="normal"){
		Ey <- Eu
	}else if(distn=="logit"){
		Ey <- exp(Eu)/(1+exp(Eu))
	}else if(distn=="probit"){
		Ey <- pnorm(Eu)
	}

	colnames(Ey) <- c("AA","Aa","aa")
	rownames(Ey) <- c("BB","Bb","bb")

	#####
	# Marginals
	#####

	# joint genotype probabilities
	p.M <- matrix(pb.V,3,1) %*% matrix(pa.V,1,3)

	# expected marginal phenotype values
	marA <- colSums(p.M*Ey)
	marB <- rowSums(p.M*Ey)
	marAB <- sum(p.M*Ey)

	#####
	# Phenotypic variance
	#####

	if(distn=="normal"){
		if(is.na(ve)){
			vy <- 1
			ve <- 1-sum(p.M*((Ey-marAB)^2))
		}else{
			vy <- sum(p.M*((Ey-marAB)^2))+ve
		}
	}else if(distn=="logit"){
#		if(!is.na(ve)){
#			warning("ve is ignored for logit model")
#		}
		ve <- (pi^2)/3
		vy <- sum(p.M*((Eu-sum(p.M*Eu))^2))+ve
	}else if(distn=="probit"){
#		if(!is.na(ve)){
#			warning("ve is ignored for probit model")
#		}
		ve <- 1
		vy <- sum(p.M*((Eu-sum(p.M*Eu))^2))+ve
	}
	
	# multiple R^2
	RR <- (vy-ve)/vy

	#####
	# Finish
	#####

	if(distn=="normal"){
		expect <- list(Ey,marA,marB,marAB)
		names(expect) <- c("Conditional","MarginA","MarginB","Unconditional")
		out <- list(expect,vy,ve)
		names(out) <- c("Expected_Values","Var_y","error_variance")
		
	}else{ # logit or probit
		expect <- list(Ey,marA,marB,marAB)
		names(expect) <- c("Penetrances","MarginA","MarginB","Prevalence")
		out <- list(expect,vy,ve)
		names(out) <- c("Expected_Values","Var_ystar","error_variance")
	}
	return(out)
}

