varex.to.noia <-
function(varex,pa=NULL,pb=NULL,distn="logit",mu=NA,ve=NA,signs=rep(1,8),totaleff=NA){

	#####
	# check for valid input
	#####

	# check for valid variance components
	if(!(length(varex)==8) || !is.numeric(varex)){
		stop("varex should be a length 8 vector of proportions of variance explained")
	}else if(any(varex<0) || any(varex>1) || sum(varex)>1 ){
		stop("Variance explained cannot be <0 or >1")
	}
	
	# check for valid allele freqs
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
	if((length(pa)==3 && (abs(sum(pb)-1)>.Machine$double.eps ^ 0.5)) || (length(pb)==3 && (abs(sum(pb)-1)>.Machine$double.eps ^ 0.5))){
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

	# check for valid distribution
	if(!(distn=="logit") && !(distn=="normal") && !(distn=="probit")){
		warning("Distribution must be 'normal', 'logit', or 'probit'. Defaulting to 'logit'.")
		distn <- "logit"
	}
	
	# force signs to be abs(signs)=1
	signs <- sign(signs)

	#####
	# initialize
	#####

	B.V <- rep(NA,9)

	if(distn=="logit"){
		if(is.na(totaleff)){
			RR <- sum(varex)
		}else{
			RR <- totaleff
		}
		if(!is.na(ve)){
			warning("ve is ignored for logit model")
		}
		ve <- (pi^2)/3
		vy <- ve/(1-RR)
	}else if(distn=="normal"){
		if(is.na(totaleff)){
			RR <- sum(varex)
		}else{
			RR <- totaleff
		}
		if(is.na(ve)){
			vy <- 1
			ve <- 1-RR
		}else{
			vy <- ve/(1-RR)
		}
	}else if(distn=="probit"){
		if(is.na(totaleff)){
			RR <- sum(varex)
		}else{
			RR <- totaleff
		}
		if(!is.na(ve)){
			warning("ve is ignored for probit model")
		}
		ve <- 1
		vy <- ve/(1-RR)
	}else{
		stop("GLM distribution must be either 'normal', 'logit', or 'probit'.")
	}

	#####
	# variance of design variables
	#####
	
	VaddA <- pa.V[2]*(1-pa.V[2])+4*pa.V[3]*(1-pa.V[2])-4*pa.V[3]^2
	VaddB <- pb.V[2]*(1-pb.V[2])+4*pb.V[3]*(1-pb.V[2])-4*pb.V[3]^2
	VdomA <- 4*pa.V[3]*pa.V[2]*pa.V[1]/VaddA
	VdomB <- 4*pb.V[3]*pb.V[2]*pb.V[1]/VaddB
	
	#####
	# Additive effects
	#####
	B.V[2] <- signs[1] * sqrt(varex[1]*vy/VaddA)
	B.V[4] <- signs[3] * sqrt(varex[3]*vy/VaddB)
	
	#####
	# dominance effects
	#####
	B.V[3] <- signs[2] * sqrt(varex[2]*vy/VdomA)
	B.V[5] <- signs[4] * sqrt(varex[4]*vy/VdomB)
	
	#####
	# interaction effects
	#####
	B.V[6] <- signs[5]*sqrt(varex[5]*vy/(VaddA*VaddB))
	B.V[7] <- signs[6]*sqrt(varex[6]*vy/(VdomA*VaddB))
	B.V[8] <- signs[7]*sqrt(varex[7]*vy/(VaddA*VdomB))
	B.V[9] <- signs[8]*sqrt(varex[8]*vy/(VdomA*VdomB))


	#####
	# if mu not given, compute so E(y)=0
	#####
	
	if(is.na(mu)){
	
		# NOIA design matrix, minus mu
		
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
			design.M <- kronecker(designB.M,designA.M)
			design.M <- cbind(design.M[,2:4],design.M[,7],design.M[,5:6],design.M[,8:9])
		
		# Genotype probabilities
		pab.M <- matrix(pb.V,nrow=3) %*% t(matrix(pa.V,nrow=3))
		
		# E(XB) with computed coefs (minus mu)
		EXB <- sum((design.M %*% as.matrix(B.V[2:9],nrow=8)) * as.vector(t(pab.M)))
		
		# set mu so E(y)=0
		B.V[1] <- -EXB
	
	
	}else{
		B.V[1] <- mu
	}

	
	####
	# Output model parameters
	####
	
	names(B.V) <- c("mu","alpha_a","delta_a",
	"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
	"alphadelta_ab","deltadelta_ab")
	return(B.V)
}

