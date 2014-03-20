pen.to.noia <-
function(pen.M,pa=NULL,pb=NULL,distn="logit"){

##########
# Computes NOIA_S regression coefficients for penetrance models
##########

	#####
	# check for valid input
	#####
	
	### Validate matrix of penetrances

	if(!is.numeric(pen.M) || !(ncol(pen.M)==3) || !(nrow(pen.M)==3) || any(pen.M<=0) || any(pen.M>=1)){
		stop("pen.M must be a 3x3 penetrance matrix with values between 0 and 1")
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

	#####
	# compute model parameters
	#####

	pen.V <- matrix(as.vector(t(pen.M)),ncol=1)
	
	# convert link funciton
	if(distn=="logit"){
		u.V <- log(pen.V/(1-pen.V))
	}else if(distn=="probit"){
		u.V <- qnorm(pen.V)
	}else{
		stop("Distribution must be either 'logit' or 'probit'.")
	}
	
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
	design.M <- cbind(design.M[,1:4],design.M[,7],design.M[,5:6],design.M[,8:9])
	
	B.V <- as.vector(solve(design.M) %*% u.V)
	names(B.V) <- c("mu","alpha_a","delta_a",
	"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
	"alphadelta_ab","deltadelta_ab")

	#####
	# Return parameters
	#####

	return(B.V)
}

