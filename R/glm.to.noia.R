glm.to.noia <- 
function(betas,model,pa=NULL,pb=NULL){

	#####
	# Validate input
	#####

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
	# Get design matrix for input, output model
	#####

	### input design matrix

	if(model=="NOIA_S"){
	
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
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
	
	}else if(model=="NOIA_F"){
	
		# univariate design matrices
		designA.M <- matrix(c(
			1,-pa.V[2]-2*pa.V[3],-pa.V[2],
			1,1-pa.V[2]-2*pa.V[3],1-pa.V[2],
			1,2-pa.V[2]-2*pa.V[3],-pa.V[2]),3,3,byrow=TRUE)
	
		designB.M <- matrix(c(
			1,-pb.V[2]-2*pb.V[3],-pb.V[2],
			1,1-pb.V[2]-2*pb.V[3],1-pb.V[2],
			1,2-pb.V[2]-2*pb.V[3],-pb.V[2]),3,3,byrow=TRUE)
	
		# combine two loci in design
		# rearrange columns for nicer presentation
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
	
	}else if(model=="G2A"){
	
		# allele freqs
		pa1 <- pa.V[3] + pa.V[2]/2
		pb1 <- pb.V[3] + pb.V[2]/2
	
		# univariate design matrices
		designA.M <- matrix(c(
			1,-2*(1-pa1),-2*(1-pa1)^2,
			1,pa1-(1-pa1),2*pa1*(1-pa1),
			1,2*pa1,-2*pa1^2),3,3,byrow=TRUE)
	
		designB.M <- matrix(c(
			1,-2*(1-pb1),-2*(1-pb1)^2,
			1,pb1-(1-pb1),2*pb1*(1-pb1),
			1,2*pb1,-2*pb1^2),3,3,byrow=TRUE)
	
		# combine two loci in design
		# rearrange columns for nicer presentation
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
		
	}else if(model=="F2"){
	
		# univariate design matrices
		designA.M <- matrix(c(
			1,-1,-1/2,
			1,0,1/2,
			1,1,-1/2),3,3,byrow=TRUE)
	
		designB.M <- designA.M
	
		# combine two loci in design
		# rearrange columns for nicer presentation
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
		
	}else if(model=="Finf"){
	
		# univariate design matrices
		designA.M <- matrix(c(
			1,-1,0,
			1,0,1,
			1,1,0),3,3,byrow=TRUE)
	
		designB.M <- designA.M
	
		# combine two loci in design
		# rearrange columns for nicer presentation
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
		
	}else if(model=="unweight"){
	
		# univariate design matrices
		designA.M <- matrix(c(
			1,-1,-1/3,
			1,0,2/3,
			1,1,-1/3),3,3,byrow=TRUE)
	
		designB.M <- designA.M
	
		# combine two loci in design
		# rearrange columns for nicer presentation
		d1.M <- kronecker(designB.M,designA.M)
		d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
	
	}else if(model=="genotype"){
	
		d1.M <- matrix(c(
			1,0,0,0,0,0,0,0,0,
			1,1,0,0,0,0,0,0,0,
			1,0,1,0,0,0,0,0,0,
			1,0,0,1,0,0,0,0,0,
			1,1,0,1,0,1,0,0,0,
			1,0,1,1,0,0,1,0,0,
			1,0,0,0,1,0,0,0,0,
			1,1,0,0,1,0,0,1,0,
			1,0,1,0,1,0,0,0,1),9,9,byrow=TRUE)
	
	}else{
		stop("Invalid input model specification. Must be one of 'NOIA_S', 'NOIA_F', 'G2A', 'F2', 'Finf', 'unweight', 'genotype'.")
	}
	
	### output design matrix

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
	d2.M <- kronecker(designB.M,designA.M)
	d2.M <- cbind(d2.M[,1:4],d2.M[,7],d2.M[,5:6],d2.M[,8:9])
	
	#####
	# Convert parameters
	#####

	B.V <- as.vector(solve(d2.M) %*% d1.M %*% as.matrix(betas,9,1))
	names(B.V) <- c("mu","alpha_a","delta_a",
	"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
	"alphadelta_ab","deltadelta_ab")

	return(B.V)
}