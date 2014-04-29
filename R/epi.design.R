epi.design <-
function(x,z,model=NULL,pa=NULL,pb=NULL){

	### check valid input

	if(!is.numeric(x) || !is.numeric(z)){
		stop("Input data must be numeric")
	}
	if(length(x)!=length(z)){
		stop("Data for SNPs x and z must be of the same dimension")
	}
	if(any(is.na(match(x,c(0,1,2,NA))))){
		stop("Invalid values in x. All values must be 0,1,2, or NA.")
	}
	if(any(is.na(match(z,c(0,1,2,NA))))){
		stop("Invalid values in z. All values must be 0,1,2, or NA.")
	}
	if(is.null(model)){
		warning("Model not specified. Defaulting to NOIA Statistical.")
		model <- "NOIA_S"
	}
	
	### use allele freqs if given, else use observed freqs
	if(is.null(pa) || missing(pa)){
		pa <- mean(x,na.rm=TRUE)/2
		pA <- 1-pa
		pAA <- mean(x==0,na.rm=TRUE)
		pAa <- mean(x==1,na.rm=TRUE)
		paa <- mean(x==2,na.rm=TRUE)
	}else if(length(pa)==1){
		if(!is.numeric(pa) || pa <= 0 || pa >= 1){
			stop("Allele frequency pa must be a numeric value between zero and 1")
		}
		pAA <- (1-pa)^2
		pAa <- 2*pa*(1-pa)
		paa <- pa^2
		pA <- 1-pa
	}else if(length(pa)==3){
		if(!is.numeric(pa) || any(pa <= 0) || any(pa >= 1)){
			stop("Genotype frequencies pa must be numeric values between zero and 1")
		}
		if(abs(sum(pa)-1)>.Machine$double.eps ^ 0.5){
			stop("Vector of genotype frequencies pa must sum to 1")
		}
		pAA <- pa[1]
		pAa <- pa[2]
		paa <- pa[3]
		pa <- paa + (pAa/2)
		pA <- 1-pa
	}else{
		stop("pa must be either a single numeric value giving allele frequencies for the first loci, or a length 3 numeric vector giving genotype frequencies")
	}
	if(is.null(pb) || missing(pb)){
		pb <- mean(z,na.rm=TRUE)/2
		pB <- 1-pb
		pBB <- mean(z==0,na.rm=TRUE)
		pBb <- mean(z==1,na.rm=TRUE)
		pbb <- mean(z==2,na.rm=TRUE)
	}else if(length(pb)==1){
		if(!is.numeric(pb) || pb <= 0 || pb >= 1){
			stop("Allele frequency pb must be a numeric value between zero and 1")
		}
		pBB <- (1-pb)^2
		pBb <- 2*pb*(1-pb)
		pbb <- pb^2
		pB <- 1-pb
	}else if(length(pb)==3){
		if(!is.numeric(pb) || any(pb <= 0) || any(pb >= 1)){
			stop("Genotype frequencies pb must be numeric values between zero and 1")
		}
		if(abs(sum(pb)-1)>.Machine$double.eps ^ 0.5){
			stop("Vector of genotype frequencies pb must sum to 1")
		}
		pBB <- pb[1]
		pBb <- pb[2]
		pbb <- pb[3]
		pb <- pbb + (pBb/2)
		pB <- 1-pb
	}else{
		stop("pb must be either a single numeric value giving allele frequencies for the first loci, or a length 3 numeric vector giving genotype frequencies")
	}

	### build design matrix	
	if(model=="NOIA_S"){
	
		# denominators for recessive effects
		denomA <- pAA+paa-(pAA-paa)^2
		denomB <- pBB+pbb-(pBB-pbb)^2
		
		# construct design matrix
		dmain.M <- data.frame(
			x_a=x-pAa-2*paa,
			x_d=c(-2*pAa*paa/denomA,4*pAA*paa/denomA,-2*pAA*pAa/denomA)[x+1],
			z_a=z-pBb-2*pbb,
			z_d=c(-2*pBb*pbb/denomB,4*pBB*pbb/denomB,-2*pBB*pBb/denomB)[z+1]
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("mu","alpha_a","delta_a",
			"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
			"alphadelta_ab","deltadelta_ab")

	}else if(model=="NOIA_F"){
	
		# construct design matrix
		dmain.M <- data.frame(
			x_a=x-pAa-2*paa,
			x_d=c(-pAa,1-pAa,-pAa)[x+1],
			z_a=z-pBb-2*pbb,
			z_d=c(-pBb,1-pBb,-pBb)[z+1]
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("R","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
	
	}else if(model=="G2A"){
		
		# construct design matrix
		dmain.M <- data.frame(
			x_a=c(-2*(1-pA),pA-(1-pA),2*pA)[x+1],
			x_d=c(-2*(1-pA)^2,2*pA*(1-pA),-2*pA^2)[x+1],
			z_a=c(-2*(1-pB),pB-(1-pB),2*pB)[z+1],
			z_d=c(-2*(1-pB)^2,2*pB*(1-pB),-2*pB^2)[z+1]
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
	
	}else if(model=="F2"){
	
		# construct design matrix
		dmain.M <- data.frame(
			x_a=x-1,
			x_d=c(-1/2,1/2,-1/2)[x+1],
			z_a=z-1,
			z_d=c(-1/2,1/2,-1/2)[z+1]
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
	
	}else if(model=="Finf"){
	
		# construct design matrix
		dmain.M <- data.frame(
			x_a=x-1,
			x_d=as.numeric(x==1),
			z_a=z-1,
			z_d=as.numeric(z==1)
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
	
	}else if(model=="unweight"){
	
		# construct design matrix
		dmain.M <- data.frame(
			x_a=x-1,
			x_d=c(-1/3,2/3,-1/3)[x+1],
			z_a=z-1,
			z_d=c(-1/3,2/3,-1/3)[z+1]
		)
		dint.M <- data.frame(
			xz_aa=dmain.M[,"x_a"] * dmain.M[,"z_a"],
			xz_da=dmain.M[,"x_d"] * dmain.M[,"z_a"],
			xz_ad=dmain.M[,"x_a"] * dmain.M[,"z_d"],
			xz_dd=dmain.M[,"x_d"] * dmain.M[,"z_d"]
		)
		design.M <- cbind(rep(1,length(x)),dmain.M,dint.M)
		colnames(design.M) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
	
	}else if(model=="genotype"){
	
		design.M <- data.frame(
			mu=rep(1,length(x)),
			gam_1=as.numeric(x==1),
			gam_2=as.numeric(x==2),
			beta_1=as.numeric(z==1),
			beta_2=as.numeric(z==2),
			i_11=as.numeric(x==1 & z==1),
			i_12=as.numeric(x==2 & z==1),
			i_21=as.numeric(x==1 & z==2),
			i_22=as.numeric(x==2 & z==2)
		)
	
	}else{
		stop("Invalid model specification. Must be one of 'NOIA_S', 'NOIA_F', 'G2A', 'F2', 'Finf', 'unweight', 'genotype'.")
	}
	
	return(design.M)
}

