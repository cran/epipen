epi.convert <-
function(modelin,params,modelout=c(modelin,"NOIA_S"),pa=NULL,pb=NULL,distn,mu=NA,ve=NA,signs=rep(1,8),totaleff=NA,verbose=TRUE){

	##########
	# Input/output models:
	# penetrance, logodds, varex, NOIA_S, NOIA_F, G2A, F2, Finf, unweight, genotype
	##########

	#######
	# Convert input model to NOIA S
	# - Identify input format
	# - Validate input
	# - Get NOIA S model params
	#######
	
	
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
	
	
	### Validate distribution
	
	if(!(distn=="logit") && !(distn=="normal") && !(distn=="probit")){
		stop("Distribution must be 'normal', 'logit', or 'probit'.")
	}

	### Validate and convert input model 
	
	if(modelin=="penetrance"){
	
		# vaidate penetrance matrix
		if(!is.numeric(params) || !(ncol(params)==3) || !(nrow(params==3)) || any(params<=0) || any(params>=1)){
			stop("For penetrance models, 'params' must be a 3x3 penetrance matrix with values between 0 and 1")
		}
		
		# restrict distributions
		if(distn=="normal"){
			stop("Penetrance models require 'logit' or 'probit' distribution.")
		}
		
		### convert to NOIA_S
		B.V <- pen.to.noia(pen.M=params, pa=pa.V, pb=pb.V, distn=distn)
		

	}else if(modelin=="logodds"){
	
		### validate log odds
		if(!is.numeric(params) || !(ncol(params)==3) || !(nrow(params)==3)){
			stop("For log odds model, 'params' must be a 3x3 matrix of log odds ratios")
		}
		
		### convert to NOIA_S
		B.V <- logodds.to.noia(logodds.M=params, pa=pa.V, pb=pb.V, distn=distn)
		

	}else if(modelin=="varex"){
	
		### validate variance components
		if(!(length(params)==8) || !is.numeric(params)){
			stop("For variance components, 'params' must be a length 8 vector of proportions of variance explained")
		}else if(any(params<0) || any(params>1) || sum(params)>1 ){
			stop("Variance explained cannot be <0 or >1")
		}
		
		### convert to NOIA_S
		B.V <- varex.to.noia(varex=params,pa=pa.V,pb=pb.V,distn=distn,mu=mu,ve=ve,signs=signs,totaleff=totaleff)
		

	}else if(modelin=="F2"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)
		

	}else if(modelin=="Finf"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)


	}else if(modelin=="G2A"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)	


	}else if(modelin=="unweight"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)
		

	}else if(modelin=="genotype"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)
		
		
	}else if(modelin=="NOIA_F"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		### convert to NOIA_S
		B.V <- glm.to.noia(params,model=modelin,pa=pa.V,pb=pb.V)
		

	}else if(modelin=="NOIA_S"){
	
		### validate betas
		if(!(length(params)==9)){
			stop("params must be a vector of the 9 parameter values for the selected linear model")
		}
		if(!is.numeric(params)){
			stop("The parameter values in 'params' must be numeric")
		}
		if(any(is.na(params))){
			stop("All parameter values must be specified in 'params'")
		}
		
		# no conversion needed
		B.V <- params
		names(B.V) <- c("mu","alpha_a","delta_a",
			"alpha_b","delta_b","alphaalpha_ab","deltaalpha_ab",
			"alphadelta_ab","deltadelta_ab")
		
		
	}else{
		stop("Invalid input model specification. Must be one of 'penetrance', 'logodds', 'varex', 'NOIA_S', 'NOIA_F', 'G2A', 'F2', 'Finf', 'unweight', 'genotype'.")
	}
	
	
	#######
	# Get desired output model(s)
	#######
	
	out <- list(NULL)
	
	### Input design matrix
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
	
	for(i in 1:length(modelout)){
	
		if(modelout[i]=="NOIA_S"){
		
			# already converted to NOIA_S as reference model
			out[[i]] <- B.V
		
		}else if(modelout[i]=="NOIA_F"){
		
			### output design matrix
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
			
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("R","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="G2A"){
		
			### output design matrix
			# allele freqs
			pa1 <- pa.V[1] + pa.V[2]/2
			pb1 <- pb.V[1] + pb.V[2]/2
		
			# univariate design matrices
			designA.M <- matrix(c(
				1,-2*(1-pa1),-2*(1-pa1)^2,
				1,-1+2*pa1,2*pa1*(1-pa1),
				1,2*pa1,-2*pa1^2),3,3,byrow=TRUE)
		
			designB.M <- matrix(c(
				1,-2*(1-pb1),-2*(1-pb1)^2,
				1,-1+2*pb1,2*pb1*(1-pb1),
				1,2*pb1,-2*pb1^2),3,3,byrow=TRUE)
		
			# combine two loci in design
			# rearrange columns for nicer presentation
			d1.M <- kronecker(designB.M,designA.M)
			d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
			
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="F2"){
		
			### output design matrix
			designA.M <- matrix(c(
				1,-1,-1/2,
				1,0,1/2,
				1,1,-1/2),3,3,byrow=TRUE)
		
			designB.M <- designA.M
		
			# combine two loci in design
			# rearrange columns for nicer presentation
			d1.M <- kronecker(designB.M,designA.M)
			d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
			
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="Finf"){
		
			### output design matrix
			designA.M <- matrix(c(
				1,-1,0,
				1,0,1,
				1,1,0),3,3,byrow=TRUE)
		
			designB.M <- designA.M
		
			# combine two loci in design
			# rearrange columns for nicer presentation
			d1.M <- kronecker(designB.M,designA.M)
			d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
			
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="unweight"){
		
			### output design matrix
			designA.M <- matrix(c(
				1,-1,-1/3,
				1,0,2/3,
				1,1,-1/3),3,3,byrow=TRUE)
		
			designB.M <- designA.M
		
			# combine two loci in design
			# rearrange columns for nicer presentation
			d1.M <- kronecker(designB.M,designA.M)
			d1.M <- cbind(d1.M[,1:4],d1.M[,7],d1.M[,5:6],d1.M[,8:9])
			
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("mu","a_a","d_a",
				"a_b","d_b","aa_ab","da_ab",
				"ad_ab","dd_ab")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="genotype"){
		
			### output design matrix
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
		
			# compute coefs
			b2.V <- as.vector(solve(d1.M) %*% design.M %*% as.matrix(B.V,9,1))
			names(b2.V) <- c("mu","gam_1","gam_2",
				"beta_1","beta_2","i_11","i_12",
				"i_21","i_22")
				
			# save
			out[[i]] <- b2.V
		
		}else if(modelout[i]=="penetrance"){
		
			if(distn=="logit"){
			
				# compute probabilities
				EXB <- design.M %*% matrix(B.V,nrow=9)
				Ey <- exp(EXB)/(1+exp(EXB))
				
				# format and save
				pen <- matrix(Ey,3,3,byrow=TRUE)
				colnames(pen) <- c("AA","Aa","aa")
				rownames(pen) <- c("BB","Bb","bb")
				out[[i]] <- pen
				
			}else if(distn=="probit"){
			
				# compute probabilities
				EXB <- design.M %*% matrix(B.V,nrow=9)
				Ey <- pnorm(EXB)
				
				# format and save
				pen <- matrix(Ey,3,3,byrow=TRUE)
				colnames(pen) <- c("AA","Aa","aa")
				rownames(pen) <- c("BB","Bb","bb")
				out[[i]] <- pen
			
			}else{
				stop("Conversion to penetrance model only possible for 'logit' or 'probit' distribution")
			}
		}else if(modelout[i]=="logodds"){
		
			if(distn=="logit"){
			
				# compute
				EXB <- design.M %*% matrix(B.V,nrow=9)
				
				# format and save
				lo <- matrix(EXB,3,3,byrow=TRUE)
				colnames(lo) <- c("AA","Aa","aa")
				rownames(lo) <- c("BB","Bb","bb")
				out[[i]] <- lo
				
			}else if(distn=="probit"){
			
				# compute probabilities
				EXB <- design.M %*% matrix(B.V,nrow=9)
				Ey <- pnorm(EXB)
				lo <- log(Ey/(1-Ey))
				
				# format and save
				pen <- matrix(lo,3,3,byrow=TRUE)
				colnames(pen) <- c("AA","Aa","aa")
				rownames(pen) <- c("BB","Bb","bb")
				out[[i]] <- pen
			}else{
				stop("Conversion to log odds model only possible for 'logit' or 'probit' distribution")
			}
		
		}else if(modelout[i]=="varex"){
		
			### variance of design variables
			VaddA <- pa.V[2]*(1-pa.V[2])+4*pa.V[3]*(1-pa.V[2])-4*pa.V[3]^2
			VaddB <- pb.V[2]*(1-pb.V[2])+4*pb.V[3]*(1-pb.V[2])-4*pb.V[3]^2
			VdomA <- 4*pa.V[3]*pa.V[2]*pa.V[1]/VaddA
			VdomB <- 4*pb.V[3]*pb.V[2]*pb.V[1]/VaddB
			
			### variance components
			varcomp.V <- c(VaddA*B.V[2]^2,VdomA*B.V[3]^2,VaddB*B.V[4]^2,VdomB*B.V[5]^2,VaddA*VaddB*B.V[6]^2,VdomA*VaddB*B.V[7]^2,VaddA*VdomB*B.V[8]^2,VdomA*VdomB*B.V[9]^2)
			
			# error variance
			if(distn=="normal"){
				if(is.na(ve)){
					vy <- 1
					ve <- 1-sum(varcomp.V)
				}else{
					vy <- sum(varcomp.V)+ve
				}
			}else if(distn=="logit"){
				if(!is.na(ve)){
					warning("ve is ignored for logit model")
				}
				ve <- (pi^2)/3
				vy <- sum(varcomp.V)+ve
			}else if(distn=="probit"){
				if(!is.na(ve)){
					warning("ve is ignored for probit model")
				}
				ve <- 1
				vy <- sum(varcomp.V)+ve
			}
			
			# proportion of total variance
			varex.V <- varcomp.V/vy
			names(varex.V) <- c("AddA","DomA","AddB","DomB","AddAxAddB","DomAxAddB","AddAxDomB","DomAxDomB")
			out[[i]] <- varex.V
		
		}else{
			stop(paste("Invalid output model: '",modelout[i],"'. Must be one of 'penetrance', 'logodds', 'varex', 'NOIA_S', 'NOIA_F', 'G2A', 'F2', 'Finf', 'unweight', 'genotype'.",sep=""))
		}
		
		names(out)[i] <- modelout[i]
	
	}
	
	#######
	# Extra info if desired (verbose=T)
	#######

	if(verbose){
		out[[i+1]] <- epi.info(betas=B.V,pa=pa.V,pb=pb.V,model="NOIA_S",distn=distn,ve=ve)
		names(out)[i+1] <- "Additional_info"
	}
	
	#######
	# Finish
	#######
	
	return(out)
	
}

# eof