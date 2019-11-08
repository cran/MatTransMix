options(show.error.messages = FALSE)
eucl.dist <- function(X, Y){

	sqrt(sum(X - Y)^2)

}

Trans.bic <- function(logl, n, M){
	return(-2 * logl + M * log(n))
}

# EM algorithm
MatTrans.init <- function(Y, K, n.start = 10, scale = 1){

	Y <- Y/scale	
	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong dimensionality T...\n")


	init <- list()
	i <- 0
	
	Psi.inv <- array(NA, c(T, T, K))
	Sigma.inv <- array(NA, c(p, p, K))
	detPsi <- rep(NA, K)
	detS <- rep(NA, K)

	
	repeat{

		s <- sample(1:n,K)
		
		mat.Y <- t(apply(Y, 3, as.matrix))
		if((p>1) || (T>1)){
			centers <- mat.Y[s,]
		}
		else{
			centers <- mat.Y[s]
		}
		D <- NULL
		if(K == 1) {D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers))}
		else{
			if((p>1) || (T>1)){
				for (k in 1:K) D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers[k,]))
			}
			else{
				for (k in 1:K) D <- cbind(D, t(mat.Y - centers[k])^2)
			}
		}
		PI <- as.matrix(D == apply(D, 1, min)) * 1
		W.result <- NULL

		iif <- 0
		for (k in 1:K){
			index <- PI[,k] == 1
			if((p>1) && (T>1)){
				var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
				var.est <- var(var.est)
				A <- try(Psi.inv[,,k] <- solve(var.est))
				detPsi[k] <- 1/det(Psi.inv[,,k])
				var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
				var.est <- var(var.est)
				B <- try(Sigma.inv[,,k] <- solve(var.est))
				detS[k] <- 1/det(Sigma.inv[,,k])

			}
			else if((p == 1) && (T > 1)){
				var.est <- var(t(Y[,,index]))
				A <- try(Psi.inv[,,k] <- solve(var.est))
				detPsi[k] <- 1/(Psi.inv[,,k])
				var.est <- var(as.vector(Y[,,index]))
				B <- try(Sigma.inv[,,k] <- solve(var.est))
				detS[k] <- 1/det(Sigma.inv[,,k])
			}
			else if(T == 1){
				var.est <- var(as.vector(Y[,,index]))
				A <- try(Psi.inv[,,k] <- solve(var.est))
				detPsi[k] <- 1/(Psi.inv[,,k])
				var.est <- var(t(Y[,,index]))
				B <- try(Sigma.inv[,,k] <- solve(var.est))
				detS[k] <- 1/det(Sigma.inv[,,k])

			}
			if((class(A) == "try-error") || (class(B) == "try-error")){iif <- 1}	
		}
		if(iif == 0){
			W.result$y <- as.vector(Y)				
			W.result$gamma1 <- PI
			W.result$invS1 <- Sigma.inv
			W.result$tau <- rep(0, K)
			W.result$Mu1 <- rep(0, p*T*K)
			W.result$invPsi1 <- Psi.inv
			W.result$detS <- detS
			W.result$detPsi <- detPsi
			W.result$scale <- scale				
	
			i <- i + 1
			init[[i]] <- W.result	

		}
		if (i == n.start) break

	}
	return(init)

}



code.convert <- function(code){
  Mu <- substr(code,1,1)
  Sigma <- substr(code,3,5)
  Psi <- substr(code,7,8)
  
  Mu.num <- ifelse(Mu == "A", 1,
	   ifelse(Mu == "G", 2, 
	   ifelse(Mu == "X", 0, -1)))

  Sigma.num <- ifelse(Sigma == "EII", 1,
               ifelse(Sigma == "VII", 2,
               ifelse(Sigma == "EEI", 3,
               ifelse(Sigma == "VEI", 4,
               ifelse(Sigma == "EVI", 5,
               ifelse(Sigma == "VVI", 6,
               ifelse(Sigma == "EEE", 7,
               ifelse(Sigma == "VEE", 8,
               ifelse(Sigma == "EVE", 9,
               ifelse(Sigma == "VVE", 10,
               ifelse(Sigma == "EEV", 11,
               ifelse(Sigma == "VEV", 12,
               ifelse(Sigma == "EVV", 13,
               ifelse(Sigma == "VVV", 14,
 	       ifelse(Sigma == "XXX", 0, -1)))))))))))))))
  
  Psi.num <- ifelse(Psi == "UI", 1,
             ifelse(Psi == "EI", 2,
             ifelse(Psi == "VI", 3,
             ifelse(Psi == "EE", 4,
             ifelse(Psi == "VE", 5,
             ifelse(Psi == "EV", 6,
             ifelse(Psi == "VV", 7, 
	     ifelse(Psi == "XX", 0, -1))))))))
  
  list(Mu = Mu, Sigma = Sigma, Psi = Psi, Mu.num = Mu.num, Sigma.num = Sigma.num, Psi.num = Psi.num)
}


code.back <- function(num){


  Mu.num <- num[1]
  Sigma.num <- num[2]
  Psi.num <- num[3]
 
  Mu <- ifelse(Mu.num == 1, "A",
	ifelse(Mu.num == 2, "G", NULL))
	
 
  Sigma <- ifelse(Sigma.num == 1, "EII",
           ifelse(Sigma.num == 2, "VII", 
           ifelse(Sigma.num == 3, "EEI", 
           ifelse(Sigma.num == 4, "VEI", 
           ifelse(Sigma.num == 5, "EVI", 
           ifelse(Sigma.num == 6, "VVI", 
           ifelse(Sigma.num == 7, "EEE", 
           ifelse(Sigma.num == 8, "VEE", 
           ifelse(Sigma.num == 9, "EVE", 
           ifelse(Sigma.num == 10, "VVE", 
           ifelse(Sigma.num == 11, "EEV", 
           ifelse(Sigma.num == 12, "VEV", 
           ifelse(Sigma.num == 13, "EVV", 
           ifelse(Sigma.num == 14, "VVV", NULL))))))))))))))
  
  Psi <- ifelse(Psi.num == 1, "UI", 
         ifelse(Psi.num == 2, "EI", 
         ifelse(Psi.num == 3, "VI", 
         ifelse(Psi.num == 4, "EE", 
         ifelse(Psi.num == 5, "VE", 
         ifelse(Psi.num == 6, "EV", 
         ifelse(Psi.num == 7, "VV", NULL)))))))
  
  list(Mu = Mu, Sigma = Sigma, Psi = Psi, Mu.num = Mu.num, Sigma.num = Sigma.num, Psi.num = Psi.num)
}





MatTrans.EM <- function(Y, initial = NULL, la = NULL, nu = NULL, model = NULL, trans = "Gaussian", la.type = 0, tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE){

	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]
	

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong dimensionality T...\n")

	if(!is.null(initial)){
		
		if(length(initial) < 1) stop("Wrong initialization...\n")

		K <- length(initial[[1]]$tau)
		if(is.null(la) && (trans != "Gaussian")){
			la <- matrix(0.5, K, p)
			cat("Initial lambda -- 0.5 \n")

		}
		if(is.null(nu) && (trans != "Gaussian")){
			nu <- matrix(0.5, K, T)
			cat("Initial nu -- 0.5 \n")
		}

		if((la.type == 0) && (trans != "Gaussian")){
			cat("Unrestricted lambda type \n")
		}
		else if((la.type == 1) && (trans != "Gaussian")){
			cat("Lambda same across all variables \n")
		}



		if(is.null(model)){
			index <- matrix(NA, 3, 14*7*2)
			Model <- rep(NA, 14*7*2)
			iter <- 0
			for(Mu.type in 1:2){
				for(Sigma.type in 1:14){
					for(Psi.type in 1:7){
																	iter <- iter + 1
						Model[iter] <- paste(code.back(c(Mu.type, Sigma.type, Psi.type))$Mu, "-", code.back(c(Mu.type, Sigma.type, Psi.type))$Sigma, "-", code.back(c(Mu.type, Sigma.type, Psi.type))$Psi, sep="")
								
						index[,iter] <- c(Mu.type, Sigma.type, Psi.type)
					}
				}
			}
		}

		else{
			
			index <- NULL
			Model <- NULL
			for(ij in 1:length(model)){
				
				code1 <- code.convert(model[ij])$Mu.num
				code2 <- code.convert(model[ij])$Sigma.num
				code3<- code.convert(model[ij])$Psi.num
			
				
				index.temp <- c(code1, code2, code3)
				
				#cat(index.temp, "\n")

 				if(any(index.temp == -1)){
					stop("model code is not identifiable...\n")
				}
				else if(any(index.temp == 0)){

					if(code1 == 0){
						code1 <- c(1,2)
					}
					if(code2 == 0){
						code2 <- seq(1,14)
					}
					if(code3 == 0){
						code3 <- seq(1,7)
					}

					
					index.temp <- t(expand.grid(x = code1, y = code2, z = code3))
					for(ij2 in 1:dim(index.temp)[2]){
						Model <- c(Model, paste(code.back(index.temp[,ij2])$Mu, "-", code.back(index.temp[,ij2])$Sigma, "-", code.back(index.temp[,ij2])$Psi, sep=""))
					}
						

				}
				else{
					Model <- c(Model, paste(code.back(index.temp)$Mu, "-", code.back(index.temp)$Sigma, "-", code.back(index.temp)$Psi, sep=""))
				}

				index <- cbind(index, index.temp)

				#cat(index, "\n")

			}

			
		}



		loglik <- rep(-Inf, dim(index)[2])
		bic <- rep(Inf, dim(index)[2]) 
		result <- list()
				
		if(trans == "Power"){
			trans.type <- 1
			cat("Power transformation \n")
		}
		else if(trans == "Manly"){

			trans.type <- 2
			cat("Manly transformation \n")

		}
		else if(trans == "Gaussian"){
			trans.type <- 0
			cat("Gaussian -- no transformation \n")

		}



		for(i in 1:length(initial)){
			if(silent != TRUE){
				cat("Initialization", i, "\n")
			}
	
			ll <- rep(0, 3)
			misc_double <- c(tol, 0.0, 0.0)
			conv <- rep(0, 3)
			id <- rep(0, n)

			r <- NULL

			for(iter in 1:dim(index)[2]){

				
				if(trans == "Power"){
					trans.type <- 1
					
				}
				else if(trans == "Manly"){

					trans.type <- 2
					

				}
				else if(trans == "Gaussian"){
					trans.type <- 0
					
				}

				#cat("trans.type", trans.type, "\n")
				misc_int <- c(p, T, n, K, max.iter, index[1,iter], index[2,iter], index[3,iter], la.type, trans.type)

				try0 <- try(temp <- .C("run_Mat_Trans_Full", y = as.double(initial[[i]]$y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(initial[[i]]$tau), la1 = as.double(as.vector(la)), nu1 = as.double(as.vector(nu)), Mu1 = as.double(initial[[i]]$Mu1), invS1 = as.double(initial[[i]]$invS1), invPsi1 = as.double(initial[[i]]$invPsi1), detS = as.double(initial[[i]]$detS), detPsi = as.double(initial[[i]]$detPsi), gamma1 = as.double(initial[[i]]$gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), scale = as.double(initial[[i]]$scale), PACKAGE = "MatTransMix"), silent = TRUE)

				 

				try1 <- try(invS <- array(temp$invS1, dim = c(p, p, K)))
				try2 <- try(invPsi <- array(temp$invPsi1, dim = c(T, T, K)))
				
				try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
				try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

				if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){

			
					if(!is.na(temp$ll[1])){	
						
						if ((temp$ll[1] > loglik[iter]) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(temp$id) > size.control) && (length(table(temp$id))==K)){


							loglik[iter] <- temp$ll[1]	
							bic[iter] <- temp$ll[2]	

							if(silent != TRUE){
								cat("Model", Model[iter], "BIC", bic[iter], "\n")
							}
											
							r <- NULL
						
							r$tau <- temp$tau
							if(trans.type != 0){
								r$la <- matrix(temp$la1, nrow = K)
								r$nu <- matrix(temp$nu1, nrow = K)
								r$LA <- array(NA, dim = c(p, T, K))
								for(k in 1:K){
									for(j in 1:p){
										for(t in 1:T){
											r$LA[j,t,k] <- r$la[k,j] + r$nu[k,t]
										}
									}
								}
							}
							else{
								r$la <- NA
								r$nu <- NA
								r$LA <- NA

							}
							r$Sigma <- Sigma
							r$Psi <- Psi
							r$detS <- temp$detS
							r$detPsi <- temp$detPsi
							r$Mu <- array(temp$Mu1, dim = c(p, T, K))
							r$gamma <- matrix(temp$gamma1, nrow = n)
							r$iter <- temp$conv[1]
							r$pars <- temp$conv[3]

							r$id <- temp$id
							r$flag <- temp$conv[2]
							r$ll <- temp$ll[1]
							r$bic <- temp$ll[2]
							
							result[[iter]] <- r
							

						}
					}

				}

			}


		}				


		find.min <- which.min(bic)
		best.result <- result[find.min]
		best.model <- Model[find.min]
		best.loglik <- loglik[find.min]
		best.bic <- bic[find.min]


		
		ret <- list(result = result, model = Model, loglik = loglik, bic = bic, best.result = best.result, best.model = best.model, best.loglik = best.loglik, best.bic = best.bic, trans = trans)
	
		class(ret) <- "MatTransMix"
	

		return(ret)

	}
	else{
		stop("Use MatTrans.init() to get initialization...\n")

	}


}






check.Sigma <- function(A){
	
	loglik <- array(A$loglik, dim= c(7, 14, 2))
	
	
	for(i in 1:2){
		for(t in 1:7){
			compare <- loglik[t,,i]	
			if(compare[2]<compare[1]){
				cat("VII>EII", compare[2], compare[1], "Mu", i, "Psi", t, "\n")
			}
			if(compare[3]<compare[1]){
				cat("EEI>EII", compare[3], compare[1], "Mu", i, "Psi", t, "\n")
			}
			if(compare[4]<compare[2]){
				cat("VEI>VII", compare[4], compare[2], "Mu", i, "Psi", t, "\n")
			}
			if(compare[4]<compare[3]){
				cat("VEI>EEI", compare[4], compare[3], "Mu", i, "Psi", t, "\n")
			}
			if(compare[5]<compare[3]){
				cat("EVI>EEI", compare[5], compare[3], "Mu", i, "Psi", t, "\n")
			}
			if(compare[6]<compare[5]){
				cat("VVI>EVI", compare[6], compare[5], "Mu", i, "Psi", t, "\n")
			}
			if(compare[6]<compare[4]){
				cat("VVI>VEI", compare[6], compare[4], "Mu", i, "Psi", t, "\n")
			}
			if(compare[7]<compare[3]){
				cat("EEE>EEI", compare[7], compare[3], "Mu", i, "Psi", t, "\n")
			}
			if(compare[8]<compare[7]){
				cat("VEE>EEE", compare[8], compare[7], "Mu", i, "Psi", t, "\n")
			}
			if(compare[8]<compare[4]){
				cat("VEE>VEI", compare[8], compare[4], "Mu", i, "Psi", t, "\n")
			}
			if(compare[9]<compare[7]){
				cat("EVE>EEE", compare[9], compare[7], "Mu", i, "Psi", t, "\n")
			}
			if(compare[9]<compare[5]){
				cat("EVE>EVI", compare[9], compare[5], "Mu", i, "Psi", t, "\n")
			}
			if(compare[10]<compare[8]){
				cat("VVE>VEE", compare[10], compare[8], "Mu", i, "Psi", t, "\n")
			}
			if(compare[10]<compare[9]){
				cat("VVE>EVE", compare[10], compare[9], "Mu", i, "Psi", t, "\n")
			}
			if(compare[10]<compare[6]){
				cat("VVE>VVI", compare[10], compare[6], "Mu", i, "Psi", t, "\n")
			}
			if(compare[11]<compare[7]){
				cat("EEV>EEE", compare[11], compare[7], "Mu", i, "Psi", t, "\n")
			}
			if(compare[12]<compare[8]){
				cat("VEV>VEE", compare[12], compare[8], "Mu", i, "Psi", t, "\n")
			}
			if(compare[12]<compare[11]){
				cat("VEV>EEV", compare[12], compare[11], "Mu", i, "Psi", t, "\n")
			}
			if(compare[13]<compare[9]){
				cat("EVV>EVE", compare[13], compare[9], "Mu", i, "Psi", t, "\n")
			}
			if(compare[13]<compare[11]){
				cat("EVV>EEV", compare[13], compare[11], "Mu", i, "Psi", t, "\n")
			}
			if(compare[14]<compare[13]){
				cat("VVV>EVV", compare[14], compare[13], "Mu", i, "Psi", t, "\n")
			}
			if(compare[14]<compare[12]){
				cat("VVV>VEV", compare[14], compare[12], "Mu", i, "Psi", t, "\n")
			}
			if(compare[14]<compare[10]){
				cat("VVV>VVE", compare[14], compare[10], "Mu", i, "Psi", t, "\n")

			}
		}
	}
	

}

check.Psi <- function(A){
	
	loglik <- array(A$loglik, dim= c(7, 14, 2))
	
	
	for(i in 1:2){
		for(j in 1:14){
			compare <- loglik[,j,i]	
			if(compare[2]<compare[1]){
				cat("EI>UI", compare[2], compare[1], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[3]<compare[2]){
				cat("VI>EI", compare[3], compare[2], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[4] <compare[2]){
				cat("EE>EI", compare[4], compare[2], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[5]<compare[3]){
				cat("VE>VI", compare[5], compare[3], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[5]<compare[4]){
				cat("VE>EE", compare[5], compare[4], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[6]<compare[4]){
				cat("EV>EE", compare[6], compare[4], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[7]<compare[6]){
				cat("VV>EV", compare[7], compare[6], "Mu", i, "Sigma", j, "\n")
			}
			if(compare[7]<compare[5]){
				cat("VV>VE", compare[7], compare[5], "Mu", i, "Sigma", j, "\n")
			}
			
		}
	}


}

check.Mu <- function(A){
	
	loglik <- array(A$loglik, dim= c(7, 14, 2))
	
	
	for(t in 1:7){
		for(j in 1:14){
			compare <- loglik[t,j,]	
			if(compare[2]<compare[1]){
				cat("G>A", compare[2], compare[1], "Sigma", j, "Psi", t, "\n")
			}	
		}
	}

	

}

