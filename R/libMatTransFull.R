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
	   ifelse(Mu == "G", 2, 0))

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
               ifelse(Sigma == "VVV", 14,0))))))))))))))
  
  Psi.num <- ifelse(Psi == "UI", 0,
             ifelse(Psi == "EI", 1,
             ifelse(Psi == "VI", 2,
             ifelse(Psi == "EE", 3,
             ifelse(Psi == "VE", 4,
             ifelse(Psi == "EV", 5,
             ifelse(Psi == "VV", 6, -1)))))))
  
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
  
  Psi <- ifelse(Psi.num == 0, "UI", 
         ifelse(Psi.num == 1, "EI", 
         ifelse(Psi.num == 2, "VI", 
         ifelse(Psi.num == 3, "EE", 
         ifelse(Psi.num == 4, "VE", 
         ifelse(Psi.num == 5, "EV", 
         ifelse(Psi.num == 6, "VV", NULL)))))))
  
  list(Mu = Mu, Sigma = Sigma, Psi = Psi, Mu.num = Mu.num, Sigma.num = Sigma.num, Psi.num = Psi.num)
}





MatTrans.EM <- function(Y, initial = NULL, id = NULL, la = NULL, nu = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, model = NULL, trans = "Power", la.type = 0, tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE){

	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong dimensionality T...\n")
	if(is.null(la)){
		la <- matrix(0.0, K, p)
	}
	if(is.null(nu)){
		nu <- matrix(0.0, K, T)
	}

	if(!is.null(initial)){
		
		if(length(initial) < 1) stop("Wrong initialization...\n")

		K <- length(initial[[1]]$tau)

		if(is.null(model)){
			index <- matrix(NA, 3, 14*7*2)
			model <- rep(NA, 14*7*2)
			iter <- 0
			for(Mu.type in 1:2){
				for(Sigma.type in 1:14){
					for(Psi.type in 0:6){
																	iter <- iter + 1
						model[iter] <- paste(code.back(c(Mu.type, Sigma.type, Psi.type))$Mu, "-", code.back(c(Mu.type, Sigma.type, Psi.type))$Sigma, "-", code.back(c(Mu.type, Sigma.type, Psi.type))$Psi, sep="")
						

						
						index[,iter] <- c(Mu.type, Sigma.type, Psi.type)
					}
				}
			}
		}
		else{
			index <- matrix(NA, 3, length(model))
			for(iter in 1:length(model)){
				code1 <- code.convert(model[iter])$Mu.num
				code2 <- code.convert(model[iter])$Sigma.num
				code3<- code.convert(model[iter])$Psi.num
				index[,iter] <- c(code1, code2, code3)
			}
		}



		loglik <- rep(-Inf, length(model))
		bic <- rep(Inf, length(model)) 
		result <- list()


		for(i in 1:length(initial)){
			if(silent != TRUE){
				cat("Initialization", i, "\n")
			}
	
			ll <- rep(0, 3)
			misc_double <- c(tol, 0.0, 0.0)
			conv <- rep(0, 3)
			id <- rep(0, n)


			for(iter in 1:dim(index)[2]){

				
				if(trans == "Power"){trans.type <- 1}
				else if(trans == "Manly"){trans.type <- 2}
				#cat("trans.type", trans.type, "\n")
				misc_int <- c(p, T, n, K, max.iter, index[1,iter], index[2,iter], index[3,iter], la.type, trans.type)

				try0 <- try(temp <- .C("run_Mat_Trans_Full", y = as.double(initial[[i]]$y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(initial[[i]]$tau), la1 = as.double(as.vector(la)), nu1 = as.double(as.vector(nu)), Mu1 = as.double(initial[[i]]$Mu1), invS1 = as.double(initial[[i]]$invS1), invPsi1 = as.double(initial[[i]]$invPsi1), detS = as.double(initial[[i]]$detS), detPsi = as.double(initial[[i]]$detPsi), gamma1 = as.double(initial[[i]]$gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), scale = as.double(initial[[i]]$scale), PACKAGE = "MatTransMix"), silent = TRUE)


				try1 <- try(invS <- array(temp$invS1, dim = c(p, p, K)))
				try2 <- try(invPsi <- array(temp$invPsi1, dim = c(T, T, K)))
				
				try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
				try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

				if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){
			
					if(!is.na(temp$ll[1])){		

						if ((temp$ll[1] > loglik[iter]) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(temp$id) > size.control) && all(table(temp$id) < n-size.control)){
							loglik[iter] <- temp$ll[1]	
							bic[iter] <- temp$ll[2]	

							if(silent != TRUE){
								cat("Model", model[iter], "BIC", bic[iter], "\n")
							}
											
							r <- NULL
						
							r$tau <- temp$tau
							r$la <- matrix(temp$la1, nrow = K)
							r$nu <- matrix(temp$nu1, nrow = K)
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
							result[[iter]] <- r

						}
					}

				}

			}


		}				

		find.min <- which.min(bic)
		best.result <- result[find.min]
		best.model <- model[find.min]
		best.loglik <- loglik[find.min]
		best.bic <- bic[find.min]


		ret <- list(result = result, model = model, loglik = loglik, bic = bic, best.result = best.result, best.model = best.model, best.loglik = best.loglik, best.bic = best.bic)

	
		class(ret) <- "MatTransMix"
	

		return(ret)

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

