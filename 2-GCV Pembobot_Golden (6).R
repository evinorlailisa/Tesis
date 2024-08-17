pembobot <- function(data, coord, param) {
  start.time <- Sys.time()
  maxit <- 10
  eps <- 5
  
  data <- data.frame(data)
  n <- nrow(data)
  y1 <- as.matrix(data[, 1])
  y2 <- as.matrix(data[, 2])
  x <- as.matrix(data[, -c(1, 2)])
  x <- cbind(1, x)
  
  # Menghitung Jarak Euclidean
  U <- coord[, 1]
  V <- coord[, 2]
  n <- nrow(x)
  p <- ncol(x)
  d <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      d[i, j] <- sqrt(((U[i] - U[j]) ^ 2) + ((V[i] - V[j]) ^ 2))
    }
  }
  d.max <- max(d)
  Wlist = list()
  
  ###### Mencari pembobot optimum dengan GCV
  ### Bandwidth Gaussian Kernel
  h_gaussian <- rep(0, n)
  GCVmin <- rep(0, n)
  
  for (i in 1:n) {
    A <- 0.0001
    B <- d.max
    iter_gaussian <- 0
    minGCV <- 0
    selisih <- 1000
    
    while ((selisih > 0.0001) && (iter_gaussian <= 1000)) {
      h_awal <- seq(A, B, by = (B - A) / 10)
      nh <- length(h_awal)  # banyaknya awalan bandwidth
      GCV <- matrix(0, 1, nh)  # membuat matrix untuk tempat nilai GCV
      colnames(GCV) <- c(h_awal)
      Wb_gaus <- matrix(0, n, n)
      
      for (k in 1:nh) {
        h <- h_awal[k]
        for (ii in 1:n) {
          for (jj in 1:n) {
            # Rumus bandwidth Gaussian kernel
            Wb_gaus[ii, jj] <- exp(-0.5 * (d[ii, jj] / h) ^ 2)
          }
        }
        
        W <- Wb_gaus
        x.gcv <- as.matrix(x[-i, ])  # menghilangkan data ke-i dari estimasi
        y1.gcv <- as.matrix(y1[-i, ])  # menghilangkan data y1 ke-i dari estimasi
        y2.gcv <- as.matrix(y2[-i, ])  # menghilangkan data y2 ke-i dari estimasi
        
        beta.gcv <- as.matrix(pop(y1.gcv, y2.gcv, x.gcv, W, param, i, maxit, eps)$param)
        y1hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[1:p]))  # menghitung y1_hat
        y2hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[(p + 1):(2 * p)]))  # menghitung y2_hat
        GCV[k] <- (n * ((abs(y1[i] - y1hat.gcv) + abs(y2[i] - y2hat.gcv)) / (n - length(beta.gcv)) ^ 2))  # menghitung GCV
      }
      hasilGCV_gaussian <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
      A0 <- A
      B0 <- B
      minGCV <- min(GCV)
      l_min <- which(GCV == minGCV)[1]
      if (l_min == 1) {
        A <- h_awal[l_min]
        B <- h_awal[l_min + 1]
      } else if (l_min == nh) {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min]
      } else {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min + 1]
      }
      selisih <- (B0 - A0) - (B - A)
      iter_gaussian <- iter_gaussian + 1
      cat('Gaussian Kernel (kota ke : ', i, ', iterasi ke : ', iter_gaussian, ', selisih : ', selisih, '\n')
    }
    hasilGCV_gaussian <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
    hasilGCV <- hasilGCV_gaussian[order(hasilGCV_gaussian[, 1]), ]
    GCVmin[i] <- hasilGCV[1, 1]
    h_gaussian[i] <- hasilGCV[1, 1]
  }
  hasilGCV_gaussian <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
  bestGCV_gaussian <- data.frame(h = 0, GCV = sum(GCVmin))  # mengambil bandwidth dengan nilai GCV paling minimum
  
  best <- data.frame(Kernel = c('Gaussian Kernel'), rbind(bestGCV_gaussian))
  colnames(best) <- c('Kernel', 'Bandwidth', 'GCV')
  best.h <- best[1, 2]
  
  Wb <- matrix(0, n, n)
  h <- h_gaussian
  for (i in 1:n) {
    for (j in 1:n) {
      # Rumus bandwidth Gaussian kernel
      Wb[i, j] <- exp(-0.5 * (d[i, j] / h[i]) ^ 2)
    }
  }
  W <- Wb
  best.h <- 'Adaptive Weight'
  
  end.time <- Sys.time()
  cat('Hasil iterasi pembobot :', '\n')
  cat(paste('Best Kernel = ', best[1, 1], ' --- Bandwidth = ', best.h, ' --- GCV = ', best[1, 3]), '\n',
      'Processing time : ', paste(round(end.time - start.time, 2), ' mins'))
  Wlist[[1]] = W
  cat('\n')
  cat('\n')
  
  
  ### Bandwidth Exponential Kernel
  start.time <- Sys.time()
  h_expo <- rep(0, n)
  GCVmin <- rep(0, n)
  
  for (i in 1:n) {
    A <- 0.0001
    B <- d.max
    iter_expo <- 0
    minGCV <- 0
    selisih <- 1000
    
    while ((selisih > 0.0001) && (iter_expo <= 1000)) {
      h_awal <- seq(A, B, by = (B - A) / 10)
      nh <- length(h_awal)  # banyaknya awalan bandwidth
      GCV <- matrix(0, 1, nh)  # membuat matrix untuk tempat nilai GCV
      colnames(GCV) <- c(h_awal)
      Wb_expo <- matrix(0, n, n)
      
      for (k in 1:nh) {
        h <- h_awal[k]
        for (ii in 1:n) {
          for (jj in 1:n) {
            # Rumus bandwidth Eksponensial kernel
            Wb_expo[ii, jj] <- sqrt(exp(-(d[ii, jj] / h) ^ 2))
          }
        }
        
        W <- Wb_expo
        x.gcv <- as.matrix(x[-i, ])  # menghilangkan data ke-i dari estimasi
        y1.gcv <- as.matrix(y1[-i, ])  # menghilangkan data y1 ke-i dari estimasi
        y2.gcv <- as.matrix(y2[-i, ])  # menghilangkan data y2 ke-i dari estimasi
        
        beta.gcv <- as.matrix(pop(y1.gcv, y2.gcv, x.gcv, W, param, i, maxit, eps)$param)
        y1hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[1:p]))  # menghitung y1_hat
        y2hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[(p + 1):(2 * p)]))  # menghitung y2_hat
        GCV[k] <- (n * ((abs(y1[i] - y1hat.gcv) + abs(y2[i] - y2hat.gcv)) / (n - length(beta.gcv)) ^ 2))  # menghitung GCV
      }
      hasilGCV_expo <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
      A0 <- A
      B0 <- B
      minGCV <- min(GCV)
      l_min <- which(GCV == minGCV)[1]
      if (l_min == 1) {
        A <- h_awal[l_min]
        B <- h_awal[l_min + 1]
      } else if (l_min == nh) {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min]
      } else {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min + 1]
      }
      selisih <- (B0 - A0) - (B - A)
      iter_expo <- iter_expo + 1
      cat('Exponential Kernel (kota ke : ', i, ', iterasi ke : ', iter_expo, ', selisih : ', selisih, '\n')
    }
    hasilGCV_expo <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
    hasilGCV <- hasilGCV_expo[order(hasilGCV_expo[, 1]), ]
    GCVmin[i] <- hasilGCV[1, 1]
    h_expo[i] <- hasilGCV[1, 1]
  }
  hasilGCV_expo <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
  bestGCV_expo <- data.frame(h = 0, GCV = sum(GCVmin))  # mengambil bandwidth dengan nilai GCV paling minimum
  
  best <- data.frame(Kernel = c('Exponential Kernel'), rbind(bestGCV_expo))
  colnames(best) <- c('Kernel', 'Bandwidth', 'GCV')
  best.h <- best[1, 2]
  
  Wb <- matrix(0, n, n)
  h <- h_expo
  for (i in 1:n) {
    for (j in 1:n) {
      # Rumus bandwidth Exponential kernel
      Wb[i, j] <- sqrt(exp(-(d[i, j] / h[i]) ^ 2))
    }
  }
  W <- Wb
  best.h <- 'Adaptive Weight'
  
  end.time <- Sys.time()
  cat('Hasil iterasi pembobot :', '\n')
  cat(paste('Best Kernel = ', best[1, 1], ' --- Bandwidth = ', best.h, ' --- GCV = ', best[1, 3]), '\n',
      'Processing time : ', paste(round(end.time - start.time, 2), ' mins'))
  Wlist[[2]] = W
  cat('\n')
  cat('\n')
  
  
  ### Bandwidth Bisquare Kernel
  start.time <- Sys.time()
  h_bi <- rep(0, n)
  GCVmin <- rep(0, n)
  
  for (i in 1:n) {
    A <- 0.0001
    B <- d.max
    iter_bi <- 0
    minGCV <- 0
    selisih <- 1000
    
    while ((selisih > 0.0001) && (iter_bi <= 1000)) {
      h_awal <- seq(A, B, by = (B - A) / 10)
      nh <- length(h_awal)  # banyaknya awalan bandwidth
      GCV <- matrix(0, 1, nh)  # membuat matrix untuk tempat nilai GCV
      colnames(GCV) <- c(h_awal)
      Wb_bi <- matrix(0, n, n)
      
      for (k in 1:nh) {
        h <- h_awal[k]
        for (ii in 1:n) {
          for (jj in 1:n) {
            # Rumus bandwidth Bisquare kernel
            if(d[ii, jj]<h){
              Wb_bi[ii, jj]=((1-(d[ii, jj] / h) ^ 2))^2
            } else {
              Wb_bi[ii, jj]=0
            }
          }
        }
        
        W <- Wb_bi
        x.gcv <- as.matrix(x[-i, ])  # menghilangkan data ke-i dari estimasi
        y1.gcv <- as.matrix(y1[-i, ])  # menghilangkan data y1 ke-i dari estimasi
        y2.gcv <- as.matrix(y2[-i, ])  # menghilangkan data y2 ke-i dari estimasi
        
        beta.gcv <- as.matrix(pop(y1.gcv, y2.gcv, x.gcv, W, param, i, maxit, eps)$param)
        y1hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[1:p]))  # menghitung y1_hat
        y2hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[(p + 1):(2 * p)]))  # menghitung y2_hat
        GCV[k] <- (n * ((abs(y1[i] - y1hat.gcv) + abs(y2[i] - y2hat.gcv)) / (n - length(beta.gcv)) ^ 2))  # menghitung GCV
      }
      hasilGCV_bi <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
      A0 <- A
      B0 <- B
      minGCV <- min(GCV)
      l_min <- which(GCV == minGCV)[1]
      if (l_min == 1) {
        A <- h_awal[l_min]
        B <- h_awal[l_min + 1]
      } else if (l_min == nh) {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min]
      } else {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min + 1]
      }
      selisih <- (B0 - A0) - (B - A)
      iter_bi <- iter_bi + 1
      cat('Bisquare Kernel (kota ke : ', i, ', iterasi ke : ', iter_bi, ', selisih : ', selisih, '\n')
    }
    hasilGCV_bi <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
    hasilGCV <- hasilGCV_bi[order(hasilGCV_bi[, 1]), ]
    GCVmin[i] <- hasilGCV[1, 1]
    h_bi[i] <- hasilGCV[1, 1]
  }
  hasilGCV_bi <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
  bestGCV_bi <- data.frame(h = 0, GCV = sum(GCVmin))  # mengambil bandwidth dengan nilai GCV paling minimum
  
  best <- data.frame(Kernel = c('Bisquare Kernel'), rbind(bestGCV_bi))
  colnames(best) <- c('Kernel', 'Bandwidth', 'GCV')
  best.h <- best[1, 2]
  
  Wb <- matrix(0, n, n)
  h <- h_bi
  for (i in 1:n) {
    for (j in 1:n) {
      # Rumus bandwidth Bisquare kernel
      if(d[i, j]<h[i]){
        Wb[i, j] <- ((1-(d[i, j] / h[i]) ^ 2))^2
      } else {
        Wb[i, j] <- 0
      }
    }
  }
  W <- Wb
  best.h <- 'Adaptive Weight'
  
  end.time <- Sys.time()
  cat('Hasil iterasi pembobot :', '\n')
  cat(paste('Best Kernel = ', best[1, 1], ' --- Bandwidth = ', best.h, ' --- GCV = ', best[1, 3]), '\n',
      'Processing time : ', paste(round(end.time - start.time, 2), ' mins'))
  Wlist[[3]] = W
  
  
  ### Bandwidth Trisquare Kernel
  start.time <- Sys.time()
  h_tri <- rep(0, n)
  GCVmin <- rep(0, n)
  
  for (i in 1:n) {
    A <- 0.0001
    B <- d.max
    iter_tri <- 0
    minGCV <- 0
    selisih <- 1000
    
    while ((selisih > 0.0001) && (iter_tri <= 1000)) {
      h_awal <- seq(A, B, by = (B - A) / 10)
      nh <- length(h_awal)  # banyaknya awalan bandwidth
      GCV <- matrix(0, 1, nh)  # membuat matrix untuk tempat nilai GCV
      colnames(GCV) <- c(h_awal)
      Wb_tri <- matrix(0, n, n)
      
      for (k in 1:nh) {
        h <- h_awal[k]
        for (ii in 1:n) {
          for (jj in 1:n) {
            # Rumus bandwidth Trisquare kernel
            if(d[ii, jj]<h){
              Wb_tri[ii, jj]=((1-(d[ii, jj] / h) ^ 3))^3
            } else {
              Wb_tri[ii, jj]=0
            }
          }
        }
        
        W <- Wb_tri
        x.gcv <- as.matrix(x[-i, ])  # menghilangkan data ke-i dari estimasi
        y1.gcv <- as.matrix(y1[-i, ])  # menghilangkan data y1 ke-i dari estimasi
        y2.gcv <- as.matrix(y2[-i, ])  # menghilangkan data y2 ke-i dari estimasi
        
        beta.gcv <- as.matrix(pop(y1.gcv, y2.gcv, x.gcv, W, param, i, maxit, eps)$param)
        y1hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[1:p]))  # menghitung y1_hat
        y2hat.gcv <- exp(x[i, ] %*% as.matrix(beta.gcv[(p + 1):(2 * p)]))  # menghitung y2_hat
        GCV[k] <- (n * ((abs(y1[i] - y1hat.gcv) + abs(y2[i] - y2hat.gcv)) / (n - length(beta.gcv)) ^ 2))  # menghitung GCV
      }
      hasilGCV_tri <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
      A0 <- A
      B0 <- B
      minGCV <- min(GCV)
      l_min <- which(GCV == minGCV)[1]
      if (l_min == 1) {
        A <- h_awal[l_min]
        B <- h_awal[l_min + 1]
      } else if (l_min == nh) {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min]
      } else {
        A <- h_awal[l_min - 1]
        B <- h_awal[l_min + 1]
      }
      selisih <- (B0 - A0) - (B - A)
      iter_tri <- iter_tri + 1
      cat('Trisquare Kernel (kota ke : ', i, ', iterasi ke : ', iter_tri, ', selisih : ', selisih, '\n')
    }
    hasilGCV_tri <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
    hasilGCV <- hasilGCV_tri[order(hasilGCV_tri[, 1]), ]
    GCVmin[i] <- hasilGCV[1, 1]
    h_tri[i] <- hasilGCV[1, 1]
  }
  hasilGCV_tri <- data.frame(GCV)  # menggabungkan nilai bandwidth awal dengan nilai GCV nya
  bestGCV_tri <- data.frame(h = 0, GCV = sum(GCVmin))  # mengambil bandwidth dengan nilai GCV paling minimum
  
  best <- data.frame(Kernel = c('Trisquare Kernel'), rbind(bestGCV_tri))
  colnames(best) <- c('Kernel', 'Bandwidth', 'GCV')
  best.h <- best[1, 2]
  
  Wb <- matrix(0, n, n)
  h <- h_tri
  for (i in 1:n) {
    for (j in 1:n) {
      # Rumus bandwidth Trisquare kernel
      if(d[i, j]<h[i]){
        Wb[i, j] <- ((1-(d[i, j] / h[i]) ^ 3))^3
      } else {
        Wb[i, j] <- 0
      }
    }
  }
  W <- Wb
  best.h <- 'Adaptive Weight'
  
  end.time <- Sys.time()
  cat('Hasil iterasi pembobot :', '\n')
  cat(paste('Best Kernel = ', best[1, 1], ' --- Bandwidth = ', best.h, ' --- GCV = ', best[1, 3]), '\n',
      'Processing time : ', paste(round(end.time - start.time, 2), ' mins'))
  Wlist[[4]] = W
  
  names(Wlist) = c('gaussian','exponential','bisquare','trisquare')
  return(Wlist)
}
