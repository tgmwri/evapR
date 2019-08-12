#' compute_coef_table
#'
#' Function to compute possible combinations of equations and their coefficients, returns parameters in three forms:
#' list of equations with every coefficient (all_coef);
#' table (ult_tab), where columns are variables and rows are equations, contains information about KGE and MRE with final equation for every row (for each combination);
#' table (koefs), which contains coefficients and equation for each combination.
#'
#' @param rad input data.frame
#'
#' @import utils
#' @export
#' @examples
#' recomended
#' options(stringsAsFactors = FALSE)

compute_coef_table <- function(rad){
  
  odhady <- c("lm", "power", "exp")
  
  sngl <- c("R", "Tw", "Ta")
  dbl <- c("H", "Tw", "Ta", "V")
  
  all_coef <- list()
  
  #======================
  #       single
  #======================
  
  koefs <- c()
  ult_tab <- data.frame()
  coupled <- data.frame(R = rep(NA,3),
                        Tw = rep(NA,3),
                        Ta = rep(NA,3))
  
  k = 0
  
  for(m in sngl){
    
    if(!all(is.na(rad[,m]))){
      
      tt <- gimme_eq(rad$vypar, rad[, m], type = "all", message = F)
      tt <- tt[-4,]
      t <- tt$equation
      
      koefs <- rbind(koefs, tt)
      
      for(i in 1:3){
        rownames(koefs)[k+i] <- paste0(m, "_", odhady[i])
      }
      
      k = k+3
      
      kge_t <- data.frame(KGE = best_tr(rad$vypar, rad[, m], type = c("lm", "power", "exp"), return = "table", cn = T)$KGE)
      mre_t <- data.frame(MRE = best_tr(rad$vypar, rad[, m], type = c("lm", "power", "exp"), return = "table", cn = T)$mre)
      final_eq <- gsub("x", m, t)
      final_eq <- gsub("y", "E", final_eq)
      
      ult <- data.frame(
        R = rep(NA,3), 
        H = rep(NA,3), 
        Tw = rep(NA,3), 
        Ta = rep(NA,3),
        V = rep(NA,3), 
        type = odhady[-4],
        KGE = kge_t, 
        MRE = mre_t, 
        final_eq)
      
      ult[, m] <- as.character(t)
      
      ult_tab <- rbind(ult_tab, ult)
      
      coupled[, m] <-  na.omit(ult_tab[m])
      
      
      for(i in 1:nrow(ult_tab)){
        f <- ult_tab$final_eq[i]
        all_coef[[f]] <- data.frame(type = ult_tab$type[i], koefs[i, c(1,2)])
      }
      
    }
    
  }
  
  #======================
  #       double
  #======================
  
  coupled2 <- coupled[rep(1:nrow(coupled), each=2), ]
  
  for(n in sngl){
    
    if(!all(is.na(rad[,n]))){
      
      pr <- best_tr(rad$vypar, rad[, n], type = c("lm", "power", "exp"), return = "data", cn = T)
      
      for(m in dbl){
        
        if(m != n & !all(is.na(rad[,m]))){
          
          t <- c()
          lnt <- c()
          kge_t <- c()
          mre_t <- c()
          
          k = nrow(koefs)
          
          for(i in 1:3){
            
            #----linearni---------------
            tt <- gimme_eq(rad$vypar-pr[, i], rad[, m], type = "lm", message = F)
            t <- c(t, tt$equation)
            koefs <- rbind(koefs, tt)
            rownames(koefs)[k+i] <- paste0(n, "_", colnames(pr)[i], "_", m, "_lm")
            h <-  pr[, i] + best_tr(rad$vypar-pr[, i], rad[, m], type = "lm", return = "data")
            kge_t <- c(kge_t, round(KGE(as.numeric(h), rad$vypar), 4))
            mre_t <- c(mre_t, round(mean(abs(as.numeric(h) - rad$vypar)/rad$vypar, na.rm = T)*100, 2))
            
            #----logaritmicka-----------
            lntt <- gimme_eq(rad$vypar-pr[, i], rad[, m], type = "log", message = F)
            lnt <- c(lnt, lntt$equation)
            koefs <- rbind(koefs, lntt)
            rownames(koefs)[k+(i+1)] <- paste0(n, "_", colnames(pr)[i], "_", m, "_log")
            h_ln <-  pr[, i] + best_tr(rad$vypar-pr[, i], rad[, m], type = "log", return = "data")
            kge_t <- c(kge_t, round(KGE(as.numeric(h_ln), rad$vypar), 4))
            mre_t <- c(mre_t, round(mean(abs(as.numeric(h_ln) - rad$vypar)/rad$vypar, na.rm = T)*100, 2))
            
            k = k+1
          }
          
          
          k = which(rownames(koefs) == paste0(n, "_lm"))
          lk = which(rownames(koefs) == paste0(n, "_lm_", m, "_lm"))
          lk_ln = which(rownames(koefs) == paste0(n, "_lm_", m, "_log"))
          
          final_eq <- c()
          j = 0 
          
          for(i in c(0,2,4)){
            
            
            if(i == 0){
              
              if(koefs$a[lk+i]>0){
                feq <- paste0("E = ", round(koefs$a[k+j], 4), n, " + ", round(koefs$a[lk+i], 4), m)
              }else{
                feq <- paste0("E = ", round(koefs$a[k+j], 4), n, " - ", abs(round(koefs$a[lk+i], 4)), m)
              }
              
              if(koefs$a[lk_ln+i]>0){
                feq2 <- paste0("E = ", round(koefs$a[k+j], 4), n, " + ", round(koefs$a[lk_ln+i], 4), " ln(", m, ")")
              }else{
                feq2 <- paste0("E = ", round(koefs$a[k+j], 4), n, " - ", abs(round(koefs$a[lk_ln+i], 4)), " ln(", m, ")")
              }
              
              
            }else{
              
              tt <- gsub("y", "E", koefs$equation[k+j])
              tt <- gsub("x", n, tt)
              
              if(koefs$a[lk+i]>0){
                feq <- paste0(tt," + ", round(koefs$a[lk+i], 4), m)
              }else{
                feq <- paste0(tt, " - ", abs(round(koefs$a[lk+i], 4)), m)
              }
              
              
              lntt <- gsub("y", "E", koefs$equation[k+j])
              lntt <- gsub("x", n, lntt)
              
              if(koefs$a[lk_ln+i]>0){
                feq2 <- paste0(lntt," + ", round(koefs$a[lk_ln+i], 4), " ln(", m, ")")
              }else{
                feq2 <- paste0(lntt, " - ", abs(round(koefs$a[lk_ln+i], 4)), " ln(", m, ")")
              }
              
              
            }
            
            s <- extractSUM(koefs$equation[k+i], koefs$equation[lk+i])
            s2 <- extractSUM(koefs$equation[k+i], koefs$equation[lk_ln+i])
            
            if(s > 0){
              feq <- paste0(feq," + ", round(abs(s), 4))
            }else{
              feq <- paste0(feq," - ", round(abs(s), 4))
            }
            
            if(s2 > 0){
              feq2 <- paste0(feq2," + ", round(abs(s2), 4))
            }else{
              feq2 <- paste0(feq2," - ", round(abs(s2), 4))
            }
            
            
            final_eq <- c(final_eq, feq, feq2)
            
            j = j+1
            
          }
          
          
          ult <- data.frame(
            R = rep(NA,6), 
            H = rep(NA,6), 
            Tw = rep(NA,6), 
            Ta = rep(NA,6),
            V = rep(NA,6), 
            type = c("lm_lm", "lm_log", "power_lm", "power_log", "exp_lm", "exp_log"),
            KGE = kge_t, 
            MRE = mre_t,
            final_eq)
          
          
          
          ult[, n] <- coupled2[,n]
          
          xlt <- c()
          for(i in 1:3){
            xlt <- c(xlt, t[i], lnt[i])
          }
          
          ult[, m] <- xlt
          
          ult_tab <- rbind(ult_tab, ult)
        }
      }
      
      end <- length(all_coef)+1
      
      for(i in end:nrow(ult_tab)){
        f <- ult_tab$final_eq[i]
        
        eq_var <- extractVAR(f)
        
        all_coef[[f]] <- data.frame(type = ult_tab$type[i], 
                                    koefs[which(rownames(koefs) == paste0(eq_var[1], "_", gsub("_log|_lm", "", ult_tab$type[i]))), c(1,2)],
                                    koefs[i, c(1,2)])
      }
    }
  }
  
  #======================
  #       triple
  #======================
  
  lk = length(koefs)+1
  trpl <- c("Tw", "V", "Ta")  
  
  for(p in trpl){
    
    for(j in 10:length(all_coef)){ 
      
      t <- tt <- c()
      lnt <- lntt <- c()
      kge_t <- mre_t <- c()
      final_eq <- c()
      
      ww <- all_coef[[j]]
      w <- extractVAR(names(all_coef[j]))
      
      if(w[2] != w[1] & w[2] != p & w[1] != p){
        
        if(!all(is.na(rad[, w[1]])) & !all(is.na(rad[, w[2]])) & !all(is.na(rad[, p]))){
          
          if(ww$type == 'lm_lm'){
            h <-  ww$a*rad[, w[1]] + ww$a.1*rad[, w[2]] + ww$b + ww$b.1
          }else if(ww$type == 'lm_log'){
            h <-  ww$a*rad[, w[1]] + ww$a.1*log(rad[, w[2]]) + ww$b + ww$b.1
          }else if(ww$type == 'power_lm'){
            h <-  ww$a*rad[, w[1]] ^ ww$b + ww$a.1*rad[, w[2]] + ww$b.1
          }else if(ww$type == 'power_log'){
            h <-  ww$a*rad[, w[1]] ^ ww$b + ww$a.1*log(rad[, w[2]]) + ww$b.1
          }else if(ww$type == 'exp_lm'){
            h <-  ww$a*exp(rad[, w[1]] * ww$b) + ww$a.1*rad[, w[2]] + ww$b.1
          }else if(ww$type == 'exp_log'){
            h <-  ww$a*exp(rad[, w[1]] * ww$b) + ww$a.1*log(rad[, w[2]]) + ww$b.1
          }
          
          pr <- best_tr(rad$vypar, rad[, w[1]], type = c("lm", "power", "exp"), return = "data", cn = T)
          k = nrow(koefs)   
          tp = gsub('_log|_lm', '', ww$type)
          
          #----linearni---------------
          
          t <- c(t, as.character(gimme_eq(rad$vypar-pr[, tp], rad[, w[2]], type = "lm", message = F)$equation))
          
          hw <- h + best_tr(rad$vypar - h, rad[, p], type = "lm", return = "data")
          ttt <- gimme_eq(rad$vypar - h, rad[, p], type = "lm", message = F)
          
          tt <- c(tt, as.character(ttt$equation))
          
          koefs <- rbind(koefs, ttt)
          
          st <- unlist(strsplit(ww$type, '_'))
          rownames(koefs)[k+1] <- paste0(w[1], "_", st[1], "_", w[2], "_" , st[2], "_", p, "_lm")
          kge_t <- c(kge_t, round(KGE(as.numeric(hw), rad$vypar), 4))
          mre_t <- c(mre_t, round(mean(abs(as.numeric(hw) - rad$vypar)/rad$vypar, na.rm = T)*100, 2))
          
          #----logaritmicka-----------  
          
          lnt <- c(lnt, as.character(gimme_eq(rad$vypar-pr[, tp], rad[, w[2]], type = "log", message = F)$equation))
          
          hw_ln <- h + best_tr(rad$vypar - h_ln, rad[, p], type = "log", return = "data")
          lnttt <- gimme_eq(rad$vypar - h, rad[, p], type = "log", message = F)
          
          lntt <- c(lntt, as.character(lnttt$equation))
          
          koefs <- rbind(koefs, lnttt)
          
          st <- unlist(strsplit(ww$type, '_'))
          rownames(koefs)[k+2] <- paste0(w[1], "_", st[1], "_", w[2], "_" , st[2], "_", p, "_log")
          kge_t <- c(kge_t, round(KGE(as.numeric(hw_ln), rad$vypar), 4))
          mre_t <- c(mre_t, round(mean(abs(as.numeric(hw_ln) - rad$vypar)/rad$vypar, na.rm = T)*100, 2))
          
          
          
          f <- names(all_coef[j])
          feq <- gsub('.{8}$', '', f)
          
          if(koefs$a[lk]>0){
            feq <- paste0(feq, " + ", round(koefs$a[lk], 4), p)
          }else{
            feq <- paste0(feq, " - ", abs(round(koefs$a[lk], 4)), p)
          }
          
          s <- extractSUM(f, koefs$equation[lk])
          
          if(s > 0){
            feq <- paste0(feq," + ", round(s, 4))
          }else{
            feq <- paste0(feq," - ", round(abs(s), 4))
          }
          
          
          feq2 <- gsub('.{8}$', '', f)
          
          if(koefs$a[lk+1]>0){
            feq2 <- paste0(feq2, " + ", round(koefs$a[lk+1], 4), " ln(", p, ")")
          }else{
            feq2 <- paste0(feq2, " - ", abs(round(koefs$a[lk+1], 4)), " ln(", p, ")")
          }
          
          s2 <- extractSUM(f, koefs$equation[lk+1])
          
          if(s2 > 0){
            feq2 <- paste0(feq2," + ", round(s2, 4))
          }else{
            feq2 <- paste0(feq2," - ", round(abs(s2), 4))
          }
          
          final_eq <- c(final_eq, feq, feq2)
          lk <- lk+2
          
          
          ult <- data.frame(
            R = rep(NA,2), 
            H = rep(NA,2), 
            Tw = rep(NA,2), 
            Ta = rep(NA,2),
            V = rep(NA,2), 
            type = c(paste0(ww$type, '_lm'), paste0(ww$type, '_log')),
            KGE = kge_t, 
            MRE = mre_t,
            final_eq)
          
          if(gsub('_lm|_log', '', ww$type) == 'lm'){
            x <- 1:2
          }else if(gsub('_lm|_log', '', ww$type) == 'power'){
            x <- 3:4
          }else{
            x <- 5:6
          }
          
          if(gsub('lm_|power_|exp_', '', ww$type) == 'lm'){
            xx <- rep(t, 2)
          }else{
            xx <- rep(lnt, 2)
          }
          
          ult[, w[1]] <- coupled2[x,w[1]]
          ult[, w[2]] <- xx
          ult[, p] <- c(tt, lntt)
          
          ult_tab <- rbind(ult_tab, ult)
        }
      }
    }
  }
  
  
  end <- length(all_coef)+1
  
  if(end < nrow(ult_tab)){
    
    for(i in end:nrow(ult_tab)){
      f <- ult_tab$final_eq[i]
      
      eq_var <- extractVAR(f)
      st <-unlist(strsplit(ult_tab$type[i], '_'))
      
      all_coef[[f]] <- data.frame(type = ult_tab$type[i], 
                                  koefs[which(rownames(koefs) == paste(eq_var[1], st[1], sep = "_", collapse = "_")), c(1,2)],
                                  koefs[which(rownames(koefs) == paste(eq_var[1:2], st[1:2], sep = "_", collapse = "_")), c(1,2)],
                                  koefs[i, c(1,2)])
    }
  }
  
  
  return(list(all_coef = all_coef, ult_tab = ult_tab, koefs = koefs))
  
}
