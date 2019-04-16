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

  tt <- gimme_eq(rad$vypar, rad[, m], type = "all", message = F)
  t <- tt$equation

  koefs <- rbind(koefs, tt)

  for(i in 1:3){
    rownames(koefs)[k+i] <- paste0(m, "_", odhady[i])
  }

  k = k+3

  kge_t <- data.frame(KGE = best_tr(rad$vypar, rad[, m], type = "all", return = "table", cn = T)$KGE)
  mre_t <- data.frame(MRE = best_tr(rad$vypar, rad[, m], type = "all", return = "table", cn = T)$mre)
  final_eq <- gsub("x", m, t)
  final_eq <- gsub("y", "E", final_eq)

  ult <- data.frame(
    R = rep(NA,3),
    H = rep(NA,3),
    Tw = rep(NA,3),
    Ta = rep(NA,3),
    V = rep(NA,3),
    type = odhady,
    KGE = kge_t,
    MRE = mre_t,
    final_eq)

  ult[, m] <- as.character(t)

  ult_tab <- rbind(ult_tab, ult)

  coupled[, m] <-  na.omit(ult_tab[m])

}

for(i in 1:nrow(ult_tab)){
  f <- ult_tab$final_eq[i]
  all_coef[[f]] <- data.frame(type = ult_tab$type[i], koefs[i, c(1,2)])
}



#======================
#       double
#======================

for(n in sngl){

  pr <- best_tr(rad$vypar, rad[, n], type = "all", return = "data", cn = T)

  for(m in dbl){

    if(m != n){

      t<-c()
      kge_t <- c()
      mre_t <- c()

      k = nrow(koefs)

      for(i in 1:3){
        t <- c(t, as.character(gimme_eq(rad$vypar-pr[, i], rad[, m], type = "lm", message = F)$equation))
        tt <- gimme_eq(rad$vypar-pr[, i], rad[, m], type = "lm", message = F)
        koefs <- rbind(koefs, tt)
        rownames(koefs)[k+i] <- paste0(n, "_", colnames(pr)[i], "_", m)
        h <-  pr[, i] + best_tr(rad$vypar-pr[, i], rad[, m], type = "lm", return = "data")
        kge_t <- c(kge_t, round(KGE(as.numeric(h), rad$vypar), 4))
        mre_t <- c(mre_t, round(mean(abs(as.numeric(h) - rad$vypar)/rad$vypar)*100, 2))
      }

      k = which(rownames(koefs) == paste0(n, "_lm"))
      lk = which(rownames(koefs) == paste0(n, "_lm_", m))

      final_eq <- c()

      for(i in 0:2){

        if(i == 0){

          if(koefs$a[lk+i]>0){
            feq <- paste0("E = ", round(koefs$a[k+i], 4), n, " + ", round(koefs$a[lk+i], 4), m)
          }else{
            feq <- paste0("E = ", round(koefs$a[k+i], 4), n, " - ", abs(round(koefs$a[lk+i], 4)), m)
          }

        }else{

          tt <- gsub("y", "E", koefs$equation[k+i])
          tt <- gsub("x", n, tt)

          if(koefs$a[lk+i]>0){
            feq <- paste0(tt," + ", round(koefs$a[lk+i], 4), m)
          }else{
            feq <- paste0(tt, " - ", abs(round(koefs$a[lk+i], 4)), m)
          }
        }

        s <- extractSUM(koefs$equation[k+i], koefs$equation[lk+i])

        if(s > 0){
          feq <- paste0(feq," + ", round(abs(s), 4))
        }else{
          feq <- paste0(feq," - ", round(abs(s), 4))
        }


        final_eq <- c(final_eq, feq)

      }

      ult <- data.frame(
        R = rep(NA,3),
        H = rep(NA,3),
        Tw = rep(NA,3),
        Ta = rep(NA,3),
        V = rep(NA,3),
        type = c("lm", "power", "exp"),
        KGE = kge_t,
        MRE = mre_t,
        final_eq)

      ult[, n] <- coupled[, n]
      ult[, m] <- t

      ult_tab <- rbind(ult_tab, ult)
    }
  }
}

end <- length(all_coef)+1

for(i in end:nrow(ult_tab)){
  f <- ult_tab$final_eq[i]

  eq_var <- extractVAR(f)

  all_coef[[f]] <- data.frame(type = ult_tab$type[i],
                              koefs[which(rownames(koefs) == paste0(eq_var[1], "_", ult_tab$type[i])), c(1,2)],
                              koefs[i, c(1,2)])
}

#======================
#       triple
#======================

dbl <- "H"
trpl <- c("Tw", "V", "Ta")

for(n in sngl){

  pr <- best_tr(rad$vypar, rad[, n], type = "all", return = "data", cn = T)

  for(m in dbl){

    for(p in trpl){

      if(m != n & m != p & n != p){

        t <- c()
        tt <- c()
        kge_t <- c()
        mre_t <- c()

        k = nrow(koefs)

        for(i in 1:3){

          t <- c(t, as.character(gimme_eq(rad$vypar-pr[, i], rad[, m], type = "lm", message = F)$equation))

          h <- pr[, i] + best_tr(rad$vypar-pr[, i], rad[, m], type = "lm", return = "data")
          hw <- h + best_tr(rad$vypar - h, rad[, p], type = "lm", return = "data")
          ttt <- gimme_eq(rad$vypar - h, rad[, p], type = "lm", message = F)

          tt <- c(tt, as.character(ttt$equation))

          koefs <- rbind(koefs, ttt)
          rownames(koefs)[k+i] <- paste0(n, "_", colnames(pr)[i], "_", m, "_", p)
          kge_t <- c(kge_t, round(KGE(as.numeric(hw), rad$vypar), 4))
          mre_t <- c(mre_t, round(mean(abs(as.numeric(hw) - rad$vypar)/rad$vypar)*100, 2))
        }

        k = which(rownames(koefs) == paste0(n, "_lm"))
        k2 = which(rownames(koefs) == paste0(n, "_lm_", m))
        lk = which(rownames(koefs) == paste0(n, "_lm_", m, "_", p))

        final_eq <- c()

        for(i in 0:2){

          if(i == 0){

            feq <- paste0("E = ", round(koefs$a[k+i], 4), n)

            if(koefs$a[k2+i]>0){
              feq <- paste0(feq, " + ", round(koefs$a[k2+i], 4), m)
            }else{
              feq <- paste0(feq, " - ", abs(round(koefs$a[k2+i], 4)), m)
            }

            if(koefs$a[lk+i]>0){
              feq <- paste0(feq, " + ", round(koefs$a[lk+i], 4), p)
            }else{
              feq <- paste0(feq, " - ", abs(round(koefs$a[lk+i], 4)), p)
            }


          }else{

            feq <- gsub("x", n, gsub("y", "E", koefs$equation[k+i]))

            if(koefs$a[k2+i]>0){
              feq <- paste0(feq, " + ", round(koefs$a[k2+i], 4), m)
            }else{
              feq <- paste0(feq, " - ", abs(round(koefs$a[k2+i], 4)), m)
            }

            if(koefs$a[lk+i]>0){
              feq <- paste0(feq, " + ", round(koefs$a[lk+i], 4), p)
            }else{
              feq <- paste0(feq, " - ", abs(round(koefs$a[lk+i], 4)), p)
            }
          }

          s <- extractSUM(koefs$equation[k+i], koefs$equation[k2+i])
          s <- s + numextract(koefs$equation[lk+i], "n")

          if(s > 0){
            feq <- paste0(feq," + ", round(abs(s), 4))
          }else{
            feq <- paste0(feq," - ", round(abs(s), 4))
          }


          final_eq <- c(final_eq, feq)

        }


        ult <- data.frame(
          R = rep(NA,3),
          H = rep(NA,3),
          Tw = rep(NA,3),
          Ta = rep(NA,3),
          V = rep(NA,3),
          type = c("lm", "power", "exp"),
          KGE = kge_t,
          MRE = mre_t,
          final_eq)

        ult[, n] <- coupled[, n]
        ult[, m] <- t
        ult[, p] <- tt

        ult_tab <- rbind(ult_tab, ult)

      }
    }
  }
}

end <- length(all_coef)+1

for(i in end:nrow(ult_tab)){
  f <- ult_tab$final_eq[i]

  eq_var <- extractVAR(f)

  all_coef[[f]] <- data.frame(type = ult_tab$type[i],
                              koefs[which(rownames(koefs) == paste0(eq_var[1], "_", ult_tab$type[i])), c(1,2)],
                              koefs[which(rownames(koefs) == paste0(eq_var[1], "_", ult_tab$type[i], "_", eq_var[2])), c(1,2)],
                              koefs[i, c(1,2)])
}

return(list(all_coef = all_coef, ult_tab = ult_tab, koefs = koefs))

}


