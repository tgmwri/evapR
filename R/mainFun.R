#' scatter_fun
#'
#' This function plots trend.
#' @param x regressors
#' @param y dependent variable
#' @param ysim fitted values
#' @param title plot title
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' scatter_fun()

scatter_fun <- function(y, x, ysim, title = "Titulek"){
  ggplot() +
    geom_point(aes(x = x, y = y)) +
    geom_line(aes(y = ysim, x = x), col = "red", size = 0.5)+
    labs(title = title) +
    theme_bw()
}



#' Best Model Function
#'
#' This function estimates coefficients for three different types of model.
#' @param data data to compute trends
#' @param x regressor
#' @param y dependent variable
#' @param type type of model to use: linear - "lm", exponential - "exp" and power trend "power". Default "all" use all three types.
#' @param return what kind of data to return:
#' "data" - fitted values;
#' "plot" - ggplot of fitted values;
#' "table" - contains information about coefficint of determination (R2), standard deviation (sd), mean absolute error (mae),
#' mean relative error (mre) and Kling-Gupta efficiency (KGE);
#' "all" - returns all;
#' @param cn whether to keep column names set by function (default is FALSE)
#'
#' @import hydroGOF purrr cowplot
#' @export
#' @examples
#' est <- best_tr(dta$evaporation, dta$radiation, cn=TRUE)
#' df <- as.data.frame(est$dta.n)
#' est$tab
#' plot_grid(plotlist = est$gg.all)
#'
best_tr <- function(y, x, type = "all", return = "all", cn = FALSE){

  dta.n <- c()
  tab <- c()


  tit <- data.frame(type = c("lm", "power", "exp"), main = c("Lineární spojnice trendu",
                                                             "Mocninná spojnice trendu",
                                                             "Exponenciální spojnice trendu"))
  if("all" %in% type){do.types <- c("lm", "power", "exp")}else{do.types = type}

  for(i in 1:length(do.types)){

    if(do.types[i] == "lm"){

      model <- lm(y~x)
      dta.n <- cbind(dta.n, coef(model)[2]*x + coef(model)[1])
      colnames(dta.n)[i] <- do.types[i]

    }else if(do.types[i] == "power"){

      model <- lm(log(y)~log(x))
      dta.n <- cbind(dta.n, exp(coef(model)[1])*(x^coef(model)[2]))
      colnames(dta.n)[i] <- do.types[i]

    }else if(do.types[i] == "exp"){

      model <- lm(log(y)~x)
      dta.n <- cbind(dta.n, exp(coef(model)[1])*exp(coef(model)[2]*x))
      colnames(dta.n)[i] <- do.types[i]

    }

    if(return %in% c("all", "table")){

      p <- data.frame(
        r2 = round(summary(model)$r.squared, 4),
        sd = round(sd(y - dta.n[, do.types[i]]), 4),
        mae = round(mae(dta.n[, do.types[i]], y), 4),
        mre = round(mean(abs(dta.n[, do.types[i]] - y)/y)*100, 2),
        KGE = round(KGE(dta.n[, do.types[i]], y), 4),
        row.names = do.types[i]
      )

      tab <- rbind(tab, p)

    }

    if(return %in% c("all", "plot")){
      gg.all <- map(do.types, ~scatter_fun(y, x, ysim = dta.n[, .x], title = tit$main[tit$type == .x]))

    }

  }

  if(return == "all"){
    if(type == "all"){
      message(paste0("Doporučuju: ", rownames(tab[tab$r2 == max(tab$r2),]), ", R2 = ", max(tab$r2)))
    }
    if(cn == F){colnames(dta.n) <- NULL}
    return(list(dta.n = dta.n, tab = tab[-1], gg.all = gg.all))

  }else if(return == "table"){
    if(type == "all"){
      message(paste0("Doporučuju: ", rownames(tab[tab$r2 == max(tab$r2),]), ", R2 = ", max(tab$r2)))
      }
    return(tab = tab[-1])

  }else if(return == "plot"){
    return(gg.all = gg.all)

  }else if(return == "data"){
    if(cn == F){colnames(dta.n) <- NULL}
    return(dta.n)
  }

}


#' Return model equation
#'
#'
#' Function that returns model equation and it's coefficients.
#' @param x regressor
#' @param y dependent variable
#' @param type type of trend to use: linear - "lm", exponential - "exp" and power trend "power" or "all". Default is "lm".
#' @param eqOnly default is FALSE, if TRUE, returns only equation without coefficients.
#' @param message whether to print equation to the console (default is TRUE).
#'
#' @export
#' @examples
#' gimme_eq(est$evaporation, eat$radiation, type = "lm", eqOnly = T)

gimme_eq <- function(y, x, type = "lm", eqOnly = FALSE, message = TRUE){

  a_f <- c()
  b_f <- c()
  eq_f <- c()

  if(type %in% c("lm", "all")){

    model <- lm(y~x)
    a <- as.numeric(coef(model)[2])
    b <- as.numeric(coef(model)[1])

    if(round(b,4) < 0){
      eq <- paste0("y = ", round(a,4), "x - ", abs(round(b,4)))
    }else{
      eq <- paste0("y = ", round(a,4), "x + ", round(b,4))
    }

    if(message == TRUE){
      message(eq)
    }

    a_f <- c(a_f, a)
    b_f <- c(b_f, b)
    eq_f <- c(eq_f, eq)

  }

  if(type %in% c("power", "all")){

    model <- lm(log(y)~log(x))
    a <- as.numeric(exp(coef(model)[1]))
    b <- as.numeric(coef(model)[2])

    eq <- paste0("y = ", round(a,4), "x^", round(b,4))

    if(message == TRUE){
      message(eq)
    }

    a_f <- c(a_f, a)
    b_f <- c(b_f, b)
    eq_f <- c(eq_f, eq)

  }

  if(type %in% c("exp", "all")){

    model <- lm(log(y)~x)
    a <- as.numeric(exp(coef(model)[1]))
    b <- as.numeric(coef(model)[2])

    eq <- paste0("y = ", round(a,4), "e^(", round(b,4), "x)")

    if(message == TRUE){
      message(eq)
    }

    a_f <- c(a_f, a)
    b_f <- c(b_f, b)
    eq_f <- c(eq_f, eq)

  }

  if(eqOnly == FALSE){
    return(data.frame(a = a_f, b = b_f, equation = eq_f))
  }

}


#' numextract
#'
#' Function to extract constant term from the equation string for further manipulations
#'
#' @param string equation
#' @param ret whether returned value is going to be string ("s", default) or numeric ("n")
#'
#' @import stringr
#' @export
#' @examples
#' numextract("E = 0.4254*R + 1.0341", "n")
#' 1.0341

numextract <- function(string, ret = "s"){
  if(ret == "n"){
    s <- str_extract(string, "[-]? \\d*\\.\\d*$")
    s <- gsub(" ", "", s)
    s <- as.numeric(s)
    return(s)
  }else if(ret == "s"){
    s <- str_extract(string, "[-]? \\d*\\.\\d*$")
    return(s)
  }
}

#' extractSUM
#'
#' Function to calculate sum of extracted constant terms from equations
#'
#' @param x first equation
#' @param y second equation
#'
#' @export
#' @examples
#' extractSUM("E = 0.4254R + 1.0341", "E = 7.5687Tvzd - 82.29")
#' -81.2559

extractSUM <- function(x, y){

  v1 <- numextract(x, "n")

  for(i in 1:length(v1)){
    if(is.na(v1[i]) == T){
      v1[i] = 0
    }
  }

  v2  <- numextract(y, "n")

  for(i in 1:length(v2)){
    if(is.na(v2[i]) == T){
      v2[i] = 0
    }
  }

  num <- cbind(v1, v2)
  num <- rowSums(num)

  return(num)
}

#' extractVAR
#'
#' Function to extract variables from equations.
#' @param e equation
#' @export
#' @examples
#' extractVAR("E = 0.0156R^1.0156 - 0.8136H + 0.0106Tw + 0.4647")
#' "R"  "H"  "Tw"

extractVAR <- function(e){
  ep <- gsub("[-]? \\d.\\d*", " ", e)
  ep <- gsub("E = ", "", ep)
  ep <- gsub("e\\^", "", ep)
  ep <- gsub("\\^", "", ep)
  ep <- gsub("\\d.\\d", "", ep)
  ep <- gsub("\\+|\\-", "", ep)
  ep <- gsub("\\(|\\)", "", ep)
  ep <- unlist(strsplit(ep, " "))
  ep <- ep[ep != ""]
  ep <- ep[ep != "0"]

  return(ep)
}
