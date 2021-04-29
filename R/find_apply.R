#' find_fun
#'
#' Function to search equations with selected variables
#'
#' @param IN string or vector of variables
#' @param E list of equations to search in
#'
#' @import utils
#' @export
#' @examples
#'
#' eq_t <- find_fun(c("Ta", "H", "R"), eq_list)
#' "E = 0.0169R - 1.0099H + 0.022Ta + 0.5102" "E = 0.0169R + 0.0374Ta - 0.4965"          "E = 0.0156R^1.0156 + 0.0336Ta - 0.4588"
#' "E = 0.0156R^1.0156 - 0.8136H + 0.6596"    "E = 0.0169R - 1.0099H + 0.829"



find_fun <- function(IN, E){
  
  # E <- readRDS(file.path(.dir,"data", "vzorecky_git.rds"))
  # E <- c("0.4254*R^1.0341", 
  #        "0.5136*R-32.894*H+26.5149",
  #        "0.5136*R+0.9423*Tw+8.035*H-22.5585",
  #        "0.5136*R+0.9423*Tw+1.5702*V-18.2897",
  #        "0.5136*R+0.9423*Tw-16.1031",
  #        "0.5136*R+1.2211*Ta-16.6161",  
  #        "1.1299*Tw^1.5123",
  #        "1.1299*Tw^1.5123-83.368*H+67.862",
  #        "1.1299*Tw^1.5123+17.52*V-95.408*H+53.074",
  #        "6.8712*Tw+15.175*V-52.491",
  #        "1.1299*Tw^1.5123+17.52*V-23.126",
  #        "7.5687*Ta-82.298*H+49.068",
  #        "2.7521*Ta^1.3064",
  #        "7.5687*Ta+16.648*V-93.739*H+35.015",
  #        "7.5687*Ta+16.648*V-39.852"
  # )
  
  var <- c("R", "H", "*.Tw.*", "*.Ta.*", "V")
  cur <- sapply(var, grepl, E, ignore.case = T)
  
  colnames(cur) <- c("R", "H", "Tw", "Ta", "V")
  
  eq_f <- c()
  
  
  for(i in 1:length(IN)){
    
    rnd <- which(rowSums(cur[,]) == 1)
    
    if(length(rnd) > 0){
      
      for(l in 1: length(rnd)){
        if(cur[rnd[l], IN[i]] == T){
          eq_f <- c(eq_f, E[rnd[l]])
        }
      }
      
    }
    
  }
  
  
  tmp_n <- c()  
  for(i in 1:nrow(cur)){
    tmp <- ifelse(all(cur[i, colnames(cur) %in% IN] == T) & all(cur[i, !(colnames(cur) %in% IN)] == F), T, F)
    tmp_n <- rbind(tmp_n, tmp)
  }
  
  if(length(E[tmp_n]) != 0){
    ifelse(!(E[tmp_n] %in% eq_f), eq_f <- c(eq_f, E[tmp_n]), eq_f <- eq_f)
  }
  
  if(length(IN) >= 2){
    inc <- combn(IN, 2)
    
    
    tmp_f <- c()
    
    for(j in 1:ncol(inc)){
      tmp_n <- c() 
      for(i in 1:nrow(cur)){
        tmp <- ifelse(all(cur[i, colnames(cur) %in% inc[, j]] == T) & 
                        all(cur[i, !(colnames(cur) %in% inc[, j])] == F), T, F)
        tmp_n <- rbind(tmp_n, tmp)
      }
      
      tmp_f <- cbind(tmp_f, tmp_n)
    }
    
    tmp_f <- as.data.frame(tmp_f)
    
    for(i in 1:ncol(inc)){
      
      colnames(tmp_f)[i] <- paste(inc[, i], collapse = "_")
      row.names(tmp_f) <- NULL
      
      k <- tmp_f[,i]
      
      if(length(E[k]) != 0){
        for(j in 1:length(E[k])){
          
          if(!(E[k][j] %in% eq_f)){
            eq_f <- c(eq_f, E[k])
          }
          
        }
        
      }
    }
    
  }
  
  return(eq_f)
  
}


#' apply_fun
#'
#' Function to make predictions based on inputed equations
#'
#' @param data data on which predictions are going to be calculated
#' @param coeff list of equation's coefficients
#' @param eq_list list of equations to apply
#'
#' @export
#' @examples
#'
#' eq_t <- find_fun(c("Ta", "H", "R"), eq_list)

apply_fun <- function(data, coeff, eq_list){

  col <- c("type", "a1", "b1", "a2", "b2", "a3", "b3")
  vysledek <- c()

  for(i in 1:length(eq_list)){

    E <- as.data.frame(coeff[eq_list[i]])
    n <- length(E)
    colnames(E) <- col[1:n]

    variables <- extractVAR(eq_list[i])

    if(E$type == "lm"){

      V <- c(0)
      k = 0

      for(j in 1:length(variables)){
        k=k+2
        v <- E[,k] * data[, variables[j]] + E[,k+1]
        V = V + v
      }

      vysledek <- cbind(vysledek, V)
      colnames(vysledek)[i] <- paste(paste(variables, collapse = '_'), E$type, sep = '_')
    }


    if(E$type == "exp"){

      V <- c(0)
      k = 0

      for(j in 1:length(variables)){
        k=k+2
        if(k == 2){
          v <- E[,k] * exp(E[,k+1] * data[, variables[j]])
        }else{
          v <- E[,k] * data[, variables[j]] + E[,k+1]
        }
        V = V + v
      }

      vysledek <- cbind(vysledek, V)
      colnames(vysledek)[i] <- paste(paste(variables, collapse = '_'), E$type, sep = '_')
    }

    if(E$type == "power"){

      V <- c(0)
      k = 0

      for(j in 1:length(variables)){
        k=k+2
        if(k == 2){
          v <- E[,k] * data[, variables[j]]^E[,k+1]
        }else{
          v <- E[,k] * data[, variables[j]] + E[,k+1]
        }
        V = V + v
      }

      vysledek <- cbind(vysledek, V)
      colnames(vysledek)[i] <- paste(paste(variables, collapse = '_'), E$type, sep = '_')
    }

  }

  return(vysledek = as.data.frame(vysledek))

}




