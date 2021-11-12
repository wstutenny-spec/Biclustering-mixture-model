#' Model selections for \code{bmm}
#'
#' fit several models for bmm along with 3 criteria values: AIC BIC and ICL
#'
#' @param df A data frame that tries to bicluster. Should only contain numeric values.
#' @param range_K All possible number of components. A vector.
#' @param range_Q All possible number of bicluster for each component. A vector
#' @param model The covaraince structure you choose, there are 16 different models belongs to
#'this family:UUU, UUG, UUD, UUC, UGU, UGG, UGD, UGC, GUU, GUG, GUD, GUC, GGU, GGG, GGD, GGC. You can
#'choose more than 1 covarance structure to do model selection.
#' @param permutation Only has effect when model contains UUU, UUG, UUD or UUC. If TRUE,
#'it assume the number of biclusters could be different for different components. If FALSE,
#'it assume the number of biclusters are the same cross all components. Default is FALSE.
#' @param iter Max iterations, defaul is 150.
#'
#'
#'@return A dataframe that contain the cov_str, K, Q, AIC, BIC, ICL values for model. There may
#'be a lot rows if large K and Q, because of lots of combinations: it is a sum of a geometric
#'series with multiplier max(Q) from 1 to max(K).
#'
#'
#'
#'
#'
#'
model_selection <- function(df, range_K, range_Q, model,permutation,iter) {

  res <- NULL


  #registerDoParallel(2)

  ##UGU##
  cov_str <- "UGU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

  ##UGD##
  cov_str <- "UGD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }

  ##UGG##
  cov_str <- "UGG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##UGC##
  cov_str <- "UGC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GUU##
  cov_str <- "GUU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)
      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  #GUD##
  cov_str <- "GUD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GUG##
  cov_str <- "GUG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GUC##
  cov_str <- "GUC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GGU##
  cov_str <- "GGU"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }


  ##GGD##
  cov_str <- "GGD"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GGG##
  cov_str <- "GGG"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }



  ##GGC##
  cov_str <- "GGC"
  if (cov_str %in% model) {
    model_select <- list()
    h <- 1
    for (K in range_K) {
      model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
        Q_K <- rep(j, K)
        paracom(df, K, Q_K, cov_str,iter)

      }
      h <- h + 1
    }
    z <- do.call(rbind, model_select)
    res <- rbind(res, z)
  }




  ##UUU##
  cov_str <- "UUU"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_K) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_K) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }




  ##UUD##
  cov_str <- "UUD"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_K) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_K) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }



  ##UUG##
  cov_str <- "UUG"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_K) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_K) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }


  ##UUC##
  cov_str <- "UUC"
  if (cov_str %in% model) {
    if(permutation==TRUE){
      model_select <- list()
      h <- 1
      for (K in range_K) {
        comb <-
          gtools::permutations(
            n = length(range_Q),
            r = K,
            v = range_Q,
            repeats.allowed = T
          )
        model_select[[h]] <-
          foreach(j = 1:dim(comb)[1], .combine = rbind) %do%  {
            Q_K <- comb[j, ]
            paracom(df, K, Q_K, cov_str,iter)

          }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
    else{
      model_select <- list()
      h <- 1
      for (K in range_K) {
        model_select[[h]] <- foreach(j = range_Q, .combine = rbind) %do%  {
          Q_K <- rep(j, K)
          paracom(df, K, Q_K, cov_str,iter)

        }
        h <- h + 1
      }
      z <- do.call(rbind, model_select)
      res <- rbind(res, z)
    }
  }

  res
}


paracom<-function(df,K,Q_K,cov_str,iter){
  init<-try(Sim_init(df,K,Q_K,cov_str))
  if(is.character(init)|isTRUE(class(init)=="try-error")){res=init}
  else{
    res<-try(bi_fa_T(df,init$mu_K,init$pie_K,init$T_K,init$D_K,init$B_K,cov_str,iter))
  }
  Q_hat=Q_K
  if(isTRUE(class(res)=="try-error")|is.character(res)) {

    f<-data.frame(name=paste(c(cov_str, " K",K," Q",Q_hat),collapse = ""),AIC=-Inf,BIC=-Inf,
                 ICL=-Inf)}
  else {
    f<-data.frame(name=paste(c(cov_str, " K",K," Q",Q_hat),collapse = ""),
                 AIC=res$AIC,BIC=res$BIC,ICL=res$ICL)
  }
  f
}










