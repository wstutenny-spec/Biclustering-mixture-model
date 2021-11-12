#' Biclustering mixture model algorithm
#'
#'
#' Main function that can do biclustering mixture model and select the best model based on BIC, AIC or ICL.
#' @param df A data frame that tries to bicluster. Should only contain numeric values.
#'
#' @param range_K When \code{single=FALSE}, the input is all possible number of components: a vector.
#' When \code{single=TRUE}, one number specify the number of components.
#' @param range_Q When \code{single=FLASE}, the input is all possible number of bicluster for each component.
#' When \code{single=TRUE}, a vector that the length equal to number of clusters and specify number of biclusters for each component.
#' @param model The covaraince structure you choose, there are 16 different models belongs to
#'this family:UUU, UUG, UUD, UUC, UGU, UGG, UGD, UGC, GUU, GUG, GUD, GUC, GGU, GGG, GGD, GGC. You can
#'choose more than 1 covarance structure to do model selection.
#' @param criteria one of AIC, BIC or ICL. The best model is depends on the criteria you choose. The default is BIC
#' @param iter Max iterations, defaul is 150.
#' @param permutation Only has effect when model contains UUU, UUG, UUD or UUC. If TRUE,
#'it assume the number of biclusters could be different for different components. If FALSE,
#'it assume the number of biclusters are the same cross all components. Default is FALSE.
#' @param single The default is FALSE. If TRUE, only run a single model, instead of model selection.
#' When runing the single model, \code{range_K} has to be one single number which specify the number of components,
#' and \code{rang_Q} has to be a vector that specify the number of biclusters for each component. Suggest to run with default.
#'
#'@return \code{best_model} contains a list of parameters: \code{z_ik} Estimated latent variable z,
#'\code{cluster} Component labels, \code{mu_K} Estimated component mean, \code{pie_K} Estimated component proportion,
#'\code{B_K} Estimated bicluster membership, \code{T_K} Estimated covariance of latent variable u,
#'\code{D_K} Estimated error covariance, \code{S_K} Estimated sparsity component covariance,
#'\code{LLK} Complete log likelihood value for each iteration, \code{ICL} ICL value, \code{BIC} BIC value,
#'\code{AIC} AIC value.
#'@return \code{all_fitted_model} display all names of fitted models in a data.frame.
#'
#'
#'@examples
#'
#'#if run model for UUU, UUG, UUD, or UUC
#'#bmm(data,2,range_Q=c(2:3),model="UUU")#it will run 2 combinations of Q=2 or 3
#'
#'#following will run 4 combinations of Q: 2 2, 3 3, 2 3 and 3 2
#'#bmm(data,2,range_Q=c(2:3),model="UUU", permutation=TRUE)
#'
#'#if run one model let range_Q be an integer
#'#bmm(data,2,2,model="GGC")
#'
#'#if run model selection let range_Q and range_K be a vector.
#'#bmm(data,c(2:3),c(2:3)) #use BIC as criteria to do models selection for all 16 models
#'#bmm(data,c(2:3),c(2:3),model=c("UUU","GUU"),criteria="AIC")#use AIC as criteria.
#'
#'#if run one single model
#'#bmm(data,range_K=3, range_Q=c(3,3,2),model=c("UUU"),single=T)
#'#notice that the different biclusters for each component only works for UUU, UUD, UUG and UUC model.
#'
#'#do not run:
#'#bmm(df,3,c(2,3,3),model="GUU",single = T)
#'
#'#run:
#'#bmm(df,3,c(3,3,3),model="GUU",single = T)#3 biclusters for each component
#'
#'
#'
#' @import mclust
#' @import foreach
#' @importFrom gtools permutations
#' @importFrom stringr str_split
#' @importFrom MASS ginv
#' @importFrom stats prcomp
#' @importFrom stats rmultinom
#' @importFrom stats cov
#' @importFrom stats kmeans
#'
#'@export
bmm<-function(df, range_K, range_Q, model, criteria,iter,permutation,single){

  #default for model to select
  if (missing(model)) {
    model <-
      c(
        "UUU",
        "UUG",
        "UUD",
        "UUC",
        "UGU",
        "UGG",
        "UGD",
        "UGC",
        "GUU",
        "GUG",
        "GUD",
        "GUC",
        "GGU",
        "GGG",
        "GGD",
        "GGC"
      )
  }
  else{
    model <- model
  }


  #default for criteria
  if (missing(criteria)) {
    criteria <- "BIC"
  } else{
    criteria <- criteria
  }


  #default for iteration
  if(missing(iter)) {
    iter = 150
  } else{
    iter = iter
  }

  #default for permutation
  if(missing(permutation)) {
    permutation <- FALSE
  } else{
    permutation <- permutation
  }

  #default for running model selection
  if(missing(single)) {
    single <- FALSE
  } else{
    single <- single
  }


  if(length(range_K)==1&length(range_Q)==range_K&length(model)==1&single==TRUE){
    b_cov <- model
    b_K <- range_K
    b_Q <- range_Q

    print(paste(c("Best model ",b_cov," K = ",b_K," Q= ",b_Q),collapse = ""))
    init <- try(Sim_init(df,b_K,b_Q,b_cov))
    res <- try(bi_fa_T(df,init$mu_K,init$pie_K,init$T_K,init$D_K,init$B_K,b_cov,iter))
    return(res)

  }
  else{
    all <- model_selection(df, range_K, unique(range_Q), model,permutation,iter)
    all <- all[order(all[, names(all) == criteria], decreasing = T), ]

    best <- all[which.max(all[, names(all) == criteria]), 1]
    best_split <- stringr::str_split(best, " ")

    b_cov <- best_split[[1]][1]
    b_K <- as.numeric(stringr::str_remove(best_split[[1]][2], "K"))
    b_Q <- as.numeric(stringr::str_split(best_split[[1]][3], "")[[1]][-1])

    print(paste(c("Best model ",b_cov," K = ",b_K," Q= ",b_Q),collapse = ""))}

  init <- try(Sim_init(df,b_K,b_Q,b_cov))
  res <- try(bi_fa_T(df,init$mu_K,init$pie_K,init$T_K,init$D_K,init$B_K,b_cov,iter))
  return(list(best_model=res,all_fitted_model=all))

}
