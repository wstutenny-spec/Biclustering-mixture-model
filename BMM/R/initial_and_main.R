#' Biclustering mixture model algorithm
#'
#' Main function that can do initialize biclustering mixture model and run one particular model for BMM
#' @param df A data frame that tries to bicluster. Should only contain numeric values.
#'
#' @param K One possible number of components.
#' @param Q_K A vector specify the number of bicluster for each component.
#' @param pie_K A vector of initial guesses of component proportion
#' @param mu_K A list of initial guess of mean vector
#' @param B_K A list of initial guess of bicluster membership for each component
#' @param T_K A list of initial guess of covariance of latent variable for each component
#' @param D_K A list of initial guess of error matrix for each component
#' @param cov_str The covaraince structure you choose, there are 16 different models belongs to
#'this family:UUU, UUG, UUD, UUC, UGU, UGG, UGD, UGC, GUU, GUG, GUD, GUC, GGU, GGG, GGD, GGC. You can
#'choose more than 1 covarance structure to do model selection.
#' @param iter Max iterations.
#'
#'
#'
#'@return z_ik Estimated latent variable z
#'@return cluster Component labels
#'@return mu_K Estimated component mean
#'@return pie_K Estimated component proportion
#'@return B_K Estimated bicluster membership
#'@return T_K Estimated covariance of latent variable u
#'@return D_K Estimated error covariance
#'@return S_K Estimated sparsity component covariance
#'@return LLK Complete log likelihood value for each iteration
#'@return ICL ICL value
#'@return BIC BIC value
#'@return AIC AIC value



####AECM function####
bi_fa_T<-function(df,mu_K,pie_K,T_K,D_K,B_K,cov_str,iter){

  p <- dim(df)[2]
  N <- dim(df)[1]
  K<-length(pie_K)
  #formate T_K
  T_K<-lapply(T_K, as.matrix)


  Q_K<-unlist(lapply(T_K, ncol))
  #U_K #basicly factors
  U_iK<-list()
  #new sigma
  S_K<-list()


  #try true labs
  # for (i in 1:N) {
  #   w_K[i,1]<-ifelse(lab[i]==1,1,0)
  #   w_K[i,2]<-ifelse(lab[i]==2,1,0)
  # }
  #

  ####based on different covariance structure, get each estimation####
  code<-strsplit(cov_str,"")[[1]]
  code_B<-code[1]
  code_T<-code[2]
  code_D<-code[3]

  if((code_T=="G"|code_T=="C")&length(unique(Q_K))!=1)
  {return("incorrect match of bicluster dimension and cov structure for T")}
  else if(code_B=="G"&length(unique(Q_K))!=1)
  {return("incorrect match of bicluster dimension and cov structure for B")}
  else{

    #initialize first guess mix density
    #use original D_K, B_K and initial pie_K, mu_K to get mix_density
    #mix_density
    mix_density <- matrix(0, nrow = N)
    for (i in 1:K) {
      COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
      mix_density <- mix_density + pie_K[i] * dmvnorm(df, mu_K[[i]],COV)
    }

    overall_loglik<-rep(-Inf,2)

    h<-2
    #overall cycle
    repeat{
      ####begin bicluster first cycle####
      #w_K
      w_K <- matrix(0, ncol = K, nrow = N)
      for (i in 1:K) {
        COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
        w_K[, i] <- pie_K[i] * dmvnorm(df, mu_K[[i]], COV)/mix_density
      }


      #pie_K
      pie_K<-colMeans(w_K)
      label<-unlist(apply(w_K, 1, which.max))
      if(length(pie_K[pie_K==0])>0|any(table(label)<3)){break}

      #mu_K
      for (i in 1:K) {
        mu_K[[i]] <- t(w_K[, i]) %*% as.matrix(df) / sum(w_K[, i])
      }

      o_g<-order(pie_K,decreasing = T)
      mu_K<-mu_K[o_g]
      pie_K<-pie_K[o_g]
      w_K<-w_K[,o_g]
      B_K<-B_K[o_g]
      D_K<-D_K[o_g]
      T_K<-T_K[o_g]
      Q_K<-Q_K[o_g]


      #mix_density
      mix_density <- matrix(0, nrow = N)
      for (i in 1:K) {
        COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
        mix_density <- mix_density + pie_K[i] * dmvnorm(df, mu_K[[i]],COV)
      }

      ####second cycle####
      #from now on, don't need to update mu_K and pie_K
      #w_K
      w_K <- matrix(0, ncol = K, nrow = N)
      for (i in 1:K) {
        COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
        w_K[, i] <- pie_K[i] * dmvnorm(df, mu_K[[i]], COV) / mix_density
      }
      pie_K<-colMeans(w_K)
      label<-unlist(apply(w_K, 1, which.max))
      if(length(pie_K[pie_K==0])>0|any(table(label)<3)){break}


      o_g<-order(pie_K,decreasing = T)
      mu_K<-mu_K[o_g]
      pie_K<-pie_K[o_g]
      w_K<-w_K[,o_g]
      B_K<-B_K[o_g]
      D_K<-D_K[o_g]
      T_K<-T_K[o_g]
      Q_K<-Q_K[o_g]


      #S_K
      w_K<-as.matrix(w_K,K,N)
      for (i in 1:K) {
        S_K[[i]] <- (t(df) - c(mu_K[[i]])) %*% diag(w_K[, i] / sum(w_K[, i])) %*%
          t(t(df) - c(mu_K[[i]]))
      }


      #U_iK
      for (i in 1:K) {
        U_iK[[i]]<-T_K[[i]]%*%t(B_K[[i]])%*%
          solvecov(B_K[[i]],T_K[[i]],D_K[[i]])%*%(t(df)-c(mu_K[[i]]))
      }

      #e2
      e2<-list()
      for (i in 1:K) {
        e2[[i]]<-T_K[[i]]%*%t(B_K[[i]])%*%solvecov(B_K[[i]],T_K[[i]],D_K[[i]])%*%
          S_K[[i]]%*%t(solvecov(B_K[[i]],T_K[[i]],D_K[[i]]))%*%B_K[[i]]%*%T_K[[i]]+
          T_K[[i]]-T_K[[i]]%*%t(B_K[[i]])%*%solvecov(B_K[[i]],T_K[[i]],D_K[[i]])%*%
          B_K[[i]]%*%T_K[[i]]
      }

      #T_K
      #E(uu^T|y)
      for (i in 1:K) {
        T_K[[i]]<-diag(diag(e2[[i]]),Q_K[i])

      }

      p_T<-sum(Q_K)
      int_T_K<-T_K

      if(code_T=="G"&length(unique(Q_K))==1){
        p_T<-Q_K[1]
        beg_T<-0
        for (i in 1:K) {
          beg_T<-int_T_K[[i]]*pie_K[i]+beg_T
        }
        for (i in 1:K) {
          T_K[[i]]<-beg_T
        }
      }
      else if(code_T=="D"){
        p_T<-K
        for (i in 1:K) {
          T_K[[i]]<-diag(mean(diag(int_T_K[[i]])),Q_K[i])
        }
      }
      else if(code_T=="C"&length(unique(Q_K))==1){
        p_T<-1
        beg_T<-0
        for (i in 1:K) {
          beg_T<-int_T_K[[i]]*pie_K[i]+beg_T
        }
        for (i in 1:K) {
          T_K[[i]]<-diag(mean(diag(beg_T)),Q_K[i])
        }
      }




      #D_K
      for (i in 1:K) {
        # D_K[[i]]<-diag(diag(S_K[[i]]-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])%*%
        #                      solvecov(B_K[[i]],T_K[[i]],D_K[[i]])%*%S_K[[i]]))


        D_K[[i]]<-diag(diag(S_K[[i]]-2*B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])%*%
                             solvecov(B_K[[i]],T_K[[i]],D_K[[i]])%*%S_K[[i]]+
                             B_K[[i]]%*%e2[[i]]%*%t(B_K[[i]])))


      }

      p_D<-K*p


      int_D_K<-D_K
      if(code_D=="G"){
        p_D<-p
        beg_D<-0
        for (i in 1:K) {
          beg_D<-int_D_K[[i]]*pie_K[i]+beg_D
        }
        for (i in 1:K) {
          D_K[[i]]<-beg_D
        }
      }
      else if(code_D=="D"){
        p_D<-K
        for (i in 1:K) {
          D_K[[i]]<-diag(mean(diag(int_D_K[[i]])),p)
        }
      }
      else if(code_D=="C"){
        p_D<-1
        beg_D<-0
        for (i in 1:K) {
          beg_D<-int_D_K[[i]]*pie_K[i]+beg_D
        }
        for (i in 1:K) {
          D_K[[i]]<-diag(mean(diag(beg_D)),p)
        }
      }




      pre_B_K<-B_K

      #change B_K
      #H2, loglikelyhood
      if(code_B=="U"){
        p_B<-K*p
        for(i in 1:K){
          for(r in 1:p){
            x<-b(r,B_K[[i]])

            #use lapply to calculate
            track_H2<-lapply(x,function(x,i){

              invd<-diag(diag(D_K[[i]])^(-1))
              invcov<-solvecov(x,T_K[[i]],D_K[[i]])
              H2<-tr(invd%*%x%*%T_K[[i]]%*%t(x)%*%invcov%*%S_K[[i]])-
                tr(invd%*%x%*%e2[[i]]%*%t(x))/2
              return(ifelse(is.na(H2),-Inf,H2))
            },i=i)



            B_K[[i]]<-x[[which.max(track_H2)]]


            # mix_density <- matrix(0, nrow = N)
            # for (j in 1:K) {
            #   COV<-B_K[[j]]%*%T_K[[j]]%*%t(B_K[[j]])+D_K[[j]]
            #   mix_density <- mix_density + pie_K[j] * dmvnorm(df, mu_K[[j]],COV)
            # }
            # print(paste("B_K",i,r,sum(log(mix_density)),track_H2))
            # print(B_K[[i]])
            # print(track_H2)

          }

        }
      }
      else if(code_B=="G"&length(unique(Q_K))==1){
        p_B<-p
        for(r in 1:p){
          x<-b(r,B_K[[1]])
          track_H2<-lapply(x,function(x){
            H2<-0
            for(i in 1:K){
              invd<-diag(diag(D_K[[i]])^(-1))
              invcov<-solvecov(x,T_K[[i]],D_K[[i]])
              H2<-H2+sum(w_K[,i])*(tr(invd%*%x%*%T_K[[i]]%*%t(x)%*%invcov%*%S_K[[i]])-
                                    tr(invd%*%x%*%e2[[i]]%*%t(x))/2)
            }
            return(ifelse(is.na(H2),-Inf,H2))
          })
          B_K[[1]]<-x[[which.max(track_H2)]]
          for(i in 1:K){
            B_K[[i]]<-B_K[[1]]
          }
          #print(c(r,unlist(track_H2)))
        }
      }


      if(code_B=="G"){
        mix_density <- matrix(0, nrow = N)
        for (i in 1:K) {
          COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
          mix_density <- mix_density + pie_K[i] * dmvnorm(df, mu_K[[i]],COV)
        }
        if(sum(log(mix_density[mix_density!=0]))>overall_loglik[h-1]){B_K<-B_K}
        else{B_K<-pre_B_K}
      }else if(code_B=="U"){
        for (i in 1:K) {
          if(all(as.vector(pre_B_K[[i]]%*%t(pre_B_K[[i]]))==as.vector(B_K[[i]]%*%t(B_K[[i]])))){
            B_K[[i]]<-pre_B_K[[i]]
          }else{B_K[[i]]<-B_K[[i]]}
        }
      }


      # print(paste("B_K",r,sum(log(mix_density)),sum(mix_density==0)))
      # print(colSums(B_K[[1]]))


      #calculate overall log likelyhood.
      mix_density <- matrix(0, nrow = N)
      for (i in 1:K) {
        COV<-B_K[[i]]%*%T_K[[i]]%*%t(B_K[[i]])+D_K[[i]]
        mix_density <- mix_density + pie_K[i] * dmvnorm(df, mu_K[[i]],COV)
      }
      overall_loglik[h]<-sum(log(mix_density[mix_density!=0]))
      if(is.na(overall_loglik[h])){overall_loglik[h]<-overall_loglik[h-1]}


      if(h>3){
        a<-(overall_loglik[h]-overall_loglik[h-1])/(overall_loglik[h-1]-overall_loglik[h-2])
        L_inf<-overall_loglik[h-1]+(overall_loglik[h]-overall_loglik[h-1])/(1-a)}else{L_inf<-Inf}

      diff<-L_inf-overall_loglik[h]
      if(is.na(diff)){diff<-Inf}
      if(abs(diff)<10^(-2)|h>iter|length(pie_K[pie_K==0])>0){break}
      h<-h+1
    }

    if(length(unique(label))<K|any(table(label)<3)|length(pie_K[pie_K==0])>0)
    {return("fewer class than defined")}
    else{
      m<-matrix(0,nrow = N,ncol=K)
      for (i in 1:N) {
        m[i,label[i]]<-1
      }
      list(
        #estimators
        z_ik=w_K,
        cluster=label,
        mu_K=mu_K,
        pie_K=pie_K,
        S_K=S_K,
        B_K=B_K,
        T_K=T_K,
        D_K=D_K,
   #loglikelyhood
        LLK=overall_loglik[-1],
        #bic
        ICL=overall_loglik[h]*2-log(N)*(K*p+K-1+p_B+p_D+p_T)+2*sum(as.vector(m)*log(as.vector(w_K))),
        BIC=overall_loglik[h]*2-log(N)*(K*p+K-1+p_B+p_T+p_D),
        AIC=overall_loglik[h]*2-2*(K*p+K-1+p_B+p_T+p_D)
      )
    }
  }
}



####Initial guess####
Sim_init <- function(df,K,Q_K,cov_str) {

  p <- dim(df)[2]
  N <- dim(df)[1]
  mu_K <- list()
  sig_K <- list()
  B_K<-list()
  #D_K, prespecified error covariance matrix
  D_K<-list()
  #T_K, diag cov matrix for factor u
  T_K<-list()

  code<-strsplit(cov_str,"")[[1]]
  code_B<-code[1]
  code_T<-code[2]
  code_D<-code[3]

  if((code_T=="G"|code_T=="C")&length(unique(Q_K))!=1)
  {return("incorrect match of bicluster dimension and cov structure for T")}
  else if(code_B=="G"&length(unique(Q_K))!=1)
  {return("incorrect match of bicluster dimension and cov structure for B")}
  else{

    ####initial values####
    res<-Mclust(df,K)
    #res<-kmeans(df,K)
    #class label
    #z<-lab
    z<-res$classification
    #z<-res$cluster
    if(length(unique(z))<K|any(table(z)<=Q_K)==TRUE){
      return("fewer class than defined")
    }
    else{
      pie_K<-table(z)/N
      o_g<-order(pie_K,decreasing = T)

      #order group
      # for (i in 1:N) {
      #   z[i]<-c(1:K)[o_g==z[i]]
      #   }

      #pie_K<-pie_K[o_g]

      #set initial for each parameter
      for (i in 1:K) {
        mu_K[[i]]<-colMeans(df[z==i,])
      }



      for(i in 1:K){
        sig_K[[i]]<-cov(df[z == i,])
      }

      #initial for T_K
      for (i in 1:K){
        fit<-prcomp(df[z==i,])
        T_K[[i]]<-diag(fit$sdev[1:Q_K[i]]^2,Q_K[i])
      }

      int_T_K<-T_K
      if(code_T=="G"&length(unique(Q_K))==1){
        beg_T<-0
        for (i in 1:K) {
          beg_T<-int_T_K[[i]]*pie_K[i]+beg_T
        }
        for (i in 1:K) {
          T_K[[i]]<-beg_T
        }
      }
      else if(code_T=="D"){
        for (i in 1:K) {
          T_K[[i]]<-diag(mean(diag(int_T_K[[i]])),Q_K[i])
        }
      }
      else if(code_T=="C"&length(unique(Q_K))==1){
        beg_T<-0
        for (i in 1:K) {
          beg_T<-int_T_K[[i]]*pie_K[i]+beg_T
        }
        for (i in 1:K) {
          T_K[[i]]<-diag(mean(diag(beg_T)),Q_K[i])
        }
      }




      #initial for D_K
      for (i in 1:K) {
        fit<-prcomp(df[z==i,])
        D_K[[i]]<-diag(diag(sig_K[[i]]-fit$rotation[,1:Q_K[i]]%*%diag(fit$sdev[1:Q_K[i]]^2,Q_K[i])%*%
                              t(fit$rotation[,1:Q_K[i]])))

      }
      int_D_K<-D_K
      if(code_D=="G"){
        beg_D<-0
        for (i in 1:K) {
          beg_D<-int_D_K[[i]]*pie_K[i]+beg_D
        }
        for (i in 1:K) {
          D_K[[i]]<-beg_D
        }
      }
      else if(code_D=="D"){
        for (i in 1:K) {
          D_K[[i]]<-diag(mean(diag(int_D_K[[i]])),p)
        }
      }
      else if(code_D=="C"){
        beg_D<-0
        for (i in 1:K) {
          beg_D<-int_D_K[[i]]*pie_K[i]+beg_D
        }
        for (i in 1:K) {
          D_K[[i]]<-diag(mean(diag(beg_D)),p)
        }
      }

      #initial for B_K
      for(i in 1:K){
        fit<-prcomp(scale(df[z==i,]))
        bk<-fit$rotation[,1:Q_K[i]]
        list_result <- lapply(split(bk,seq(NROW(bk))),function(x,i){
          y<-rep(0,Q_K[i])
          y[which.max(x)]<-1
          y
        },i=i)
        B_K[[i]]<- do.call(rbind,list_result)

      }
      if(code_B=="G"&length(unique(Q_K))==1){
        for (i in 1:K) {
          B_K[[i]]<-B_K[[1]]
        }
      }


      #reorder
      # for (i in 1:K) {
      #   if(Q_K[i]!=1){
      #   index<-apply(B_K[[i]],2,sum)
      #   B_K[[i]]<-B_K[[i]][,order(index,decreasing = T)]
      #   T_K[[i]]<-diag(diag(T_K[[i]])[order(index,decreasing = T)])
      #   }
      # }

      list(mu_K=mu_K,
           pie_K=pie_K,
           T_K=T_K,
           D_K=D_K,
           B_K=B_K)
    }
  }
}

####invert factor structure covariance matrix####
solvecov<-function(b,t,d){
  p<-dim(t)[1]
  invd<-diag(diag(d)^(-1))
  bt<-b%*%t^(1/2)
  invd-invd%*%bt%*%solve(diag(1,p,p)+t(bt)%*%invd%*%bt)%*%t(bt)%*%invd
}






####trace function####
tr<-function(x){sum(diag(x))}



####update B_K function####
b<-function(r,B_K){
  #define the number of column for B_K
  q<-length(B_K[r,])
  B_test<-vector('list',q)
  #create q identical B_K
  for(e in 1:q){
    B_test[[e]]<-B_K
    B_test[[e]][r,]<-0
    B_test[[e]][r,e]<-1
  }
  B_test
}




# #function to store parameter to table
# store_table<-function(para,name,procid){
#   output.filename <- paste("simulationOutput_",name,".csv",sep='')
#   write.table(t(para),file=output.filename,row.names=paste("simulation",procid),
#               col.names=F,sep=",",append = TRUE)
# }



# #function start randomly generate B matrix
# rand_B<-function(p,Q_K){
#   B_K<-list()
#   for (i in 1:length(Q_K)) {
#     B_K[[i]]<-t(rmultinom(p,1,rep(1/Q_K[i],Q_K[i])))
#   }
#   B_K
# }










