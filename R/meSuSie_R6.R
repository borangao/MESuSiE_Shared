require(progress)
require(R6)
#' meSuSie core function
#'
#' @param  R_mat_list A list of length N ancestry with each elment being correlation matrix with dimension p*p, the column name of the correlation matrix should match to the order of SNP name in summary_stat_list
#' @param  summary_stat_list A list of length N ancestry with each elment being summary statistics. The minimum requirement of summary statistics contains columns of SNP, Beta, Se, Z, and N. The order of the SNP should match the order of the correlation matrix 
#' @param  L Number of effects supposed within the region
#' @return An R6 object with pip, credible sets, and other features of the fine-mapping result. 
#' @export
meSuSie_core<-function(R_mat_list,summary_stat_list,L,residual_variance=NULL,prior_weights=NULL,optim_method ="optim",estimate_residual_variance =F,max_iter =100){
  
  cat("*************************************************************\n
  Multiple Ancestry Sum of Single Effect Model (meSuSie)          \n
   Visit http://www.xzlab.org/software.html For Update            \n
            (C) 2022 Boran Gao, Xiang Zhou                        \n
              GNU General Public License                          \n
*************************************************************") 
  
  time_start<-Sys.time()
  cat("\n# Start data processing for sufficient statistics \n")
  meSuSieData_obj<-meSuSieData$new(R_mat_list,summary_stat_list)
  
  n_snp = nrow(summary_stat_list[[1]])
  n_ancestry = length(summary_stat_list)
  if(is.null(prior_weights)){
    prior_weights = rep(1/n_snp,n_snp)
  }
  if(is.null(residual_variance)){
    residual_variance = rep(1,n_ancestry)
  }
  cat("# Create meSuSie object \n")
  meSuSieObject_obj<-meSuSieObject$new(n_snp,L,n_ancestry,residual_variance,prior_weights,optim_method,estimate_residual_variance,max_iter)
  
  cat("# Start data analysis \n")
    pb = progress_bar$new(format = paste("(Iteration = :iteration) :elapsed"),clear = TRUE,total = max_iter,show_after = 0)
  
  n_iter = 0
  for (iter in 1:max_iter) {
    
    comp_residual<-meSuSieObject_obj$compute_residual(meSuSieData_obj,meSuSieObject_obj)
    
    
    for (l_index in seq(1,L,1)) {
      
      comp_residual = comp_residual + t(meSuSieObject_obj$Xr[[l_index]])
      
      SER_res<-single_effect_regression(comp_residual,meSuSieData_obj$XtX.diag, meSuSieObject_obj,l_index) 
      
      meSuSieObject_obj$par_update(SER_res,l_index)
      
      meSuSieObject_obj$compute_KL(SER_res,meSuSieObject_obj$compute_SER_posterior_loglik(meSuSieData_obj,comp_residual,l_index),l_index)
      
      meSuSieObject_obj$compute_Xr(meSuSieData_obj,SER_res,l_index)
      
      comp_residual = comp_residual - t(meSuSieObject_obj$Xr[[l_index]])
      
      
    }

    pb$tick(tokens = list(iteration = iter))
    
    updated_sigma2 = meSuSieObject_obj$update_residual_variance(meSuSieData_obj,iter)
    
    if((meSuSieObject_obj$ELBO[iter+1] - meSuSieObject_obj$ELBO[iter])<0.001){
      break
    }
    if(meSuSieObject_obj$estimate_residual_variance==TRUE){
      meSuSieObject_obj$sigma2 =  updated_sigma2
    }
    n_iter = n_iter + 1
    #check convergence and update sigma2
  }

  cat("\n# Data analysis is done, and now generates result \n\n")
  ###Use function in Utility to output result
  meSuSieObject_obj$get_result(meSuSie_get_cs(meSuSieObject_obj,R_mat_list),meSusie_get_pip(meSuSieObject_obj))
  meSuSieObject_obj$mesusie_summary(meSuSieData_obj)

  time_end<-Sys.time()
  cat(c("\n# Total time used for the analysis:",paste0(round(as.numeric(difftime(time_end,time_start,units=c("mins"))),2)," mins\n")))
  return(meSuSieObject_obj)
}


meSuSieObject <- R6Class("meSuSieObject",public = list(
  initialize = function(p,L,N_ancestry, residual_variance,prior_weights,estimate_prior_method,estimate_residual_variance,max_iter){
    
    self$alpha = matrix(1/p,nrow = L,ncol = p)
    self$mu1 = rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
    self$mu2 = rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
    self$Xr     = rep(list(matrix(0,nrow = N_ancestry,ncol=p)),L)
    
    self$KL     = rep(as.numeric(NA),L)
    self$lbf    = rep(as.numeric(NA),L)
    self$lbf_variable = vector("list",L)
    self$ELBO = rep(NA,max_iter)
    self$ELBO[1] = -Inf
    self$sigma2 = residual_variance
    
    
    self$V      = rep(list(matrix(0,ncol = N_ancestry,nrow=N_ancestry)),L)
    self$pi     = prior_weights
    
    self$estimate_prior_method = estimate_prior_method
    self$estimate_residual_variance = estimate_residual_variance
    
    self$L = L
    self$nSNP = p
    self$nancestry = N_ancestry
    
    self$cs = list()
    self$pip = rep(as.numeric(NA),p)
  },    
  
  compute_Xb = function (meSuSie_Data,b){
    lapply(1:length(meSuSie_Data$XtX_list),function(x){
      meSuSie_Data$XtX_list[[x]]%*%b[x,]
    })},
  
  compute_residual = function (meSuSie_Data,meSuSie_Obj) {
    residual<-Reduce(cbind,meSuSie_Data$Xty_list)-t(Reduce("+",meSuSie_Obj$Xr))
    return(residual)
  },
  compute_Xr = function(meSuSie_Data,SER_res,l_index){
    self$Xr[[l_index]] = t(Reduce(cbind,lapply(1:self$nancestry,function(ancestry_index)meSuSie_Data$XtX_list[[ancestry_index]]%*%(SER_res$mu1_multi[,ancestry_index]*SER_res$alpha))))
  },
  
  par_update = function(SER_res,l_index){
    self$alpha[l_index,]<-SER_res$alpha
    self$mu1[[l_index]]<-SER_res$mu1_multi
    self$mu2[[l_index]]<-SER_res$mu2_multi
    self$lbf_variable[[l_index]]<-SER_res$lbf_multi
    self$V[[l_index]]<-SER_res$V
  },
  compute_SER_posterior_loglik = function(meSuSie_Data,comp_residual,l_index){
    return(Reduce("+",lapply(1:self$nancestry,function(x){
      sum(-0.5/self$sigma2[[x]]*(-2*comp_residual*self$mu1[[l_index]][,x]*self$alpha[l_index,]+meSuSie_Data$XtX.diag[[x]]*self$mu2[[l_index]][,x]*self$alpha[l_index,]))
    })))
    
  },
  compute_KL = function(SER_res,value,l_index){
    
    self$KL[l_index] = -SER_res$loglik+value
    
  },
  
  update_residual_variance = function(meSuSie_Data,niter){
    
    B_1 = lapply(1:self$L,function(x){
      apply(self$mu1[[x]],2,function(y)y*t(self$alpha[x,]))
    })
    
    B_1_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(B_1,function(x)x[,y])))
    
    
    BXXB = lapply(1:self$nancestry,function(x){
      sum((t(B_1_ancestry[[x]])%*% meSuSie_Data$XtX_list[[x]])*t(B_1_ancestry[[x]]))
    })  ####{E(BX)}^2
    
    betabar = Reduce("+",B_1)
    
    BbarXXBbar = lapply(1:self$nancestry,function(x){
      sum(betabar[,x]*(meSuSie_Data$XtX_list[[x]]%*%betabar[,x]))
    }) ###{E(BbarX)}^2
    
    BbarXty = lapply(1:self$nancestry,function(x){
      2*sum(betabar[,x]*meSuSie_Data$Xty_list[[x]])
    }) ###{E(BbarX)}^2
    
    B_2 = lapply(1:self$L,function(x){
      apply(self$mu2[[x]],2,function(y)y*t(self$alpha[x,]))
    })
    
    B_2_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(B_2,function(x)x[,y])))
    
    XXB2 = lapply(1:self$nancestry,function(x){
      sum(meSuSie_Data$XtX.diag[[x]]*B_2_ancestry[[x]])
    })
    
    updated_sigma = lapply(1:self$nancestry,function(x){
      
      ((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))/meSuSie_Data$N_list[[x]]
    })
    
    
    self$ELBO[niter+1] = Reduce(sum,lapply(1:self$nancestry,function(x){
      -meSuSie_Data$N_list[[x]]/2*log(2*pi*self$sigma2[x])-(1/(2*self$sigma2[x]))*((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))
    }))-sum(self$KL)
    
    return( unlist(updated_sigma))
  },
  
  get_result = function(cs, pip){
    self$cs = cs
    self$pip = pip
  },
  mesusie_summary = function(meSuSie_Data){
    
    cat(c(paste0("Potential causal SNPs with PIP > 0.5: "),meSuSie_Data$Summary_Stat[[1]]$SNP[which(self$pip>0.5)],"\n\n"))
    cat("Credible sets for effects: \n")
    print(self$cs)
    cat("\n Use meSusie_plot_pip() for Mahattan and PIP Plot")

  }
  
),lock_objects = F
)


##############################################################
#  Assume var(y) = 1 and causal snps contribute negligible variance
#  therefore se(\beta) = var(y)*diag(xtx) => diag(xtx) = var(y)/se(beta)
#                          or
#  R^2 = z^2/(z^2+N-2) => sigma^2 = var(y)*(N-1)/(z^2+N-2)=>diag(xtx)=sigma^2/se(beta)^2
#
#  xtx = sqrt(diag(xtx))%*%R%*%sqrt(diag(xtx))
#  xty = diag(xtx)\beta
#
#############################################################
meSuSieData <- R6Class("meSuSieData",public = list(
  initialize = function(X,Y,var_y = 1){
    self$R<- X
    self$Summary_Stat <-Y
    self$var_y = var_y
    
    
    self$Name_list<-as.list(names(self$R))
    names(self$Name_list)<-names(self$R)
    self$N_ancestry<-length(X)
    
    self$XtX.diag<-self$XtX_diag(self$Summary_Stat,self$Name_list)
    self$XtX_list<-self$XtX_pro(self$R,self$XtX.diag,self$Name_list)
    self$Xty_list<-self$Xty_pro(self$Summary_Stat,self$XtX.diag,self$Name_list) ##diag(xtx)^*betahat
    
    self$N_list<-lapply(self$Summary_Stat,function(x)median(x$N))
    self$yty_list<-lapply(self$N_list,function(x)return(self$var_y*(x-1)))
    
    return(self)},
  ##First compute XtX.diag
  XtX_diag = function(Summary_Stat,Name_list){
    return( lapply(Name_list,function(x){
      R2 = (summary_stat_list[[x]]$Z^2)/(summary_stat_list[[x]]$Z^2+summary_stat_list[[x]]$N-2)
      sigma2 = self$var_y*(1-R2)*(summary_stat_list[[x]]$N-1)/(summary_stat_list[[x]]$N-2)
      return(sigma2/(Summary_Stat[[x]]$Se)^2)
    }))
  },
  
  
  ###Process XtX
  XtX_pro = function(R,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      return(diag(sqrt(XtX.diag[[x]]))%*%R[[x]]%*%diag(sqrt(XtX.diag[[x]])))
    }))
  },
  
  ###Process XtY
  Xty_pro = function(Summary_Stat,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      XtX.diag[[x]]*Summary_Stat[[x]]$Beta
    }))
    
  }
),
lock_objects = F)

single_effect_regression<-function(XtR,XtX.diag, meSuSieObject_obj,l_index){
  
  N_ancestry = meSuSieObject_obj$nancestry
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))
  
  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))
  
  betahat = shat2 * Xty_standardized
  
  if(meSuSieObject_obj$estimate_prior_method =="optim"){
    
    opt_par<-pre_optim(N_ancestry,-30,10)
    
    update_V<-optim(opt_par$inital_par,fn = loglik_cpp_R6,gr=NULL,betahat,shat2,meSuSieObject_obj$pi,opt_par$nancestry,opt_par$diag_index,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
  }
  
  V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
  
  #multivariate_out<-multivariate_regression(Xty_standardized,shat2,V_mat) 
  multivariate_out<-mvlmm_reg(betahat,shat2,V_mat) 
  
  lbf<-multivariate_out$lbf
  lbf[is.na(lbf)]<-0
  
  softmax_out<-compute_softmax(lbf,meSuSieObject_obj$pi)
  
  return(list(alpha = softmax_out$alpha_wmulti,mu1_multi = multivariate_out$post_mean,mu2_multi = multivariate_out$post_mean2,lbf_multi = lbf,V = V_mat,loglik =softmax_out$loglik))
}


