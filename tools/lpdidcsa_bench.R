# ______________________________________________________________________________
# Lpdidcsa timing Script -------------------------------------------------------
# ______________________________________________________________________________

## Libraries ----
library(bench)
library(data.table)
library(lpdidcsa)
library(did)
library(fixest)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Small database ----

all <- NULL


# Creation of a toy data 
#Number of individuals
n_indiv <- 5000L

#Number of firms
n_firms <- 100L


toydata <- sim_staggered_panel(n=n_indiv,n_firms=n_firms,
                               t_min_treat = 5L,
                               t_max_treat = 15L)

format(object.size(toydata),units="Mb")
toydata <- toydata[,list(id,t,firm_id,female,treat,log_earnings,g_i)]
toydata[,g_i:=fifelse(is.na(g_i),0,g_i)]



n_it <- 10

## Format wide -------------------------------
par_type_horizon <- "wide"
par_one_reg <- F



## Wide horizon database ---------------------
result <- bench::mark({df  <- lpdidcsa_data(toydata, 
                                          unit = "id", 
                                          time = "t",
                                          dependent = "log_earnings",
                                          treat = "treat",
                                          type_horizon= par_type_horizon
)},
iterations = n_it,
check = FALSE  
)

size_df <- format(object.size(df),units="Mb")
size_df
colnames(df)
nrow(df)

speed <- data.table(
  method = "data wide",
  format = "wide",
  one_reg = NA,
  median_time = result$median,  
  min_time = min(unlist(result$time)),
  max_time = max(unlist(result$time)),
  iterations = result$n_itr,
  n_rows=nrow(df),
  n_cols=ncol(df),
  size=size_df)
speed

all <- rbind(all,speed,fill=T)
write.csv(all,"bench.csv",row.names = F)


# Function for benchmarking the methods with  bench::mark & tryCatch
bench_method <- function(df, method, par_type_horizon, par_one_reg, control = "female", min_iterations = 3) {
  tryCatch({
    result <- bench::mark(
      {
        res <- lpdidcsa(
          data=df, 
          unit = "id", time = "t",
          dependent = "log_earnings",
          meth = method,
          controls = control,
          type_horizon = par_type_horizon,
          one_reg = par_one_reg
        )
      },
      iterations = min_iterations,
      check = FALSE,
      memory = FALSE
    )
    
    # Simplified summary
    list(
      median_time = result$median,  
      min_time = min(unlist(result$time)),
      max_time = max(unlist(result$time)),
      iterations = result$n_itr   
    )
  }, error = function(e) {
    message(paste("Error for the method", method, ":", e$message))
    return(NULL)
  })
  
}



## Wide ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df <- rbindlist(results, idcol = "method")
results_df$format <- "wide"
results_df$one_reg <- FALSE

all <- rbind(all,results_df,fill=T)
all
write.csv(all,"bench.csv",row.names = F)
gc()

### Wide method adj----
methods <- c("lpdid_adj")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",  
    par_one_reg = FALSE,
    min_iterations = 2
  )
}

# Display results
results_df_adj <- rbindlist(results, idcol = "method")
print(results_df_adj[, .(method, median_time,  iterations)])

results_df_adj$format <- "wide"
results_df_adj$one_reg <- FALSE

all <- rbind(all,results_df_adj,fill=T)
all
write.csv(all,"bench.csv",row.names = F)
gc()


### Wide interacted ----
methods <- c("lpdid_ipw","lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",  # ou "long"
    control="female:i(t)",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df_int <- rbindlist(results, idcol = "method")
print(results_df_int[, .(method, median_time,  iterations)])

results_df_int$format <- "wide"
results_df_int$one_reg <- FALSE

all <- rbind(all,results_df_int,fill=T)
all

write.csv(all,"bench.csv",row.names = F)

## CSA ipw ----
time_csa_ipw <- bench::mark({csa_ipw_det <- att_gt(yname="log_earnings",
                                                   tname="t",
                                                   idname="id",
                                                   gname="g_i",
                                                   base_period="universal",
                                                   control_group = c("notyettreated"),
                                                   est_method="ipw",
                                                   bstrap = FALSE,
                                                   cband = FALSE,
                                                   xformla=~female,
                                                   data=toydata)
csa_ipw <- aggte(csa_ipw_det, type = "dynamic")},
iterations = n_it,
check = FALSE,
memory = FALSE
)    
rm(csa_ipw,csa_ipw_det)

time_csa_ipw

speed <- data.table(
  method = "csa ipw",
  format = NA,
  one_reg = NA,
  median_time = time_csa_ipw$median,  
  min_time = min(unlist(time_csa_ipw$time)),
  max_time = max(unlist(time_csa_ipw$time)),
  iterations = time_csa_ipw$n_itr)
speed

all <- rbind(all,speed,fill=T)
all
write.csv(all,"bench.csv",row.names = F)

## CSA reg ----
time_csa_reg <- bench::mark({csa_reg_det <- att_gt(yname="log_earnings",
                                                   tname="t",
                                                   idname="id",
                                                   gname="g_i",
                                                   base_period="universal",
                                                   control_group = c("notyettreated"),
                                                   est_method="reg",
                                                   bstrap = FALSE,
                                                   cband = FALSE,
                                                   xformla=~female,
                                                   data=toydata)
csa_reg <- aggte(csa_reg_det, type = "dynamic")},
iterations = n_it,
check = FALSE,
memory = FALSE)    
rm(csa_reg,csa_reg_det)

speed <- data.table(
  method = "csa reg",
  format = NA,
  one_reg = NA,
  median_time = time_csa_reg$median,  
  min_time = min(unlist(time_csa_reg$time)),
  max_time = max(unlist(time_csa_reg$time)),
  iterations = time_csa_reg$n_itr)
speed

all <- rbind(all,speed,fill=T)
all
write.csv(all,"bench.csv",row.names = F)



gc()

rm(df)
## Long horizon database ---------------------
result <- bench::mark({df_l  <- lpdidcsa_data(toydata, 
                                              unit = "id", 
                                              time = "t",
                                              dependent = "log_earnings",
                                              treat = "treat",
                                              type_horizon= "long"
)},
iterations = n_it,
check = FALSE  # Désactive la vérification des résultats identiques
)

size_df_l <- format(object.size(df_l),units="Mb")
size_df_l
colnames(df_l)

speed <- data.table(
  method = "data long",
  format = "long",
  one_reg = NA,
  median_time = result$median,  
  min_time = min(unlist(result$time)),
  max_time = max(unlist(result$time)),
  iterations = result$n_itr,
  n_rows=nrow(df_l),
  n_cols=ncol(df_l),
  size=size_df_l)
speed

all <- rbind(all,speed,fill=T)
all

write.csv(all,"bench.csv",row.names = F)



## Long multiple regressions ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking méthode:", method, "\n")
  results[[method]] <- bench_method(
    df = df_l,
    method = method,
    par_type_horizon = "long",  # ou "long"
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df_l <- rbindlist(results, idcol = "method")
print(results_df_l[, .(method, median_time,  iterations)])

results_df_l$format <- "long"
results_df_l$one_reg <- FALSE

all <- rbind(all,results_df_l,fill=T)
all

write.csv(all,"bench.csv",row.names = F)
gc()

## Long and 1 regressions ----
n_it <- 3
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw","lpcsa","lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking méthode:", method, "\n")
  results[[method]] <- bench_method(
    df = df_l,
    method = method,
    par_type_horizon = "long",
    par_one_reg = TRUE,
    min_iterations = n_it
  )
}

# Display results
results_df_l1 <- rbindlist(results, idcol = "method")
print(results_df_l1[, .(method, median_time,  iterations)])

results_df_l1$format <- "long"
results_df_l1$one_reg <- TRUE

all <- rbind(all,results_df_l1,fill=T)
all

write.csv(all,"bench.csv",row.names = F)


# Large database ----
rm(df_l)
rm(toydata)
gc()
all <- NULL

# Creation of a toy data 
#Number of individuals
n_indiv <- 200000L


#Number of firms
n_firms <- 10000L

toydata <- sim_staggered_panel(n=n_indiv,n_firms=n_firms,
                               t_min_treat = 5L,
                               t_max_treat = 15L)

format(object.size(toydata),units="Mb")
toydata <- toydata[,list(id,t,firm_id,female,treat,log_earnings,g_i)]
toydata[,g_i:=fifelse(is.na(g_i),0,g_i)]

n_it <- 5

## Format wide -------------------------------
par_type_horizon <- "wide"
par_one_reg <- F



## Wide horizon database ---------------------
result <- bench::mark({df  <- lpdidcsa_data(toydata, 
                                            unit = "id", 
                                            time = "t",
                                            dependent = "log_earnings",
                                            treat = "treat",
                                            type_horizon= par_type_horizon
)},
iterations = n_it,
check = FALSE  
)

size_df <- format(object.size(df),units="Mb")
size_df
colnames(df)
nrow(df)

speed <- data.table(
  method = "data wide",
  median_time = result$median,  
  min_time = min(unlist(result$time)),
  max_time = max(unlist(result$time)),
  iterations = result$n_itr,
  n_rows=nrow(df),
  n_cols=ncol(df),
  size=size_df)
speed

all <- rbind(all,speed,fill=T)
write.csv(all,"bench_big.csv",row.names = F)

## Wide ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df <- rbindlist(results, idcol = "method")

all <- rbind(all,results_df,fill=T)
all
write.csv(all,"bench_big.csv",row.names = F)
gc()


### Wide interacted ----
methods <- c("lpdid_ipw","lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",  # ou "long"
    control="female:i(t)",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

results_df_int <- rbindlist(results, idcol = "method")
print(results_df_int[, .(method, median_time,  iterations)])

all <- rbind(all,results_df_int,fill=T)
all

write.csv(all,"bench_big.csv",row.names = F)

n_it <- 2 

## CSA ipw ----
time_csa_ipw <- bench::mark({csa_ipw_det <- att_gt(yname="log_earnings",
                                                   tname="t",
                                                   idname="id",
                                                   gname="g_i",
                                                   base_period="universal",
                                                   control_group = c("notyettreated"),
                                                   est_method="ipw",
                                                   bstrap = FALSE,
                                                   cband = FALSE,
                                                   xformla=~female,
                                                   data=toydata)
csa_ipw <- aggte(csa_ipw_det, type = "dynamic")},
iterations = n_it,
check = FALSE,
memory = FALSE
)    
rm(csa_ipw,csa_ipw_det)

time_csa_ipw

speed <- data.table(
  method = "csa ipw",
  median_time = time_csa_ipw$median,  
  min_time = min(unlist(time_csa_ipw$time)),
  max_time = max(unlist(time_csa_ipw$time)),
  iterations = time_csa_ipw$n_itr)
speed

all <- rbind(all,speed,fill=T)
all
write.csv(all,"bench_big.csv",row.names = F)

## CSA reg ----
time_csa_reg <- bench::mark({csa_reg_det <- att_gt(yname="log_earnings",
                                                   tname="t",
                                                   idname="id",
                                                   gname="g_i",
                                                   base_period="universal",
                                                   control_group = c("notyettreated"),
                                                   est_method="reg",
                                                   bstrap = FALSE,
                                                   cband = FALSE,
                                                   xformla=~female,
                                                   data=toydata)
csa_reg <- aggte(csa_reg_det, type = "dynamic")},
iterations = n_it,
check = FALSE,
memory = FALSE)    
rm(csa_reg,csa_reg_det)

speed <- data.table(
  method = "csa reg",
  median_time = time_csa_reg$median,  
  min_time = min(unlist(time_csa_reg$time)),
  max_time = max(unlist(time_csa_reg$time)),
  iterations = time_csa_reg$n_itr)
speed

all <- rbind(all,speed,fill=T)
all
write.csv(all,"bench_big.csv",row.names = F)
gc()


# Very large database ----
rm(df)
rm(toydata)
gc()
all <- NULL

# Creation of a toy data 
#Number of individuals
n_indiv <- 2000000L

#Number of firms
n_firms <- 100000L


toydata <- sim_staggered_panel(n=n_indiv,n_firms=n_firms,
                               t_min_treat = 5L,
                               t_max_treat = 15L)

format(object.size(toydata),units="Mb")
toydata <- toydata[,list(id,t,firm_id,female,treat,log_earnings,g_i)]
toydata[,g_i:=fifelse(is.na(g_i),0,g_i)]

n_it <- 1

## Format wide -------------------------------
par_type_horizon <- "wide"
par_one_reg <- F



## Wide horizon database ---------------------
result <- bench::mark({df  <- lpdidcsa_data(toydata, 
                                            unit = "id", 
                                            time = "t",
                                            dependent = "log_earnings",
                                            treat = "treat",
                                            type_horizon= par_type_horizon
)},
iterations = n_it,
check = FALSE  
)

size_df <- format(object.size(df),units="Mb")
size_df
colnames(df)
nrow(df)

speed <- data.table(
  method = "data wide",
  median_time = result$median,  
  min_time = min(unlist(result$time)),
  max_time = max(unlist(result$time)),
  iterations = result$n_itr,
  n_rows=nrow(df),
  n_cols=ncol(df),
  size=size_df)
speed

all <- rbind(all,speed,fill=T)
write.csv(all,"bench_very_big.csv",row.names = F)

## Wide ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df <- rbindlist(results, idcol = "method")

all <- rbind(all,results_df,fill=T)
all
write.csv(all,"bench_very_big.csv",row.names = F)
gc()


## Wide lpcsa ----
methods <- c("lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",
    par_one_reg = FALSE,
    min_iterations = n_it
  )
}

# Display results
results_df <- rbindlist(results, idcol = "method")

all <- rbind(all,results_df,fill=T)
all
write.csv(all,"bench_very_big.csv",row.names = F)
gc()

rm(list=ls())
gc()
