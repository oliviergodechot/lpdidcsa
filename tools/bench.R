# ______________________________________________________________________________
# Lpdidcsa timing Script 
# ______________________________________________________________________________


library(bench)
library(data.table)
library(lpdidcsa)
library(did)
library(fixest)

# Parameters ----


# List of  methods to test
# Beware lpdid_adj is very slow. Consider dropping it, if you test this script on a large dataframe.

methods <- c("lpdid", "lpdid_rw", "lpdid_adj","lpdid_ipw", "lpcsa", "lpcsa_ipw")
# methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")

#Number of individuals
no_indiv <- 5000L

#Number of firms
no_firms <- 100L

# Creation of a toy data 
toydata <- sim_staggered_panel(n=no_indiv,n_firms=no_firms,
                               t_min_treat = 2L,
                               t_max_treat = 19L)

format(object.size(toydata),units="Mb")
toydata <- toydata[,list(id,t,firm_id,female,treat,log_earnings,g_i)]
toydata[,g_i:=fifelse(is.na(g_i),0,g_i)]



# Format wide -------------------------------
par_type_horizon <- "wide"
par_one_reg <- F

min_iterations <- 2

## Wide horizon database ---------------------
time <- bench::mark({df  <- lpdidcsa_data(toydata, 
                                          unit = "id", 
                                          time = "t",
                                          dependent = "log_earnings",
                                          treat = "treat",
                                          type_horizon= par_type_horizon
)},
iterations = min_iterations,
check = FALSE  
)
df[,g_i:=NULL]
time$median


size_df <- format(object.size(df),units="Mb")
size_df
colnames(df)

## Long horizon database ---------------------
time <- bench::mark({df_l  <- lpdidcsa_data(toydata, 
                                          unit = "id", 
                                          time = "t",
                                          dependent = "log_earnings",
                                          treat = "treat",
                                          type_horizon= "long"
)},
iterations = min_iterations,
check = FALSE  # Désactive la vérification des résultats identiques
)
df_l[,g_i:=NULL]
time$median


size_df <- format(object.size(df_l),units="Mb")
size_df
colnames(df_l)



nb_it <- 10

method <- "lpdid_ipw"

# Function for benchmarking the methods with  bench::mark & tryCatch
bench_method <- function(df, method, par_type_horizon, par_one_reg, min_iterations = 3) {
  tryCatch({
    result <- bench::mark(
      {
        res <- lpdidcsa(
          data=df, 
          unit = "id", time = "t",
          dependent = "log_earnings",
          meth = method,
          controls = "female",
          type_horizon = par_type_horizon,
          one_reg = par_one_reg
        )
      },
      iterations = min_iterations,
      check = FALSE  
    )
    
    # Simplified summary
    list(
      method = method,
      median_time = result$median,  
      memory = result$mem_alloc,  
      iterations = result$n_itr   
    )
  }, error = function(e) {
    message(paste("Error for the method", method, ":", e$message))
    return(NULL)
  })
}


# Wide ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",
    par_one_reg = FALSE,
    min_iterations = nb_it
  )
}

# Display results
results_df <- rbindlist(results, idcol = "method")
print(results_df[, .(method, median_time, memory, iterations)])


# Wide method adj----
methods <- c("lpdid_adj")
results <- list()
for (method in methods) {
  cat("Benchmarking method:", method, "\n")
  results[[method]] <- bench_method(
    df = df,
    method = method,
    par_type_horizon = "wide",  # ou "long"
    par_one_reg = FALSE,
    min_iterations = 2
  )
}

# Display results
results_df_adj <- rbindlist(results, idcol = "method")
print(results_df_adj[, .(method, median_time, memory, iterations)])


# CSA ipw ----
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
                    iterations = nb_it,
                    check = FALSE  )    

time_csa_ipw


# CSA reg ----
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
iterations = nb_it,
check = FALSE  )    

time_csa_reg




# Long multiple regressions ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw", "lpcsa", "lpcsa_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking méthode:", method, "\n")
  results[[method]] <- bench_method(
    df = df_l,
    method = method,
    par_type_horizon = "long",  # ou "long"
    par_one_reg = FALSE,
    min_iterations = nb_it
  )
}

# Display results
results_df_l <- rbindlist(results, idcol = "method")
print(results_df_l[, .(method, median_time, memory, iterations)])


nb_it <- 3

# Long and 1 regressions ----
methods <- c("lpdid", "lpdid_rw", "lpdid_ipw")
results <- list()
for (method in methods) {
  cat("Benchmarking méthode:", method, "\n")
  results[[method]] <- bench_method(
    df = df_l,
    method = method,
    par_type_horizon = "long",
    par_one_reg = TRUE,
    min_iterations = nb_it
  )
}

# Display results
results_df_l1 <- rbindlist(results, idcol = "method")
print(results_df_l1[, .(method, median_time, memory, iterations)])
str(results_df)


# Aggregation of results -----

# did csa results
csa <- rbind( setDT(time_csa_ipw[,c("median","mem_alloc","n_itr")]),
              setDT(time_csa_reg[,c("median","mem_alloc","n_itr")]))
colnames(csa) <- c("median_time","memory","iterations")
csa
csa$method <- c("did:ipw","did:reg")
csa
csa$format <- "wide"

results_df$format <- "wide"
results_df$one_reg <- FALSE

results_df_l$format <- "long"
results_df_l$one_reg <- FALSE

results_df_l1$format <- "long"
results_df_l1$one_reg <- TRUE

all <- rbind(setDT(results_df),setDT(results_df_l),setDT(results_df_l1),csa,fill=T)
all <- all[,-2]

all

write.csv(all,"bench.csv",row.names = F)
