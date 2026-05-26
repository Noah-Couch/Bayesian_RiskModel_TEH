
########################## Alternative Methods to TEH ##########################

Alternative_Methods <- function(Trial_Data, Noise.Vars, Het_Type) {
  
  ### Calculating Interaction p-value ------------------------------------------
  
  noise_str <- paste(grep("^X", names(Trial_Data), value = TRUE), collapse = " + ")
  formula_str <- paste("Outcome ~ Treatment * Z4 + Z1 + Z2 + Z3 +", noise_str)
  
  
  mod.PgV <- lm(as.formula(formula_str),
                data = Trial_Data)
  
  pval.PgV <- summary(mod.PgV)$coefficients[dim(summary(mod.PgV)$coefficients)[1],4]
  
  if(Het_Type == "linComb") {
    formula_str2 <- paste("Outcome ~ Treatment * Z2 + Z1 + Z3 + Z4 +", noise_str)
    
    
    mod.PgV2 <- lm(as.formula(formula_str2),
                  data = Trial_Data)
    
    pval.PgV2 <- summary(mod.PgV)$coefficients[7+Noise.Vars,4]
    
    pval.PgV <- min(pval.PgV, pval.PgV2)
  }
  
  
  
  ### Calculating SDIR p-value -------------------------------------------------
  
  ### SDIR
  var.Trt  <- var(filter(Trial_Data, Treatment == 1)$Outcome) 
  var.Ctrl <- var(filter(Trial_Data, Treatment == 0)$Outcome) 
  sdir <- sqrt( abs( var.Trt - var.Ctrl ) )
  
  ### P-value from F-test
  sdir.pval <- var.test(filter(Trial_Data, Treatment == 1)$Outcome,
                        filter(Trial_Data, Treatment == 0)$Outcome,
                        ratio = 1,
                        alternative = "greater",
                        conf.level  = 0.95)$p.value
  
  return(list(Intx_pval = pval.PgV,
              SDIR = sdir,
              SDIR_pval = sdir.pval))
}
