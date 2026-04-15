
########################## Alternative Methods to TEH ##########################

Alternative_Methods <- function(Trial_Data, Noise.Vars) {
  
  ### Calculating Interaction p-value ------------------------------------------
  
  noise_str <- paste(grep("^X", names(Trial_Data), value = TRUE), collapse = " + ")
  formula_str <- paste("Outcome ~ Treatment * Z4 + Z1 + Z2 + Z3 +", noise_str)
  
  
  mod.PgV <- lm(as.formula(formula_str),
                data = Trial_Data)
  
  pval.PgV <- summary(mod.PgV)$coefficients[7+Noise.Vars,4]
  
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
