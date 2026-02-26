
########################## Alternative Methods to TEH ##########################

Alternative_Methods <- function(Trial_Data) {
  
  ### Calculating Interaction p-value ------------------------------------------
  
  mod.PgV <- lm(Outcome ~ Treatment * ProgVar + Age + Gender + Outcome_BL,
                data = Trial_Data)
  
  pval.PgV <- summary(mod.PgV)$coefficients[7,4]
  
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
