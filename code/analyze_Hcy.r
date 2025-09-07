library(readxl)
library(dplyr)
library(ggplot2)
library(stats)

# 1. Read Excel file
df <- read.csv('temp/cleaned_data.csv')

# Baseline Hcys associations
model <- lm(Hcy ~ AOE + DD + DCI + Sex + BW + COMTIuse + VitB6_Pyridoxal + VitB12 + FA + LD + I(as.factor(HY)), data = df)
summary(model)

log_model <- lm(logHcy ~ AOE + I(AOE^2) + DD + BW + DCI + Sex + COMTIuse + logVitB6 + logVitB12 + logFA + LD + I(as.factor(HY)), data = df)
summary(log_model)

# stratified analysis by LD
range = 150
ld_centers <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)
results <- data.frame(ld_center = numeric(), 
        N = numeric(), N_user = numeric(), N_non_user = numeric(),
        estimate = numeric(), sd = numeric(), p = numeric())

for (ld_center in ld_centers) {
    print(paste("Analyzing LD center:", ld_center))
    data_range = df[(df$LD>(ld_center-range))&(df$LD<(ld_center+range)),] 
    model_range <- lm(Hcy ~ AOE + DD + DCI + Sex + BW + COMTIuse + VitB6_Pyridoxal + VitB12 + FA + LD + I(as.factor(HY)), 
        data = data_range)
    k = summary(model_range)
    n = nrow(model_range$model)
    n_user = sum(model_range$model$COMTIuse, na.rm=TRUE)
    n_non_user = n - n_user
    kc = k$coefficients
    est = kc['COMTIuseTRUE','Estimate']
    sd = kc['COMTIuseTRUE','Std. Error']
    p = kc['COMTIuseTRUE','Pr(>|t|)']        
    results <- rbind(results, data.frame(ld_center = ld_center, N = n, N_user = n_user, N_non_user = n_non_user, estimate = est, sd = sd, p = p))
}
results
# Plot ld_center vs est with confidence interval ribbon
results <- results %>%
    mutate(ci_lower = estimate - 1.96 * sd,
           ci_upper = estimate + 1.96 * sd)

ggplot(results, aes(x = ld_center, y = est)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "skyblue", alpha = 0.4) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    xlab("LD_LEDD center") +
    ylab("Estimate (COMTIuse effect)") +
    ggtitle("Estimate of COMTIuse effect by LD_LEDD center\nwith 95% Confidence Interval") +
    theme_minimal()

# # VitB6_Pyridoxal associations
# model2 <- lm(logB6 ~ AOE + DD + BW + Sex + I(as.factor(HY)) + LD_LEDD + COMTIuse, data = dfers)
# summary(model2)


# Local estimate for each LD value (±range), stratified by COMTIuse
range <- 150
ld_seq <- seq(min(400, na.rm=TRUE), 1000, length.out=100)
pred_list <- list()

for (ld_center in ld_seq) {
  df_sub <- df %>%
    filter(LD > (ld_center - range), LD < (ld_center + range))
  if (nrow(df_sub) < 10) next  # skip if too few samples

  # Set covariates to mean/mode in local window
  mean_vals <- df %>%
    summarise(
      AOE = mean(AOE, na.rm=TRUE),
      DD = mean(DD, na.rm=TRUE),
      BW = mean(BW, na.rm=TRUE),
      DCI = mean(DCI, na.rm=TRUE),
      Sex = mean(Sex, na.rm=TRUE),
      VitB6_Pyridoxal = mean(VitB6_Pyridoxal, na.rm=TRUE),
      VitB12 = mean(VitB12, na.rm=TRUE),
      FA = mean(FA, na.rm=TRUE),
      HY = names(sort(table(HY), decreasing=TRUE))[1]
    )

  for (comti in c(TRUE, FALSE)) {
    newdata <- data.frame(
      LD = ld_center,
      COMTIuse = comti,
      AOE = mean_vals$AOE,
      DD = mean_vals$DD,
      BW = mean_vals$BW,
      DCI = mean_vals$DCI,
      Sex = mean_vals$Sex,
      VitB6_Pyridoxal = mean_vals$VitB6_Pyridoxal,
      VitB12 = mean_vals$VitB12,
      FA = mean_vals$FA,
      HY = as.factor(mean_vals$HY)
    )
    model_local <- lm(Hcy ~ AOE + DD + DCI + Sex + BW + COMTIuse + VitB6_Pyridoxal + VitB12 + FA + LD + I(as.factor(HY)), data = df_sub)
    pred <- predict(model_local, newdata=newdata, se.fit=TRUE)
    ci_lower <- pred$fit - 1.96 * pred$se.fit
    ci_upper <- pred$fit + 1.96 * pred$se.fit
    pred_list[[length(pred_list)+1]] <- data.frame(
      LD=ld_center, COMTIuse=comti, Hcy_pred=pred$fit, ci_lower=ci_lower, ci_upper=ci_upper
    )
  }
}

pred_df <- do.call(rbind, pred_list)

ggplot(pred_df, aes(x=LD, y=Hcy_pred, color=COMTIuse, fill=COMTIuse)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.2, color=NA) +
  xlab("LD") +
  ylab("Estimated Hcy") +
  ggtitle("Estimated Hcy by LD (local window)\nStratified by COMTIuse with 95% Confidence Interval") +
  theme_minimal()
ggsave("report/ld_hcy_local_estimate.jpg")

