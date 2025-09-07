library(readxl)
library(dplyr)
library(ggplot2)
library(stats)

# 1. Read Excel file
df <- read.csv('temp/cleaned_data.csv')

# 2. List columns with missing value ratio > 0.1
missing_ratio <- sapply(df, function(x) mean(is.na(x)))
cols_complete <- names(missing_ratio[missing_ratio == 0])
cols_near_complete <- names(missing_ratio[(missing_ratio > 0)&(missing_ratio <0.1)])
cols_missing_gt_01 <- names(missing_ratio[missing_ratio > 0.1])
cat("\nColumns with no missing values:\n", cols_complete, "\n")
cat("\nColumns with missing ratio between 0 and 0.1:\n", cols_near_complete, "\n")
cat("\nColumns with missing ratio > 0.1:\n", cols_missing_gt_01, "\n")

df$COMTIuse <- !is.na(df$`COMT-Istart`)
cols = c(
    "ID", 
    "Sex",
    "HY", 
    "COMTIuse",
    "AOE", 
    "LD", 
    "LD_LEDD",
    "Total_LEDD",
    "LDBW",
    "Dopa_LEDD",
    "DCI", 
    "AST", 
    "ALT", 
    "LDH", 
    "Cre", 
    "Hb", 
    "VitB6_Pyridoxal",
    "logB6",
    "VitB12",
    "FA",
    "BW", 
    "DD", 
    "CK",
    "Alb",
    "Hcy",
    "gGTP",
    "Alb",
    "ChE",
    "VitB1"
)
cat("\nSelected columns:\n", cols, "\n")

df <- df %>% select(all_of(cols))

# Create a histgram of numeric variables other than ID and save it to report/histograms.pdf
numeric_cols <- sapply(df, is.numeric)
cols_numeric <- names(numeric_cols[numeric_cols])
cols_numeric <- setdiff(cols_numeric, "ID")
pdf("report/histograms.pdf")
for (col in cols_numeric) {
    p <- ggplot(df, aes_string(x = col)) +
        geom_histogram(bins = 30, na.rm = TRUE) +
        ggtitle(paste('Histogram of', col)) +
        xlab(col) +
        ylab('Frequency')
    print(p)
}
dev.off() 

# filter outliers (> 4 SD from mean) for numeric columns
df_no_outliers <- df
for (col in cols_numeric) {
  mean_col <- mean(df[[col]], na.rm = TRUE)
  sd_col   <- sd(df[[col]], na.rm = TRUE)
  high_thresh <- mean_col + 4 * sd_col
  low_thresh <- mean_col - 4 * sd_col
  # count the number of outliers
  outlier_idx <- (df[[col]] < low_thresh) | (df[[col]] > high_thresh)
  num_outliers <- sum(outlier_idx, na.rm = TRUE)
  if (num_outliers > 0) {
      cat(paste("Column:", col, num_outliers, 'outliers', "High Threshold:", round(high_thresh, 2), "Low Threshold:", round(low_thresh, 2), "\n"))
  }
  df_no_outliers[[col]] <- ifelse(outlier_idx, NA, df_no_outliers[[col]])
}

# draw histogram again after outlier removal
pdf("report/histograms_4sd_removed.pdf")
for (col in cols_numeric) {
    p <- ggplot(df_no_outliers, aes_string(x = col)) +
        geom_histogram(bins = 30, na.rm = TRUE) +
        ggtitle(paste('Histogram of', col, 'after outlier removal')) +
        xlab(col) +
        ylab('Frequency')
    print(p)
}
dev.off()

# Hcys associations
model <- lm(log(Hcy) ~ AOE + I(AOE^2) + DD + BW + DCI + Sex + COMTIuse + VitB6_Pyridoxal + VitB12 + FA + LDBW + I(as.factor(HY)), data = df_no_outliers)
summary(model)

model <- lm(Hcy ~ AOE + DD + DCI + Sex + BW + COMTIuse + VitB6_Pyridoxal + VitB12 + FA + LD + I(as.factor(HY)), data = df_no_outliers)
summary(model)


# scatter plot of Hcy vs FA
plot(Hcy~LD, data=df_no_outliers)

# stratified analysis by LD
range = 150
ld_centers <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)
results <- data.frame(ld_center = numeric(), 
        N = numeric(), N_user = numeric(), N_non_user = numeric(),
        estimate = numeric(), sd = numeric(), p = numeric())

for (ld_center in ld_centers) {
    print(paste("Analyzing LD center:", ld_center))
    data_range = df_no_outliers[(df_no_outliers$LD>(ld_center-range))&(df_no_outliers$LD<(ld_center+range)),] 
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
# model2 <- lm(logB6 ~ AOE + DD + BW + Sex + I(as.factor(HY)) + LD_LEDD + COMTIuse, data = df_no_outliers)
# summary(model2)


# Local estimate for each LD value (±range), stratified by COMTIuse
range <- 150
ld_seq <- seq(min(400, na.rm=TRUE), 1000, length.out=100)
pred_list <- list()

for (ld_center in ld_seq) {
  df_sub <- df_no_outliers %>%
    filter(LD > (ld_center - range), LD < (ld_center + range))
  if (nrow(df_sub) < 10) next  # skip if too few samples

  # Set covariates to mean/mode in local window
  mean_vals <- df_no_outliers %>%
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

