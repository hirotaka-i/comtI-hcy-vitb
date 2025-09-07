library(readxl)
library(dplyr)
library(ggplot2)
library(stats)

# 1. Read Excel file
df <- read_excel('data/Vitamine B6 data.xlsx', na = c('NA', 'NaN', ' ', '', 'n.a.'))

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
df$logHcy = log(df$Hcy)
df$logFA = log(df$FA)
df$logVitB12 = log(df$VitB12)
df$logVitB6 = log(df$VitB6_Pyridoxal)

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

# save the cleaned data
write.csv(df_no_outliers, 'temp/cleaned_data.csv', row.names = FALSE)