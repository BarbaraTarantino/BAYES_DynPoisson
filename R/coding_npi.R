# Create data for Bayesian models, combining count data and NPI covariates

rm(list = ls(all = TRUE)); graphics.off(); cat("\014")
setwd("/Users/barbaratarantino/Desktop/BAYES_DynPoisson")

library(tidyverse)
library(gtools)
library(ppcor)
library(stringi)
library(corrplot)
library(caret)
library(dplyr)
library(tidyr)

# Data ------------------------------------------------------------------------
{
country = c("France", "Germany", "Italy", "Spain", "United States")
df_final = c()

for (c in country){
  print(c)
if (c == "Italy"){
  Date_initial = "2020-02-27"
  
  data <-read_csv("WHO-COVID-19-global-data.csv") %>%
    filter(Country  == c,
           Date_reported > Date_initial) %>%
    arrange(-desc(Date_reported)) %>%
    rename(Number = New_cases,
           Date = Date_reported,
           location_name = Country,
           location = Country_code) %>%
    mutate(Number = abs(Number)) %>%
    dplyr::select(location_name, location,Date,Number) %>%
    dplyr::mutate(Number = ifelse(Number == 0, NA, Number)) %>%
    tidyr::fill(Number, .direction = c("down"))
  
} else if (c == "France" | c == "Germany" | c == "Spain") {
  Date_initial = "2020-03-05"
  
  data <-read_csv("WHO-COVID-19-global-data.csv") %>%
    filter(Country  == c,
           Date_reported > Date_initial) %>%
    arrange(-desc(Date_reported)) %>%
    rename(Number = New_cases,
           Date = Date_reported,
           location_name = Country,
           location = Country_code) %>%
    mutate(Number = abs(Number)) %>%
    dplyr::select(location_name, location,Date,Number) %>%
    dplyr::mutate(Number = ifelse(Number == 0, NA, Number)) %>%
    tidyr::fill(Number, .direction = c("down"))
  
} else if (c == "United States"){
  Date_initial = "2020-03-08"
  
  data <-read_csv("WHO-COVID-19-global-data.csv") %>%
    filter(Country  == "United States of America",
           Date_reported > Date_initial) %>%
    arrange(-desc(Date_reported)) %>%
    rename(Number = New_cases,
           Date = Date_reported,
           location_name = Country,
           location = Country_code) %>%
    mutate(Number = abs(Number)) %>%
    dplyr::select(location_name, location,Date,Number) %>%
    dplyr::mutate(Number = ifelse(Number == 0, NA, Number)) %>%
    tidyr::fill(Number, .direction = c("down"))
  
}

y <- data %>%
  dplyr::select(Date,Number) %>%
  mutate(Date = as.POSIXct.Date(Date, format = "%Y-%m-%d"))

if (c == "United States"){
  
OxCGRT_latest <- read.csv("OxCGRT_latest.csv")
xreg  <- OxCGRT_latest %>%
  transform(OxCGRT_latest, Date = as.Date(as.character(Date), "%Y%m%d")) %>%
  filter(CountryName  == c,
         Date > Date_initial, 
         Jurisdiction == "NAT_TOTAL") %>%
  filter(Date < tail(data$Date,1)+1 ) %>%
  arrange(-desc(Date)) %>%
  dplyr::select(C1_School.closing, C2_Workplace.closing, 
         C3_Cancel.public.events , C4_Restrictions.on.gatherings, 
         C5_Close.public.transport, C6_Stay.at.home.requirements, 
         C7_Restrictions.on.internal.movement, C8_International.travel.controls, 
         H1_Public.information.campaigns, H2_Testing.policy, H3_Contact.tracing, 
         H4_Emergency.investment.in.healthcare, H5_Investment.in.vaccines, 
         H6_Facial.Coverings, H7_Vaccination.policy, H8_Protection.of.elderly.people) 
} else {
  
  OxCGRT_latest <- read.csv("OxCGRT_latest.csv")
  xreg  <- OxCGRT_latest %>%
    transform(OxCGRT_latest, Date = as.Date(as.character(Date), "%Y%m%d")) %>%
    filter(CountryName  == c,
           Date > Date_initial) %>%
    filter(Date < tail(data$Date,1)+1 ) %>%
    arrange(-desc(Date)) %>%
    dplyr::select(C1_School.closing, C2_Workplace.closing, 
                  C3_Cancel.public.events , C4_Restrictions.on.gatherings, 
                  C5_Close.public.transport, C6_Stay.at.home.requirements, 
                  C7_Restrictions.on.internal.movement, C8_International.travel.controls, 
                  H1_Public.information.campaigns, H2_Testing.policy, H3_Contact.tracing, 
                  H4_Emergency.investment.in.healthcare, H5_Investment.in.vaccines, 
                  H6_Facial.Coverings, H7_Vaccination.policy, H8_Protection.of.elderly.people) 
  
}

# Variable selection
# Binary coding scheme

sum(is.na(xreg))
vars = colnames(xreg)
for (i in vars){
  xreg <- xreg %>%
    fill(all_of(i),.direction= "down")
}
sum(is.na(xreg))

y = data %>%
  dplyr::select(Date,Number)

df_tot_2 <- xreg %>%
  mutate(C1_School.closing = recode(C1_School.closing,
                                    "1" = 0,
                                    "2" = 1,
                                    "3" = 1),
         C2_Workplace.closing = recode(C2_Workplace.closing,
                                       "1" = 0,
                                       "2" = 0,
                                       "3" = 1),
         C3_Cancel.public.events = recode(C3_Cancel.public.events,
                                          "1" = 0,
                                          "2" = 1),
         C4_Restrictions.on.gatherings = recode(C4_Restrictions.on.gatherings,
                                                "1" = 0,
                                                "2" = 0,
                                                "3" = 1, #
                                                "4" = 1),
         C5_Close.public.transport = recode(C5_Close.public.transport,
                                            "1" = 1,
                                            "2" = 1),
         C6_Stay.at.home.requirements = recode(C6_Stay.at.home.requirements,
                                               "1" = 0,
                                               "2" = 1, 
                                               "3" = 1),
         C7_Restrictions.on.internal.movement = recode(C7_Restrictions.on.internal.movement,
                                          "1" = 0,
                                          "2" = 1),
         C8_International.travel.controls = recode(C8_International.travel.controls,
                                         "1" = 0,
                                         "2" = 0, #
                                         "3" = 1,#
                                         "4" = 1),
         H2_Testing.policy = recode(H2_Testing.policy,
                                    "1" = 0,
                                    "2" = 0, 
                                    "3" = 1),
         H3_Contact.tracing = recode(H3_Contact.tracing,
                                     "1" = 0,
                                     "2" = 1),
         H6_Facial.Coverings = recode(H6_Facial.Coverings,
                                      "1" = 0,
                                      "2" = 0, 
                                      "3" = 1,# 
                                      "4" = 1),
          H7_Vaccination.policy = recode( H7_Vaccination.policy,
                                        "1" = 0,
                                        "2" = 1, 
                                        "3" = 1,# 
                                        "4" = 1,
                                        "5" = 1))
df_tot_2 <- df_tot_2 %>%
  dplyr::select(C1_School.closing, C2_Workplace.closing, C3_Cancel.public.events,
                C4_Restrictions.on.gatherings, C5_Close.public.transport,
                C6_Stay.at.home.requirements, C7_Restrictions.on.internal.movement, 
                C8_International.travel.controls, H2_Testing.policy, H3_Contact.tracing, 
                H6_Facial.Coverings,  H7_Vaccination.policy ) %>%
  dplyr::rename(c1 = C1_School.closing, c2 = C2_Workplace.closing, c3 = C3_Cancel.public.events,
                c4 = C4_Restrictions.on.gatherings,  c5 = C5_Close.public.transport,
                c6 = C6_Stay.at.home.requirements, c7 = C7_Restrictions.on.internal.movement,
                c8 = C8_International.travel.controls, h2 = H2_Testing.policy, h3 = H3_Contact.tracing,
                h6 = H6_Facial.Coverings, h7 =  H7_Vaccination.policy )

df_tot_2 <- df_tot_2 %>% 
  mutate(c1_c2 = case_when((c1 == 0 & c2 == 0 ) ~ 0,
                           (c1 == 0 & c2 == 1 |
                              c1 == 1 & c2 == 0 ) ~ 0, 
                           (c1 == 1 & c2 == 1 ) ~ 1),
         c3_c4 = case_when((c3 == 0 & c4 == 0 ) ~ 0,
                           (c3 == 0 & c4 == 1 |
                              c3 == 1 & c4 == 0 ) ~ 0, 
                           (c3 == 1 & c4 == 1 ) ~ 1),
         c5_c6_c7 = case_when((c5 == 0 & c6 == 0 & c7== 0) ~ 0,
                              (c5 == 1 & c6 == 0 & c7 == 0 | 
                                 c5 == 0 & c6 == 1 & c7 == 0 |
                                 c5 == 0 & c6 == 0 & c7 == 1 ) ~ 0,
                              (c5 == 1 & c6 == 1 & c7 == 0  | 
                                 c5 == 0 & c6 == 1 & c7 == 1  | 
                                 c5 == 1 & c6 == 0 & c7 == 1) ~ 1,
                              (c5 == 1 & c6 == 1 & c7 == 1 ) ~ 1),
         h2_h3_h6 = case_when((h2 == 0 & h3 == 0 & h6 == 0) ~ 0,
                              (h2 == 1 & h3 == 0 & h6 == 0 | 
                                 h2 == 0 & h3 == 1 & h6 == 0 |
                                 h2 == 0 & h3 == 0 & h6 == 1 ) ~ 0,
                              (h2 == 1 & h3 == 1 & h6 == 0  | 
                                 h2 == 0 & h3 == 1 & h6 == 1  | 
                                 h2 == 1 & h3 == 0 & h6 == 1) ~ 1,
                              (h2 == 1 & h3 == 1 & h6 == 1 ) ~ 1))
df_tot_2 <- df_tot_2 %>%
  dplyr::select(c1_c2, c3_c4, c5_c6_c7, c8, h2_h3_h6, h7)

df <- y %>%
  cbind(df_tot_2) %>%
  mutate(Country = c) %>%
  dplyr::select(Country,Date,Number, starts_with(c("c", "h")))
  
# data <- df %>%
#   pull(Number) %>% as.vector()
# 
# xreg <- df %>% dplyr::select(starts_with(c("c", "h")))   

df_final <- df_final %>%
  rbind(df)

}

write.csv(df_final, "df_bayes_2_jhu.csv", row.names=F)

}

# NPI exploratory analysis -----------------------------------------------------
{
# Import
df <- read.csv("df_bayes.csv")

c = "Italy"

df = df %>%
  dplyr::filter(Country == c) %>%
  dplyr::select(!Country)

# Contigency table for NPIs' freq 
# Remove constant NPI (frequency one of the two cat < 10)

df_cross <- df %>% 
  dplyr::select(starts_with(c("c", "h"))) %>%
  gather(key,value)
cross <- data.frame(table(df_cross$key,df_cross$value))
cross <- cross %>% 
  spread(Var2, Freq) 

rem_npi <- cross %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(check1 = cross$'0' / cross$'1', 
         check2 = cross$'1' / cross$'0') %>%
  filter_at(vars(check1,check2), any_vars(. < 0.1)) %>%
  dplyr::select(Var1) %>%
  pull() %>% as.vector()

df = df %>% dplyr::select(!all_of(rem_npi))

# Pearson correlation btw count data and NPIs

y = df %>% dplyr::select(Number)
metrics_total = c()
vars_df = df %>%
  dplyr::select(-c(Number, Date))
vars = as.vector(colnames(vars_df))

for (i in vars){
  
  policy_int = df %>%
    dplyr::select(i) 
  
  cor_3 = cor.test(x=unlist(y), y=unlist(policy_int), alternative=NULL, method = "pearson")
  pearson = data.frame(pearson=paste(round(cor_3$estimate,4), stars.pval(cor_3$p.value), sep = " "))
  rownames(pearson) <- i
  
  metrics_total = metrics_total %>%
    rbind(pearson) 
  
}

rem_npi = metrics_total %>% 
  as.data.frame() %>%
  rownames_to_column %>%
  dplyr::rename(NPI = rowname, 
                Pearson = pearson) %>%
  mutate(sig = grepl("\\*$",Pearson)) %>%
  filter(sig == "FALSE") %>%
  pull(NPI) %>% as.vector()

df = df %>% dplyr::select(!all_of(rem_npi))

# Correlation btw policy interventions

df2 <- df %>%
  dplyr::select(-c(Date,Number)) %>%
  filter(complete.cases(.))

# drop_na; na_omit
# df2 <- df2[complete.cases(df2),]

correlations <- cor(df2)
p.correlation <- pcor(df2)
partialcorrelations <- p.correlation$estimate
rownames(partialcorrelations) <- rownames(correlations)
colnames(partialcorrelations) <- colnames(correlations)

# matrix of the p-value of the correlation

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(df2)
p.mat.partial <- p.correlation$p.value

stars = stars.pval(p.mat.partial)
stars = stars %>%
  as.data.frame() 
colnames(stars) <- colnames(partialcorrelations)

partialcorrelations_2 = partialcorrelations %>%
  as.data.frame() %>%
  mutate_all(function(x) round(x,3)) 

var.names = colnames(partialcorrelations_2)
var = vector("list", length(var.names))
for (k in 1:length(var.names)){
  var[[k]]= as.data.frame(stri_join(partialcorrelations_2[,c(var.names[k])],stars[,c(var.names[k])],sep=" "))
  colnames(var[[k]]) <- c(var.names[k])
}


grid_tot <-bind_cols(var) 
suppressWarnings(grid_tot[upper.tri(grid_tot)] <- "")
grid_tot <- grid_tot %>%
  mutate(NPI = colnames(grid_tot)) %>%
  dplyr::select(NPI, starts_with("c")) 

# Highly correlated vars 

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

corr <- flattenCorrMatrix(correlations, p.mat)
corr_tot <- corr %>% 
  subset(cor >= 0.5 & p <= 0.05) 

pcorr <- flattenCorrMatrix(partialcorrelations, p.mat.partial)
pcorr_sig <- pcorr %>% 
  subset(cor >= 0.5 & p <= 0.05)

pcorr_tot <- pcorr %>% 
  subset(cor >= 0.5)

corr_var <- union(unique(pcorr_tot$row), unique(pcorr_tot$column))
df_pcorr <- df2[,corr_var]

p.correlation <- pcor(df_pcorr); partialcorrelations <- p.correlation$estimate
rownames(partialcorrelations) <- colnames(df_pcorr)
colnames(partialcorrelations) <- colnames(df_pcorr)
p.mat.partial <- p.correlation$p.value

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(partialcorrelations, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat.partial, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
highlyCorrelated_part <- findCorrelation(partialcorrelations, cutoff=0.7, names = TRUE)
print(highlyCorrelated_part)

df = df %>%
  dplyr::select(!all_of(highlyCorrelated_part))

}
