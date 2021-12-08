install.packages("anomaly")
library("anomaly")
data <- read.csv("SP500.csv")

prices <- data$S.P.500
dates <- data$Effective.date

X = diff( log(prices) )

out_capa <- capa( X, beta_tilde = 20, beta=120 )

plot(out_capa)

out_scapa <- scapa.uv( X, beta_tilde = 20, beta=120 ) 

plot( out_scapa )