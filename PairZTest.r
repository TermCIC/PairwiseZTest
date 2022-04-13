# Project: Evaluation on ten rearing materials of BSF
# based on growth limit, growth rate, and the delayed time of pupation
# Function: Pairwise Z tests on coefficnets of logistic growth functions
# Author: Chun-I Chiu; tim77918@msn.com


# ======================================================================#
# The function of Z test, get p-value from normal distribution
# Four inputs
# B1 is the coefficient, and SE1 are coefficient  from model 1
# B2 and SE2 are of coefficient from model 2
z.test <- function(B1, B2, SE1, SE2) {
  distance <- (B1 - B2)
  z.score <- distance / (SE1 ^ 2 + SE2 ^ 2) ^ 0.5
  
  if (z.score > 0) {
    p.value <- pnorm(1 - z.score)
  }
  else{
    p.value <- pnorm(z.score)
  }
  
  return(c(distance, z.score, p.value))
}

# Pairwise comparison & grouping
require(agricolae)
pairwise.z.test <-
  function(label, means, SE) {
    Npair <- length(label)
    try(if (length(means) != Npair |
            length(SE) != Npair)
      stop('Unequal length among inputs!'))
    frame <- matrix(1, nrow = Npair, ncol = Npair)
    rownames(frame) <- colnames(frame) <- label
    for (i in 1:Npair) {
      for (j in 1:Npair) {
        getValue <- z.test(means[j], means[i], SE[j], SE[i])[3] * (Npair - 1)
        frame[i, j] <- frame[j, i] <- as.numeric(getValue)
      }
    }

    group <-
      orderPvalue(
        treatment = label,
        means = means,
        alpha = 0.05,
        pvalue = frame,
        console = TRUE
      )
    return(group)
  }

# ======================================================================#
# Now input your data here
Labels = c(
  'Kitchen leftover',
  'Chicken manure',
  'Pineapple',
  'Soybean',
  "Wet brewer's grain",
  'Rice',
  'Chicken meat',
  'Rice + Soybean',
  'Rice + Chicken meat',
  "Rice + Wet brewer's grain"
)
L = c(5.905,
      4.3583,
      5.9097,
      5.5929,
      5.1468,
      5.7406,
      5.3807,
      5.3491,
      5.516,
      5.9423)
SEL = c(0.0713,
        0.0727,
        0.1884,
        0.1092,
        0.0506,
        0.0478,
        0.0753,
        0.0522,
        0.1025,
        0.0854)
k = c(0.4199,
      0.2627,
      0.1292,
      0.5037,
      0.6279,
      0.3126,
      0.4103,
      0.6091,
      0.4925,
      0.5286
)
SEk = c(0.0242,
        0.018,
        0.0062,
        0.0373,
        0.0316,
        0.0132,
        0.0286,
        0.0365,
        0.0362,
        0.0298
)
c = c(4.242,
      3.1745,
      2.6208,
      3.4158,
      4.787,
      3.5735,
      3.3304,
      3.9042,
      3.0618,
      4.5598
      
)
SEc = c(0.2356,
        0.2044,
        0.0714,
        0.2333,
        0.2343,
        0.1474,
        0.2264,
        0.2333,
        0.2091,
        0.2411
)
tl = c(20.20481067,
       24.16825276,
       40.56965944,
       13.56283502,
       15.2476509,
       22.86308381,
       16.23397514,
       12.81956986,
       12.43370558,
       17.25236474
)
SEtl = c(1.617168855,
         2.272421229,
         2.238704418,
         1.366325233,
         1.070421361,
         1.34959605,
         1.580630419,
         1.084882699,
         1.247502806,
         1.333459022
)
td = c(6.795189331,
       9.83174724,
       -4.569659443,
       3.437164979,
       4.7523491,
       13.13691619,
       9.76602486,
       10.18043014,
       4.566294416,
       1.747635263
)
SEtd = c(1.617168855,
         2.272421229,
         2.238704418,
         1.366325233,
         1.070421361,
         1.34959605,
         1.580630419,
         1.084882699,
         1.247502806,
         1.333459022
)

pairwise.z.test(Labels, L, SEL)
pairwise.z.test(Labels, k, SEk)
pairwise.z.test(Labels, c, SEc)
pairwise.z.test(Labels, tl, SEtl)
pairwise.z.test(Labels, td, SEtd)
