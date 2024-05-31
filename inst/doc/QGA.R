## ----setup, include = FALSE---------------------------------------------------
library(QGA)

## -----------------------------------------------------------------------------
KnapsackProblem <- function(solution,eval_func_inputs) {
  solution <- solution - 1
  items <- eval_func_inputs[[1]]
  maxweight <- eval_func_inputs[[2]]
  # Fitness
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(tot_items)
}

## -----------------------------------------------------------------------------
items <- as.data.frame(list(Item = paste0("item",c(1:500)),
                            weight = rep(NA,500)))
set.seed(1234)
items$weight <- rnorm(500,mean=200,sd=80)
head(items)


## -----------------------------------------------------------------------------
sum(items$weight)
maxweight = sum(items$weight) / 5
maxweight

## -----------------------------------------------------------------------------
popsize = 20
generation_max = 500
nvalues_sol = 2
Genome = nrow(items)
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.025
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome+1)
mutation_rate_end = 2/(Genome+1)
mutation_flag = TRUE
plotting = FALSE
verbose = FALSE
progress = FALSE
eval_fitness = KnapsackProblem
eval_func_inputs = list(items,maxweight)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  knapsackSolution  <- QGA(popsize,
#                  generation_max,
#                  nvalues_sol,
#                  Genome,
#                  thetainit,
#                  thetaend,
#                  pop_mutation_rate_init,
#                  pop_mutation_rate_end,
#                  mutation_rate_init,
#                  mutation_rate_end,
#                  mutation_flag,
#                  plotting,
#                  verbose,
#                  progress,
#                  eval_fitness,
#                  eval_func_inputs)

## ----eval=TRUE, echo=FALSE, include=FALSE-------------------------------------
load("knapsackSolution.RData")

## -----------------------------------------------------------------------------
QGA:::plot_Output(knapsackSolution [[2]])

## -----------------------------------------------------------------------------
best <- knapsackSolution[[1]] - 1
sum(best)

## -----------------------------------------------------------------------------
sum(items$weight[best])
maxweight

## -----------------------------------------------------------------------------
TravellingSalesman <- function(solution,distance) {
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  # Fitness function
  l = l + distance[solution[1],solution[length(solution)]]
  # Penalization
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- -(penal+l)
  return(cost)
}

## -----------------------------------------------------------------------------
cities <- read.csv("cities.csv")
ncities <- 9
cities <- cities[c(1:ncities),]
cities

## -----------------------------------------------------------------------------
distance <- as.matrix(dist(cities[,c(2:3)]))
distance

## ----eval=FALSE---------------------------------------------------------------
#  popsize = 20
#  Genome = nrow(cities)
#  nvalues_sol = nrow(cities)
#  set.seed(4321)
#  TSPsolution <- QGA(popsize,
#                  generation_max = 1000,
#                  nvalues_sol,
#                  Genome,
#                  thetainit = 3.1415926535 * 0.01,
#                  thetaend = 3.1415926535 * 0.01,
#                  # pop_mutation_rate_init = 1/(popsize + 1),
#                  # pop_mutation_rate_end = 1/(popsize + 1),
#                  # mutation_rate_init = 1/(Genome + 1),
#                  # mutation_rate_end = 1/(Genome + 1),
#                  mutation_flag = FALSE,
#                  plotting = FALSE,
#                  verbose = FALSE,
#                  progress = FALSE,
#                  eval_fitness = TravellingSalesman,
#                  eval_func_inputs = distance)
#  

## ----eval=TRUE, echo=FALSE, include=FALSE-------------------------------------
load("TSPsolution.RData")

## -----------------------------------------------------------------------------
QGA:::plot_Output(TSPsolution[[2]])

## -----------------------------------------------------------------------------
solution <- TSPsolution[[1]]
cities$city[solution]
cities_tsp <- cities[solution,]
plot(y~x,data=cities_tsp)
polygon(cities_tsp$x,cities_tsp$y,border="red")
text(x = cities_tsp$x, y = cities_tsp$y, labels = cities_tsp$city, cex=.75)
title("Best path")

## -----------------------------------------------------------------------------
clustering <- function(solution, eval_func_inputs) {
  maxvalue <- 5
  penalfactor <- 2
  df <- eval_func_inputs[[1]]
  vars <- eval_func_inputs[[2]]
  # Fitness function
  fitness <- 0
  for (v in vars) {
    cv <- tapply(df[,v],solution,FUN=sd) / tapply(df[,v],solution,FUN=mean)
    cv <- ifelse(is.na(cv),maxvalue,cv)
    fitness <- fitness + sum(cv)
  }
  # Penalization on unbalanced clusters
  b <- table(solution)/nrow(df)
  fitness <- fitness + penalfactor * (sum(abs(b - c(rep(1/(length(b)),length(b))))))
  return(-fitness)
}

## -----------------------------------------------------------------------------
data(iris)
vars <- colnames(iris)[1:4]
vars

## ----eval=FALSE---------------------------------------------------------------
#  nclust = 3
#  popsize = 20
#  Genome = nrow(iris)
#  set.seed(1234)
#  solutionQGA <- QGA(popsize,
#                  generation_max = 1500,
#                  nvalues_sol = nclust,
#                  Genome,
#                  thetainit = 3.1415926535 * 0.1,
#                  thetaend = 3.1415926535 * 0.05,
#                  pop_mutation_rate_init = 1/(popsize + 1),
#                  pop_mutation_rate_end = 1/(popsize + 1),
#                  mutation_rate_init = 1/(Genome + 1),
#                  mutation_rate_end = 1/(Genome + 1),
#                  mutation_flag = TRUE,
#                  plotting = FALSE,
#                  verbose = FALSE,
#                  progress = FALSE,
#                  eval_fitness = clustering,
#                  eval_func_inputs = list(iris, vars))

## ----eval=TRUE, echo=FALSE, include=FALSE-------------------------------------
load("CLUSTsolution.RData")

## -----------------------------------------------------------------------------
QGA:::plot_Output(CLUSTsolution[[2]])

## -----------------------------------------------------------------------------
solution <- CLUSTsolution[[1]]
table(solution)

## -----------------------------------------------------------------------------
iris$cluster <- solution
xtabs( ~ Species + cluster, data=iris)

