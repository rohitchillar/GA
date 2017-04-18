##################################################
###       Genetic Algorithm in Engineering Processes
###       Name:  Rohit Chillar
###       18-4-2017
###       Github link for the complete code is https://github.com/rohitchillar/GA
##################################################

library("datasets")
library("locfit")
library("MASS")
library("BioPhysConnectoR")

#############################
# Primary Function
#############################

# Carries out variable selection optimization using auxiliary functions
# User should enter filename and dependent variable of choice between < >
fileName <- "<filename>" 
dependentVariable <- "<dependent variable>"
dataFile <- readDataFile(fileName, sep="", dependentVariable)

select <- function(initialPopulationSize, numberGenerations, regressionType, criterion="AIC", mutationProbability, replacementPercentage) { 
     parentPopulation <- initiate(dataParameter, initialPopulationSize)
     for(i in 1:numberGenerations) { 
          scoredParent <- runModel(parentPopulation, dependentVariable, regressionType, criterion)
          sortedParent <- mat.sort(scoredParent, ncol(scoredParent), decreasing=FALSE)
          offspringPopulation <- breeding(sortedParent)
          offspringPopulation <- offspringPopulation[,-ncol(offspringPopulation)] 
          mutatedOffspring<- mutation(offspringPopulation)
          scoredOffspring <- runModel(mutatedOffspring, dependentVariable, regressionType, criterion)
          sortedOffspring <- mat.sort(scoredOffspring, ncol(scoredOffspring), decreasing=FALSE)
          newGeneration <- replacement(sortedParent, sortedOffspring, replacementPercentage)
          parentPopulation <- newGeneration[-dim(newGeneration)]
          
     }
     scoredFinalGeneration <- runModel(parentPopulation, dependentVariable, regressionType, criterion)
     sortedFinalGeneration <- mat.sort(scoredFinalGeneration, ncol(scoredFinalGeneration), decreasing = FALSE)
     # Skip over any models that have 0 scores. Due to mutation and crossover, some models may 
     # end up with 0 parameters. This prevents those models from being declared the best.
     index = 1
     while (sortedFinalGeneration[index,dim(sortedFinalGeneration)[2]]==0){
       index = index + 1
     }
     print(sortedFinalGeneration[index,])
}


#############################
# Auxiliary Functions
#############################

#--------------------
# Reading function
#--------------------
# Call before the select function.
# Note that if your data file is NOT comma separated, you will need to change sep to the correct type.

readDataFile <- function(fileName, sep=",", dependentVariable) {
     dataOriginal <- read.table(fileName, sep=sep, header=TRUE, stringsAsFactors=TRUE)
     colnum <- which(colnames(dataOriginal) == dependentVariable)
     dataParameter <- dataOriginal[-colnum]
     dataFileName <- list(dataOriginal, dataParameter)
     return(dataFileName)
}

#---------------------
# Intitiation function
#---------------------
# Called from the primary function. 
# Its purpose is to generate the initial bolean matrix used as parent population for breeding. 
# All zero rows are excluded, since the purpose of the function is for variable selection. 

initiate <- function(data, populationSize){
     columnLength = dim(data)[2]
     #Create a matrix of 0 with desired dimension
     parentPopulation <- matrix(0, nrow=populationSize, ncol=columnLength)
     colnames(parentPopulation) <- c(names(data))
     boolean <- 0:1
     #Using a for loop, fill each row of the matrix with a vector of booleans
     for (i in 1:populationSize) {
          parentPopulation[i,] <- sample(boolean, columnLength, replace=T)
     }
     #Make sure that no all-zero row exist in the matrix
     allZeros <- function(x) {
          y <- rep(0, columnLength)
          isTRUE(all.equal(x,y))
     }
     
     logVec <- sapply(1:populationSize, function(x) allZeros(parentPopulation[x,]))
     while(length(which(logVec == T)) >0) {
          parentPopulation[which(logVec==T),] <- sample(boolean, columnLength, replace=T)
          logVec <- sapply(1:populationSize, function(x) allZeros(parentPopulation[x,]))
     }
     return(as.data.frame(parentPopulation))
}

#--------------------
# Fitness function
#--------------------

# Called from within the runModels function. Determines model fit based on choice of AIC or BIC criteria
# and LM or GLM.
fitness <- function(model, data, dependentVariable, regressionType, criterion="AIC") {
     # Vector to contain the independent variables to be included in the model (excludes dependent).
     includedCols = c()
     # Loop through full dataset and pull out column names NOT equal to dependent variable.
     for (colIndex in 1:length(model)) {
          colName = names(model[colIndex])
          if (model[colIndex] == 1 && colName != dependentVariable) {
               includedCols = c(includedCols, colName)
          }
     }
     if (length(includedCols) == 0) {
          return(0)
     }
     # Create a formula to be used in the LM and GLM functions
     formula <- as.formula(paste(dependentVariable, "~", paste(includedCols, collapse="+")))
     if(regressionType!="LM" && regressionType!="GLM") {
          print("Please use only LM or GLM as regression type")
     } else {
          if(criterion!="AIC" && criterion!="BIC") {
               print("Please use only AIC or BIC as criterion type")
          } else if(criterion=="AIC") {
               if(regressionType=="LM") {
                    lm1 <- lm(formula, data = data)
                    score <- AIC(lm1)
                    return(score)
               } else {
                    glm1 <- glm(formula, data = data)
                    score <- AIC(glm1)
                    return(score)
               }
          } else {
               if(regressionType=="LM") {
                    lm1 <- lm(formula, data = data)
                    score <- BIC(lm1)
                    return(score)
               } else {
                    glm1 <- glm(formula, data = data)
                    score <- BIC(glm1)                    
                    return(score)
               }
          }
     }
     return(0)
}

#Primary function called when calculating fitness of a dataframe of models
runModel <- function(models, dependentVariable, regressionType, criterion="AIC") {
     # Add a column of 0s to the binary models (population) file. 
     models <- cbind(models, 0)
     # Loop through all models (each row) in the population matrix.Populate the newly added column
     # with appropriate fitness scores.
     for (modelIndex in 1:nrow(models)) { 
          model = models[modelIndex, ]
          modelFit <- fitness(model, dataOriginal, dependentVariable, regressionType, criterion)
          model[length(model)] <- modelFit
          models[modelIndex,] = model
     }
     return(models)
}

# Crossover function

# Realizes the crossover operation from parents to form new offspring. Splitting is realized at a random place 
crossover <- function(chromosome1, chromosome2){
  if (length(chromosome1)!=length(chromosome2)){
    print("Chromosomes are not of the same size")
  }
  else{
    #find the split spot  
    splitPositionMin = 1
    splitPositionMax = length(chromosome1)
    splitPosition = floor(runif(1, splitPositionMin, splitPositionMax))
    print(paste("Splitting occurs after position",splitPosition))  
    
    #split each chromosome
    preChromosome1 = chromosome1[1:splitPosition]
    postChromosome1 = chromosome1[(splitPosition+1):length(chromosome1)]
    preChromosome2 = chromosome2[1:splitPosition]
    postChromosome2 = chromosome2[(splitPosition+1):length(chromosome1)]
    
    #recombine chromosomes
    newChromosome1 = c(preChromosome1, postChromosome2)
    return(newChromosome1)
  } 
}

# Mutation function

mutation <- function(chromosome, mutationProbability = 0.01){
     n = mutationProbability     
     mutationFun <- function(chromosome, mutationProbability = n){          
          #generate associated uniform random variable for each locus
          mutationLocus = runif(length(chromosome),0,1)          
          #mutation occurs if r.v. < mutationProbability
          #find the location of mutation
          mutationOccur = mutationLocus < mutationProbability         
          #return the final result
          if (length(which(mutationOccur==T)) > 0) {
               print(paste("Mutation occurred at position", which(mutationOccur==T)))
          }        
          return((mutationOccur + chromosome) %% 2)
     }
     mutated <- as.data.frame(t(sapply(1:dim(chromosome)[1], function(x) mutationFun(chromosome[x,]))))
     return(mutated)
}


# Breeding function

matching <- function(parentPopulation, selectionProbability){
  #get the index of the winning first parent
  indexParent1 = sample(x=1:length(selectionProbability), size=1, replace=TRUE, prob=selectionProbability)
  print("First parent is")
  print(indexParent1)
  parent1 = parentPopulation[indexParent1, ]
  
  #get a random second parent
  indexParent2 = indexParent1
  while (indexParent2==indexParent1){
    indexParent2 = sample(1:length(selectionProbability), size=1, replace=TRUE)
  }
  print("Second parent is")
  print(indexParent2)
  parent2 = parentPopulation[indexParent2, ]
  
  parents <- rbind(parent1, parent2)
  return(parents)
}

breeding <- function(parentPopulation){
  parentPopulation = as.matrix(parentPopulation)
  finalColumn = dim(parentPopulation)[2]
  ranking = parentPopulation[, finalColumn]
  P = dim(parentPopulation)[2]
  selectionProbability = (2 * ranking)/(P*(P+1))
  
  index=1
  print("Breeding number")
  print(index)
  parents = matching(parentPopulation, selectionProbability)
  offspringPopulation = crossover(parents[1,], parents[2,])
  
#the breeding operation is repeated until the offspring population size equals the parent population size
  while (index<dim(parentPopulation)[1]){
    index = index +1
    print("Breeding number")
    print(index)
    parents = matching(parentPopulation, selectionProbability)
    parent1 = parents[1,]
    parent2 = parents[2,]
    #do the crossover with the two selected parents
    offsprings = crossover(parent1, parent2)
    offspringPopulation = rbind(offspringPopulation, offsprings)
  }
  
  print("The final offspring population is")
  return(offspringPopulation)
}


# Replacement function

replacement <- function(rankedParent, rankedOffspring, replacementPercentage=0.4) {
     
     #Extract the desired percentage of the fittest Parent and Offspring population
     #floor() and ceiling() make sure that we are always returned with matrix of the same dimension
     topParentPopulation <- head(rankedParent, floor((1-replacementPercentage)*nrow(rankedParent)))
     topOffspringPopulation <- head(rankedOffspring,ceiling(replacementPercentage*nrow(rankedOffspring)))
     newPopulation <- rbind(topParentPopulation, topOffspringPopulation)
     return(newPopulation[,1:dim(newPopulation)[2]])
}


# Testing 

# Function performs formal testing 

#-----------------
# Test 1
#-----------------
fileName <- "cars.dat" 
dependentVariable <- "MPG"
dataFile <- readDataFile(fileName, sep="", dependentVariable)
# Copy the original data file so updates can be made while preserving original.
dataOriginal <- as.data.frame(dataFile[1])
# Remove the first column, as it is a list of unique strings in this data set, and thus not relevant.
dataOriginal <- dataOriginal[-1]
dataParameter <- as.data.frame(dataFile[2])
# Remove column containing dependent variable
dataParameter[1] <- NULL

# Test a LM using AIC criteria
AICGAtest <- select(initialPopulationSize=25, numberGenerations=10, regressionType = "LM", replacementPercentage=.5, mutationProbability = .01)

# Test a GLM using BIC criteria
BICGAtest <- select(initialPopulationSize=25, numberGenerations=10, regressionType = "GLM", criterion = "BIC", replacementPercentage=.6)

testing <- function(criterion, data=dataOriginal, dependentVariable, regressionType) {
  formulatest <- as.formula(paste(dependentVariable, "~", paste(names(dataParameter), collapse="+")))
  if(regressionType!="LM" && regressionType!="GLM") {
    print("Please use only LM or GLM as regression type")
  } else {
    if(criterion!="AIC" && criterion!="BIC") {
      print("Please use only AIC or BIC as criterion type")
    } else if(criterion=="AIC") {
      if(regressionType=="LM") {
        lm2 <- lm(formulatest, data = dataOriginal)
        step <- stepAIC(lm2, direction="both")
        return(step)
      } else {
        glm2 <- glm(formulatest, data = data)
        step <- stepAIC(glm2, direction="both")
        return(step)
      } 
      # Note here that the AIC function is bering used to calculate BIC, with the value of k
      # manipulated appropriately.
    } else {
      if(regressionType=="LM") {
        lm2 <- lm(formulatest, data = data)
        step <- stepAIC(lm2, direction="both", k=log(nrow(dataOriginal)))
        return(step)
      } else {
        glm2 <- glm(formulatest, data = data)
        step <- stepAIC(glm2, direction="both", k=log(nrow(dataOriginal)))
        return(step)
      }
    }
  }
}

testing(criterion="AIC", regressionType="LM", dependentVariable="MPG", data=dataOriginal)
AICGAtest

# Test 2

testing(criterion="BIC", regressionType="GLM", dependentVariable="MPG", data=dataOriginal)
BICGAtest

# Plotting

# Create plots of the tests in order to see at what generation convergence occurred.

selectPlot <- function(initialPopulationSize, numberGenerations, regressionType, criterion="AIC", mutationProbability, replacementPercentage) { 
     parentPopulation <- initiate(dataParameter, initialPopulationSize)
     scoredParent <- runModel(parentPopulation, dependentVariable, regressionType, criterion)
     sortedParent <- mat.sort(scoredParent, ncol(scoredParent), decreasing=FALSE)
     if(length(which(sortedParent[ncol(sortedParent)] == 0)) > 0){
     parentPlot <- sortedParent[-which(sortedParent[ncol(sortedParent)] == 0),]}
     else{parentPlot <- sortedParent}
     parentPlotdat <- cbind(rep(0,nrow(parentPlot)),parentPlot[ncol(sortedParent)])
     plot(parentPlotdat, xlim = c(0,numberGenerations+1), xlab = "Generation", ylab = "Fitness Score", main = "Result of GA", pch=20)
     for(i in 1:numberGenerations) { 
          scoredParent <- runModel(parentPopulation, dependentVariable, regressionType, criterion)
          sortedParent <- mat.sort(scoredParent, ncol(scoredParent), decreasing=FALSE)
          offspringPopulation <- breeding(sortedParent)
          offspringPopulation <- offspringPopulation[,-ncol(offspringPopulation)] 
          mutatedOffspring<- mutation(offspringPopulation)
          scoredOffspring <- runModel(mutatedOffspring, dependentVariable, regressionType, criterion)
          sortedOffspring <- mat.sort(scoredOffspring, ncol(scoredOffspring), decreasing=FALSE)
          newGeneration <- replacement(sortedParent, sortedOffspring, replacementPercentage)
          parentPopulation <- newGeneration[-dim(newGeneration)]
          if(length(which(newGeneration[ncol(newGeneration)] == 0))> 0) {
          plotNewGeneration <- newGeneration[-which(newGeneration[ncol(newGeneration)] == 0),]}
          else{plotNewGeneration <- newGeneration}
          plotdat <- cbind(rep(i,nrow(plotNewGeneration)), plotNewGeneration[ncol(newGeneration)])
          points(plotdat,pch=20)
     }
}

AICGAtestplot <- selectPlot(initialPopulationSize=25, numberGenerations=10, regressionType = "LM", replacementPercentage=.5, mutationProbability = .01)
BICGAtestplot <- selectPlot(initialPopulationSize=25, numberGenerations=10, regressionType = "GLM", criterion = "BIC", replacementPercentage=.6)
