Hw1
================
Brandon Turner
2/4/2021

### 1) What is the mean and variance for the loaded dice?

``` r
# mean
diceNumbers = seq(1,6)
diceProb = c(rep(0.1,5), 0.5)
diceMean = sum(diceNumbers * diceProb)

#variance
diceVariance = 0
for (i in 1:6){
  diceVariance = diceVariance + (diceProb[i] * ((diceNumbers[i] - diceMean)**2))
}

diceMean
```

    ## [1] 4.5

``` r
diceVariance
```

    ## [1] 3.25

### 2) Make a function in R that “rolls” this dice; return a vector containing the rolls. So if I call: myRolls &lt;- rollLoadedDie(10000) I would get a vector of size 10,000 that contains the rolls of my loaded die.

``` r
rollLoadedDie <- function(rolls){
  outputVector = c()
  for (i in 1:rolls){
    outputVector[i] = sample(1:6,1,prob = diceProb)
  }
  return(outputVector)
}
rollLoadedDie(10)
```

    ##  [1] 6 6 5 6 6 6 6 5 2 3

### 3) Make a histogram of some large number of rolls. Do the rolls of the loaded die approximate a uniform distribution?

``` r
par(mfrow=c(1,2))
hist(rollLoadedDie(10000))
hist(runif(10000))
```

![](HW1_files/figure-markdown_github/unnamed-chunk-3-1.png)

Clearly the histogram does not look uniformly distributed

### 4) Modify the code on Slide \#58 of lecture \#2 so that the means vs. trial size plots are from the loaded die. Generate these plots a few times. How many rolls appear to be necessary to get convergence on the expected values for the mean and variance?

``` r
trialSizes = c(5, 10, 15, 20, 25, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 100000)
means = c()
variances = c()

# trials
for (i in 1:length(trialSizes)){
  rolls = c()
  # make a vector of size trialSize[i]
  for (j in 1:trialSizes[i]){
    rolls[j] = sample(1:6, 1, prob = diceProb)
  }
  # mean/var of those trials
  means[i] = mean(rolls)
  variances[i] = var(rolls)
}

# plots 
par(mfrow=c(1,2))
plot(log10(trialSizes), means)
lines(log10(trialSizes), rep(diceMean, length(trialSizes)))
plot(log10(trialSizes), variances)
lines(log10(trialSizes), rep(diceVariance, length(trialSizes)))
```

![](HW1_files/figure-markdown_github/unnamed-chunk-4-1.png)

With the eye check it looks like it takes about 500-1000 rolls for the sample mean/var of the trials to converge on the expected mean and var
