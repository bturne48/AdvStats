HW5
================
Brandon Turner
4/12/2021

## 1) We will utilize our RNA seq dataset of E. Coli genes from mice. Read and normalize the counts table ( “nc101\_scaff\_dataCounts.txt “ into R).

``` r
# read in data
myData = read.table('/Users/bturne48/Documents/GitHub/AdvStats/HW5/longitdunalRNASeqData/nc101_scaff_dataCounts.txt', header = T, row.names = 1, sep ='\t', colClasses = c('character', rep('numeric', 11)))

# remove rare genes
myData = myData[apply(myData,1,median) > 5, ]

# normalize
myData_norm = myData
for (i in 1:ncol(myData)){
  colSum = sum(myData[,1])
  myData_norm[,i] = myData_norm[,i]/colSum
}
head(myData_norm)
```

    ##                    D2_01        D2_02        D2_03       W12_01       W12_02
    ## NC101_00003 5.680896e-04 5.007477e-05 6.906864e-05 1.674915e-04 0.0007304009
    ## NC101_00004 5.180148e-05 3.453432e-06 1.036030e-05 1.208701e-05 0.0001433174
    ## NC101_00005 2.503738e-04 9.151595e-05 3.971447e-05 2.106594e-04 0.0004679401
    ## NC101_00006 8.806252e-04 1.381373e-04 9.324267e-05 2.693677e-04 0.0012380554
    ## NC101_00007 8.633580e-06 5.180148e-06 0.000000e+00 2.244731e-05 0.0000431679
    ## NC101_00008 2.417403e-05 6.906864e-06 3.453432e-06 1.381373e-05 0.0001122365
    ##                   W12_03       w20_01       w20_02       w20_03       w20_04
    ## NC101_00003 2.434670e-04 8.098298e-04 5.905369e-04 1.526417e-03 1.231149e-03
    ## NC101_00004 5.870835e-05 1.312304e-04 7.770222e-05 2.261998e-04 2.192929e-04
    ## NC101_00005 2.210197e-04 1.761250e-04 9.842282e-05 2.313800e-04 1.933922e-04
    ## NC101_00006 3.332562e-04 8.443642e-04 3.936913e-04 9.255198e-04 8.409107e-04
    ## NC101_00007 1.726716e-06 1.554044e-05 3.453432e-06 5.180148e-05 1.554044e-05
    ## NC101_00008 2.417403e-05 3.591569e-04 1.415907e-04 5.214683e-04 4.161386e-04
    ##                   w20_05
    ## NC101_00003 8.978924e-05
    ## NC101_00004 1.122365e-04
    ## NC101_00005 2.244731e-05
    ## NC101_00006 9.842282e-05
    ## NC101_00007 1.899388e-05
    ## NC101_00008 2.417403e-05

## 2) For every row in the normalized spreadsheet, run three t-tests ( “day 2” vs. “week 12”, “day 2” vs. “week 18” and “week 12” vs. “week 18”. At a p &lt; .05 threshold fill in the following table:

``` r
# uncorrected ounts
y1_uncor = 0
y1_pval = c()
y2_uncor = 0
y2_pval = c()
y3_uncor = 0
y3_pval = c()
for (i in 1:nrow(myData_norm)){
  # day 2 vs week 12
  x1 = t.test(myData[i,1:3], myData[i,4:6])
  y1_pval[i] = x1$p.value
  if (x1$p.value < 0.05){
    y1_uncor = y1_uncor + 1
  }
  
  # day 2 vs week 18
  x2 = t.test(myData[i,1:3], myData[i,7:11])
  y2_pval[i] = x2$p.value
  if (x2$p.value < 0.05){
    y2_uncor = y2_uncor + 1
  }
  
  # week 12 vs week 18
  x3 = t.test(myData[i,4:6], myData[i,7:11])
  y3_pval[i] = x3$p.value
  if (x3$p.value < 0.05){
    y3_uncor = y3_uncor + 1
  }
  
}
uncorrected_res = c(y1_uncor, y2_uncor, y3_uncor)


# FDR corrected
fdr_corrected_res = c(sum(p.adjust(y1_pval, 'fdr') < 0.1), sum(p.adjust(y2_pval, 'fdr') < 0.1), sum(p.adjust(y3_pval, 'fdr') < 0.1))


# Bon corrected
threshold = 0.05/length(myData_norm)
bonf_corrected_res = c(sum(p.adjust(y1_pval, "bonferroni") <= threshold), sum(p.adjust(y2_pval, "bonferroni") <= threshold), sum(p.adjust(y3_pval, "bonferroni") <= threshold))


# results df 
endTable = data.frame(uncorrected_res, fdr_corrected_res, bonf_corrected_res)
endTable
```

    ##   uncorrected_res fdr_corrected_res bonf_corrected_res
    ## 1              70                 0                  0
    ## 2            2317              2540                  0
    ## 3             671                 1                  0

## 3) Make histograms of all the uncorrected p-values for each of the three companions. Are any of the distributions uniform?

``` r
hist(y1_pval, main ='Day 2 vs Week 12')
```

![](HW5_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
hist(y2_pval, main ='Day 2 vs Week 18')
```

![](HW5_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
hist(y3_pval, main ='Week 12 vs Week 18')
```

![](HW5_files/figure-markdown_github/unnamed-chunk-3-3.png)

None of the histrograms appear to show a unif dist

## 4) Based on these data, when is the biggest shift in the transcriptome? Which samples are most different from one another?

Day 2 vs Week 18 shows the most values in column 1 of the histogram having the highest count.
