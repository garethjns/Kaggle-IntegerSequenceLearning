# Kaggle Integer Sequence Learning
[Data and description](https://www.kaggle.com/c/integer-sequence-learning) | [Solutions post](http://blog.kaggle.com/2016/11/21/integer-sequence-learning-competition-solution-write-up-team-1-618-gareth-jones-laurent-borderie/)

- Final position: 17/269 

# Scripts

**classifyingSequences.ipynb** - Notebook proposing automatic classification of basic properties of sequences, which can be used to guide which approaches are more likely to work for a sequence. 

**BestRun.R** - Final submission before merging in alternative public kernel solutions. See below and [Kaggle forums](https://www.kaggle.com/c/integer-sequence-learning/forums/t/24971/solutions) for details.


# Aim
Predict the next value in integer sequences from the [Online Encyclopedia of Integer Sequences](http://oeis.org/).

# Methods
Sequentially test sequences with different solvers, ordered approximately by assumed reliability. 

## Solvers
Sequence solvers were applied in the order listed, based on assumed reliability. More detailed descriptions of each below. The brackets show the number of apparent solutions provided by each solver.

 - [Common differences](http://www.purplemath.com/modules/nextnumb.htm) (**2820**)
 - Common differences with variable step size (**2750** reliable and **9246** "dodgy")
 - Pattern search (**1785**)
 - Pattern search on common difference levels (**435**)
 - [Recurrence relation](https://en.wikipedia.org/wiki/Recurrence_relation) (**9849**)
 - Linear fitting using previous points (**37143**)
 - Non-linear fitting using previous points (**9000**)
 - Borrowed fallbacks (**<4000**)
 - [Mode-fallback](https://www.kaggle.com/wcukierski/integer-sequence-learning/mode-benchmark/run/255053/code) (**40196**)

### Common differences (diffTableSolve, diffTablePredict)
Take the difference between each adjacent term in the sequence, if the differences are all the same, the next term can be predicted. It's a special case of the recurrence relation (and presumably produces redundant predictions when both are used).

For example: 
```
Sequence = [2, 4, 6, 8, 10]
First differences = [4-2, 6-4, 8-6, 10-8] = [2, 2, 2, 2]

Next term = 10 + 2 = 12
```

If the differences aren't the same, continue down: 
``` 
Sequence = [2, 4, 7, 11, 16] 
1st diffs = [4-2, 7-4, 11-7, 16-11] = [2, 3, 4, 5] 
2nd diffs = [3-2, 4-3, 5-4] = [1, 1, 1]

Next term = 16 + 5 + 1 = 22
```

Each level of differences decreases in length by 1, and false positives are a risk when too few values are available at level to be sure the level is constant.

### Common differences with variable spacing (diffTableSolve2, diffTablePredict2)
An extension of the method of common differences is to take differences over a step size > 1, instead of from only adjacent terms.

For example, with a step of 2: 
```
Sequence = [2, 6, 8, 12] 
1st diffs = [8-2, 12-6] = [6, 6]

Next term  = ((n+1)-step) + diff = 8 + 6 = 14
```

In this case, the length of each level reduces by the step size, making false positives are a greater risk with this approach. "Dodgy" solutions were ones that were possible solutions, but were shorter than some threshold. If no solutions were proposed by pattern search or pattern search on difference levels, the dodgy solution was used.

### Pattern search (patternSearch)
- Starting with the second half of a sequence, pattern search compares it to the first half.
 - If the proportion of matched values is greater than some threshold, the match is used to extend the sequence.
 - If there's no match, one terms are dropped sequentially from the middle of the sequence (not the end) and compared to the (growing) start of the sequence.

### Pattern search on common difference levels (diffTablePattern)
This does a pattern search on each common difference level. In theory it might be able to find patterns in different levels, even if the difference levels never converge to a constant value... Maybe.

### Recurrence relation (RRSolve)
Based on [https://www.kaggle.com/ncchen/integer-sequence-learning/recurrence-relation/notebook](https://www.kaggle.com/ncchen/integer-sequence-learning/recurrence-relation/notebook)

### Linear fitting (fitModPP)
Attempts to fits polynomial to sequence using rolling window of n previous points. The last point of the sequence is held out and multiple functions fit. The held out point is used to assess the fits and determine the best, if accuracy is above a certain threshold, the winning function is refit on the entire sequence. The next (unknown term) is then predicted from this fit.

### Mode-fallback
If no solution was proposed by a solver, the mode of the sequence was used, as per the [benchmark](https://www.kaggle.com/wcukierski/integer-sequence-learning/mode-benchmark/run/255053).

# Conclusions
The application of specific solvers is effective to some extent, but still leaves a large proportion of the sequences unsolved (ie. mode fallback used) or incorrectly predicted.

## Limitations

### False positives
False positives were a significant danger. Fit an infinite number of non-linear functions and an infinite number will perfectly describe the sequence, but it doesn't mean it's correct one, and if it isn't, it's very unlikely to have any predictive value. Competition scoring was binary accuracy - the predicted term exactly correct or wrong.

### Generality

Looking at the number of mode-fallbacks (~40,000) and linear fits (~37,000) it's clear that the majority of the sequences are constructed by still-unknown functions (unknown in a world where the OEIS doesn't exist, that is).

For the linear fits, 10496 of the fits were scored perfectly, meaning ~27,000 only represent a polynomial estimations of another function. In addition, an unknown proportion of the 10496 "perfect" fits will be correct by chance, meaning only <1/3rd of the ~37,000 linear fits are likely to be true linear polynomial functions.

It might be possible to implement more solvers that would find more of these unknown sequences.

## Improvements

### Reliability
Solvers were applied on the basis of assumed reliability - this was measured (roughly) by how their implementation and priority affected overall score. A better approach would be to run each solver individually on a set of data and assess in isolation.

### Sequence classification
It's also likely that solver reliability varies as a function of sequence class. In fact, a sensible classification mechanism for a sequence is the best approach to predict the next term (although this is outside the scope of the Kaggle challenge). It might be sensible to determine certain basic properties of sequences first, then cluster the sequences in an unsupervised fashion. Solver reliability could then be assessed per-cluster and solvers applied in a guided fashion depending on the clusters properties. See this **classifyingSequences.ipynb** / [here](https://www.kaggle.com/garethjns/integer-sequence-learning/classifying-tagging-sequences) for further discussion and examples.

[Sequence classification](images/Figure2.png "Logo Title Text 1")

