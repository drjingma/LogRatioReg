# testing sequential correlation tests
rm(list=ls())
alpha = 0.05

which.cors = 2 # 1-7
which.cards = 2 # 1-5
# when 2, 2, don't match up

cardinality.ests_options = matrix(
  c(
    c(5, 27, 28),
    c(27, 5, 28),
    c(28, 27, 5)
  ), byrow = TRUE, ncol = 3
)
cardinality.ests = cardinality.ests_options[which.cards, ]
cors_options = matrix( # let r1 = rmax = 0.5, but don't rank the others
  #   let 0.1 away be similar, 0.2 away be different
  c(
    c(0.5, 0.5, 0.5), # (0, 0, 0) -> (r1, r2, r3)
    # c(), # (0, 0, 1) -> (r1, r2), (r1, r3) -- not possible!
    c(0.4, 0.3, 0.5), # (0, 1, 0) -> (r1, r2), (r2, r3)
    c(0.3, 0.4, 0.5), # (1, 0, 0) -> (r1, r3), (r2, r3)
    c(0.5, 0.5, 0.3), # (0, 1, 1) -> (r1, r2), (r3)
    c(0.3, 0.5, 0.5), # (1, 0, 1) -> (r1, r3), (r2)
    c(0.3, 0.3, 0.5), # (1, 1, 0) -> (r1), (r2, r3)
    c(0.1, 0.5, 0.3) # (1, 1, 1) -> (r1), (r2), (r3)
  ), byrow = TRUE, ncol = 3
)
cors = cors_options[which.cors, ]

# method 1: truth table method #################################################
max.pair = which.max(abs(cors))
selected.pair = max.pair
z.scores = 0.5 * log((1 + abs(cors)) / (1 - abs(cors)))
sigma.z1.minus.z2 = 0.1 #sqrt(2 / (n - 3))

other.pairs = (1:3)[-max.pair]
# tests (I), (II), (III)
test.stats = c(
  z.scores[other.pairs] - z.scores[max.pair], 
  diff(z.scores[other.pairs])
) / sigma.z1.minus.z2
p.values = rep(NA, length(test.stats))
for(i in 1:length(test.stats)){ # since lower.tail doesn't take a vector
  p.values[i] = pnorm(q = test.stats[i], lower.tail = test.stats[i] < 0)
}
signif.diffs = p.values <= alpha
if(all(signif.diffs[c(1, 2)])){ # if (I) & (II) == 1, 
  # choose Bmax
} else if(all(!(signif.diffs[c(1, 2)]))){ # if(I) & (II) == 0, 
  # choose sparsest of the 3 balances
  selected.pair = which.min(cardinality.ests)
} else if(sum(signif.diffs[c(1, 2)]) == 1){ # if (I) xor (II) == 1 (other 0),
  # check which is bigger, Bmax or the other balance
  # (either B1 or B2)
  sim.other.pair = other.pairs[!signif.diffs[c(1, 2)]]
  if(cardinality.ests[sim.other.pair] < cardinality.ests[max.pair]){
    selected.pair = sim.other.pair
    if(!signif.diffs[3]){ # if (III) == 0, 
      notsim.other.pair = other.pairs[signif.diffs[c(1, 2)]]
      # check which is bigger between B1 and B2
      if(cardinality.ests[notsim.other.pair] < 
         cardinality.ests[sim.other.pair]){
        selected.pair = notsim.other.pair
      }
    }
  }
}
method1.selected.pair = selected.pair

# method 2: ranked sequential tests ############################################
max.pair = which.max(abs(cors))
selected.pair = max.pair
z.scores = 0.5 * log((1 + abs(cors)) / (1 - abs(cors)))
sigma.z1.minus.z2 = 0.1 #sqrt(2 / (n - 3))

# rank the absolute correlations and comptue their z.scores
ranks = rank(-abs(cors), ties.method = "first")
# desc.order = order(abs(cors), decreasing = TRUE)
# ordered tests, where |r1| >= |r2| >= |r3|
test.stats = c(
  z.scores[ranks == 1] - z.scores[ranks == 2], # test(|r1|, |r2|)
  z.scores[ranks == 2] - z.scores[ranks == 3], # test(|r2|, |r3|)
  z.scores[ranks == 1] - z.scores[ranks == 3] # test(|r1|, |r3|)
) / sigma.z1.minus.z2
p.values = rep(NA, length(test.stats))
for(i in 1:length(test.stats)){ # since lower.tail doesn't take a vector
  p.values[i] = pnorm(q = test.stats[i], lower.tail = test.stats[i] < 0)
}
signif.diffs = p.values <= alpha
# if |r1| & |r2| are different, choose B1 (with max |r|, i.e. |r1|).
if(!signif.diffs[1]){ 
  # if |r1| & |r2| are similar, check if |r2| & |r3| are similar
  if(signif.diffs[2]){
    # if |r1| & |r2| are similar, but |r2| & |r3| are different,
    #   choose the sparsest among B1 & B2
    if(cardinality.ests[ranks == 2] < cardinality.ests[ranks == 1]){
      selected.pair = which(ranks == 2)
    }
  } else {
    # if |r1| & |r2| are similar, and |r2| & |r3| are similar, 
    #   check if |r1| & |r3| are similar --
    #   -- if |r1| & |r3| are different, choose the sparsest among B1 & B2
    if(cardinality.ests[ranks == 2] < cardinality.ests[ranks == 1]){
      selected.pair = which(ranks == 2)
    }
    if(!signif.diffs[3]){
      # if |r1| & |r2| are similar, and |r2| & |r3| are similar, 
      #   and also |r1| & |r3| are different, 
      #   choose the sparsest among B1, B2, & B3
      if(cardinality.ests[ranks == 3] < cardinality.ests[selected.pair]){
        selected.pair = which(ranks == 3)
      }
    }
  }
}
method2.selected.pair = selected.pair

t(data.frame(
  Method =c(1, 2),
  Sel.Pair = c(method1.selected.pair, method2.selected.pair), 
  Sel.CardB = c(cardinality.ests[method1.selected.pair], cardinality.ests[method2.selected.pair])
))


