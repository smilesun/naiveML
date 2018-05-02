set.seed(1L)
library("R6")
MQueue = R6Class("MQueue",
  public = list(
    queue = NULL,
    len = NULL,
    front = NULL,
    initialize = function() {
      self$queue = list()
      self$len = 0L
      self$front = 1L
    },

    get = function(cur) {
      return(self$queue[[cur]])
    },
    enqueue = function(ele) {
      self$len = length(self$queue)
      self$queue[[self$len + 1L]] = ele
      self$len = length(self$queue)
    },

    dequeue = function() {
      if (self$front <= self$len) {
        self$front = self$front + 1L
        return(self$queue[[self$front - 1L]])
      }
    },

    update = function(ele) {
      front = self$front
      if (front > 1L) {
        self$queue[[front - 1L]] = c(self$queue[[front - 1L]], ele)
        return(TRUE)
      }
      return(FALSE)
    },

    isEmpty = function() {
      return(self$front > self$len)
    }
    )
  )

Mtree = R6Class("Mtree",
  public = list(
  queue = NULL,
  minInst = NULL,
  minLeafsize = NULL,
  X = NULL,
  y = NULL,
  initialize = function(X, y) {
    self$queue = MQueue$new()
    self$minInst = 20L
    self$minLeafsize = 6L
    self$X = X
    self$y = y
  },

  gini = function(ind) {
    y = self$y[ind]
    mlevels = levels(y)
    n = length(y)
    p.list = lapply(mlevels, function(l) {
      (sum(y == l) / n) ^ 2
})
    1 - Reduce(sum, p.list)
  },

  splitOneColumn = function(colname, ind) {
    X = self$X
    x = X[ind, colname]
    all.cut = unique(x)
    res.list = lapply(all.cut, function(c) {
      r = ind[(X[ind, colname] > c)]
      l = ind[(X[ind, colname] <= c)]
      gini = self$gini(r) + self$gini(l)
      if (length(r) < self$minLeafsize || length(l) < self$minLeafsize) {
        return(NA)
      }
      list(l = l, r = r, gini = gini)
    })
    gini.list = lapply(res.list, function(x) {
      if (is.na(x)) return(NA)
      x$gini
    })
    if (all(is.na(gini.list))) return(NA)
    nd = which.min(gini.list)
    return(list(l = res.list[[nd]]$l, r = res.list[[nd]]$r, gini = gini.list[[nd]], fun = function(x, cut, colname) {
      x[colname] > cut}, env = list(cut = all.cut[nd], colname = colname)))
  },

  msplit = function(ind) {
    if (length(ind) < self$minInst) return(NULL)
    xns = names(self$X)
    s.list = lapply(xns, self$splitOneColumn, ind = ind)
    b.list = lapply(s.list, function(x) {
      if (is.na(x)) return(NA)
      x$gini
    })
    id = which.max(b.list)
    return(list(l = s.list[[id]]$l, r = s.list[[id]]$r, fun = s.list[[id]]$fun, env = s.list[[id]]$env))
  },

  train = function() {
    flag = TRUE
    X = self$X
    y = self$y
    ind = 1:nrow(X)
    self$queue$enqueue(list(ind = ind, target = y[ind]))
    while (flag) {
      if (self$queue$isEmpty()) return(flag)
      ele = self$queue$dequeue()
      ind = ele$ind
      res = self$msplit(ind)
      if (is.null(res)) {
        next
      }
      ele = list(lind = res$l, rind = res$r, fun = res$fun, lpos = self$queue$front, rpos = self$queue$front + 1L, env = res$env)
      self$queue$update(ele)
      self$queue$enqueue(list(ind = res$l, target = y[res$l]))
      self$queue$enqueue(list(ind = res$r, target = y[res$r]))
    }
  },

  mode = function(x) {
      ta = table(x)
      tam = max(ta)
      if (all(ta == tam))
          mod = NA
      else
          if (is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
      else
          mod = names(ta)[ta == tam]
      return(mod)
  },

  majorvote = function(x) {
    self$mode(x)
  },

  pred.scalar = function(x) {
    cur = 1L
    while (TRUE) {
      ele = self$queue$get(cur)
      fun = ele$fun
      env = ele$env
      if (is.null(fun)) {
        return(self$majorvote(ele$target))
      }
      cur = ifelse(eval(fun(x, env$cut, env$colname)), ele$rpos, ele$lpos)
    }
  },

  pred = function(X) {
    apply(X, 1, self$pred.scalar)
  }
  )
)

library(mlr)
trainset = subsetTask(pid.task, subset = 1:760)
testset = subsetTask(pid.task, subset = 761:768)
res = getTaskData(trainset, target.extra = TRUE)
mt = Mtree$new(res$data, res$target)
mt$train()
test = getTaskData(testset, target.extra = TRUE)
test = test$data
mt$pred.scalar(test[1, ])
mpred = mt$pred(test)
lrn = makeLearner("classif.rpart")
mod = train(lrn, trainset)
hpred = predict(mod, testset)
hpred$data

mRF = R6Class("mRF",
  public = list(
    ntree = NULL,
    mtry = NULL,
    X = NULL,
    y = NULL,
    ncols = NULL,
    model = NULL,
    initialize = function(X, y, ntree, mtry) {
      self$ntree = ntree
      self$mtry = mtry
      self$X = X
      self$y = y
      self$ncols = ncol(self$X)
    },

    mode = function(x) {
      ta = table(x)
      tam = max(ta)
      if (all(ta == tam))
          mod = NA
      else
          if (is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
      else
          mod = names(ta)[ta == tam]
      return(mod)
  },



    sample = function() {
      sample.int(n = self$ncols, size = self$mtry)
    },

    train = function() {
      self$model = lapply(1:self$ntree, function(i) {
        cols = self$sample()
        mt = Mtree$new(self$X[, cols], self$y)
        mt$train()
        list(model = mt, cols = cols)
    })
    },

    pred.scalar = function(x) {
      res.list = lapply(self$model, function(m) {
        m$model$pred.scalar(x[m$cols])
    })
      self$mode(unlist(res.list))   # FIXME: factor
    },

    pred = function(X) {
      apply(X, 1L, self$pred.scalar) 
    }

    )
  )

# mtry is the the number of features subsampled
lrn = makeLearner("classif.randomForest", ntree = 10, mtry = floor(sqrt(8)))
trainset = subsetTask(pid.task, subset = 1:760)
# I made the first 760 instances to be training and the rest to be test
testset = subsetTask(pid.task, subset = 761:768)
mod = train(lrn, trainset)
pred = predict(mod, testset)
pred$data
data = getTaskData(pid.task)  # you can use this dataset for proofing your code



library(mlr)
trainset = subsetTask(pid.task, subset = 1:760)
testset = subsetTask(pid.task, subset = 761:768)
res = getTaskData(trainset, target.extra = TRUE)
mt = mRF$new(res$data, res$target, ntree = 10, mtry = 3)
test = getTaskData(testset, target.extra = TRUE)
test = test$data
mt$train()
mt$pred.scalar(test[1, ])
mpred = mt$pred(test)
mpred
