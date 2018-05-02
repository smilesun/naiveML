set.seed(1L)
library("R6")

# Queue
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

# Tree
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
      lp = length(l) / (length(l) + length(r))
      rp = length(r) / (length(l) + length(r))
      gini = rp * self$gini(r) + lp * self$gini(l)
      if (length(r) < self$minLeafsize || length(l) < self$minLeafsize) {
        return(NA)
      }
      list(l = l, r = r, gini = gini)
    })
    gini.list = lapply(res.list, function(x) {
      if (any(is.na(x))) return(NA)
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
      if (any(is.na(x))) return(NA)
      x$gini
    })
    id = which.min(b.list)
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

# Modified tree which allows for sampling features at each split
RMtree = R6Class("RMtree",
  inherit = Mtree,
  public = list(
  msplit = function(ind) {
    if (length(ind) < self$minInst) return(NULL)
    xns = names(self$X)
    xns = sample.int(n = length(xns), size = self$mtry)   # key difference here
    # at each split, only random sampled columns are used, this is different from a normal tree
    s.list = lapply(xns, self$splitOneColumn, ind = ind)
    b.list = lapply(s.list, function(x) {
      if (any(is.na(x))) return(NA)
      x$gini
    })
    id = which.min(b.list)
    return(list(l = s.list[[id]]$l, r = s.list[[id]]$r, fun = s.list[[id]]$fun, env = s.list[[id]]$env))
  }
    )
  )
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

    train = function() {
      self$model = lapply(1:self$ntree, function(i) {
        ind = sample.int(nrow(self$X), replace = TRUE)  # generate a bootstrap sampling of rows for eah tree
        X = self$X[ind, ]
        y = self$y[ind]
        mt = Mtree$new(X, y)
        mt$train()
        list(model = mt)
    })
    },

    pred.scalar = function(x) {
      res.list = lapply(self$model, function(m) {
        m$model$pred.scalar(x)
    })
      self$mode(unlist(res.list))   # FIXME: only works with factor now
    },

    pred = function(X) {
      apply(X, 1L, self$pred.scalar)
    }

    )
  )

## mlr data loading
library(mlr)
trainset = subsetTask(pid.task, subset = 1:760)
testset = subsetTask(pid.task, subset = 761:768)
train = getTaskData(trainset, target.extra = TRUE)
test = getTaskData(testset, target.extra = TRUE)
test = test$data


# Tree
mt = Mtree$new(train$data, train$target)
mt$train()
mt$pred(test)

## Random Forest
mrf = mRF$new(train$data, train$target, ntree = 10, mtry = 3)
mrf$train()
mrf$pred(test)

## mlr tree and Random Forest

lrn = makeLearner("classif.rpart")
mod = train(lrn, trainset)
pred = predict(mod, testset)
pred$data
lrn = makeLearner("classif.randomForest", ntree = 10, mtry = floor(sqrt(8)))
mod = train(lrn, trainset)
pred = predict(mod, testset)
pred$data
