# Internal helper: locate a file shipped in inst/extdata
.get_extdata_path <- function(fname, subdir = NULL) {
  if (is.null(subdir)) {
    f <- system.file("extdata", fname, package = "FastGibbsSamplers")
  } else {
    f <- system.file("extdata", subdir, fname, package = "FastGibbsSamplers")
  }
  if (identical(f, "")) {
    stop(
      "Cannot find file in inst/extdata: ", fname, "\n",
      "Place it under inst/extdata/", if (!is.null(subdir)) paste0(subdir, "/"), ".",
      call. = FALSE
    )
  }
  f
}

# Internal helper: safe model.matrix power expansion
.mm_power <- function(df, power = 2) {
  df <- as.data.frame(df)
  if (power == 2) {
    stats::model.matrix(~ .^2, data = df)[, -1, drop = FALSE]
  } else if (power == 3) {
    stats::model.matrix(~ .^3, data = df)[, -1, drop = FALSE]
  } else {
    stop("power must be 2 or 3.", call. = FALSE)
  }
}

generate_data <- function(dataset_name) {

  # ---- simulated datasets (unchanged) ----
  if (dataset_name == "sim1") {
    n <- 30000; p <- 30
    x <- matrix(stats::rnorm(n * p), n, p)
    y <- stats::rnorm(n)
    vbeta0 <- rep(c(-3, 3), 10)
    inds <- sort(sample.int(p, length(vbeta0)))
    y <- drop(y + x[, inds] %*% vbeta0)
  }

  if (dataset_name == "sim2") {
    n <- 500; p <- 400
    x <- matrix(stats::rnorm(n * p), n, p)
    y <- stats::rnorm(n)
    vbeta0 <- rep(c(-3, 3), 10)
    inds <- sort(sample.int(p, length(vbeta0)))
    y <- drop(y + x[, inds] %*% vbeta0)
  }

  if (dataset_name == "sim3") {
    n <- 100; p <- 500
    x <- matrix(stats::rnorm(n * p), n, p)
    y <- stats::rnorm(n)
    vbeta0 <- rep(c(-3, 3), 10)
    inds <- sort(sample.int(p, length(vbeta0)))
    y <- drop(y + x[, inds] %*% vbeta0)
  }

  if (dataset_name == "sim4") {
    n <- 50; p <- 5000
    x <- matrix(stats::rnorm(n * p), n, p)
    y <- stats::rnorm(n)
    vbeta0 <- rep(c(-3, 3), 10)
    inds <- sort(sample.int(p, length(vbeta0)))
    y <- drop(y + x[, inds] %*% vbeta0)
  }

  # ---- external-data generators (use inst/extdata instead of here("data",...)) ----
  if (dataset_name == "qtl") {
    res <- generate_qtl_data()
    y <- res$y
    x <- res$x
  }

  if (dataset_name == "covid") {
    res <- generate_covid_data()
    y <- res$y
    x <- res$x
    s <- apply(x, 2, sum)
    inds <- which(s == 0)
    if (length(inds) > 0) x <- x[, -inds, drop = FALSE]
  }

  if (dataset_name == "riboflavin") {
    res <- generate_riboflavin_data()
    y <- res$y
    x <- res$x
  }

  if (dataset_name == "cookie") {
    res <- generate_cookie_data()
    y <- res$y
    x <- res$x
  }

  # ---- Kakadu (from your shipped CSV) ----
  if (dataset_name == "Kakadu2") {
    f <- .get_extdata_path("Kakadu.csv")  # put this under inst/extdata/
    dat <- utils::read.csv(f, stringsAsFactors = TRUE)

    y <- as.vector(dat$income)
    x <- dat[, c(2:21, 23)]
    x <- .mm_power(x, power = 2)
  }

  if (dataset_name == "Kakadu3") {
    f <- .get_extdata_path("Kakadu.csv")
    dat <- utils::read.csv(f, stringsAsFactors = TRUE)

    y <- as.vector(dat$income)
    x <- dat[, c(2:21, 23)]
    x <- .mm_power(x, power = 3)
  }

  # ---- diabetes (CRAN dataset from lars) ----
  if (dataset_name == "diabetes") {
    if (!requireNamespace("lars", quietly = TRUE)) {
      stop("Dataset 'diabetes' requires the 'lars' package.", call. = FALSE)
    }
    data("diabetes", package = "lars")
    y <- lars::diabetes$y
    x <- lars::diabetes$x
  }

  if (dataset_name == "diabetes2") {
    if (!requireNamespace("lars", quietly = TRUE)) {
      stop("Dataset 'diabetes2' requires the 'lars' package.", call. = FALSE)
    }
    data("diabetes", package = "lars")
    y0 <- lars::diabetes$y
    x0 <- lars::diabetes$x

    norm <- normalize(y0, x0, scale = TRUE)
    x <- .mm_power(data.frame(x = norm$mX), power = 2)
    y <- norm$vy
  }

  # ---- energy (your shipped CSV) ----
  if (dataset_name == "energy2") {
    f <- .get_extdata_path("energydata_complete.csv") # move to inst/extdata/
    dat <- utils::read.csv(f)

    y0 <- dat$Appliances
    x0 <- dat[, -c(1, 2), drop = FALSE]

    norm <- normalize(y0, x0, scale = TRUE)
    x <- .mm_power(data.frame(x = norm$mX), power = 2)
    y <- norm$vy
  }

  if (dataset_name == "energy3") {
    f <- .get_extdata_path("energydata_complete.csv")
    dat <- utils::read.csv(f)

    y0 <- dat$Appliances
    x0 <- dat[, -c(1, 2), drop = FALSE]

    norm <- normalize(y0, x0, scale = TRUE)
    x <- .mm_power(data.frame(x = norm$mX), power = 3)
    y <- norm$vy
  }

  # ---- Crime (your shipped RData) ----
  if (dataset_name == "Crime") {
    f <- .get_extdata_path("comData.Rdata")  # move to inst/extdata/
    e <- new.env(parent = emptyenv())
    base::load(f, envir = e)

    # expects X and Y objects in that RData (as in your current code) :contentReference[oaicite:2]{index=2}
    X <- e$X
    Y <- e$Y

    mX <- stats::na.omit(X)
    mX <- mX[, -which(colnames(mX) %in% c("ownHousQrange", "rentUpperQ")), drop = FALSE]
    y <- Y[, "murders"]
    x <- mX
  }

  if (dataset_name == "Crime2") {
    f <- .get_extdata_path("comData.Rdata")
    e <- new.env(parent = emptyenv())
    base::load(f, envir = e)
    X <- e$X
    Y <- e$Y

    mX <- stats::na.omit(X)
    mX <- mX[, -which(colnames(mX) %in% c("ownHousQrange", "rentUpperQ")), drop = FALSE]
    y <- Y[, "murders"]
    x <- .mm_power(data.frame(x = mX), power = 2)
  }

  # ---- Hitters (CRAN dataset from ISLR) ----
  if (dataset_name == "Hitters") {
    if (!requireNamespace("ISLR", quietly = TRUE)) {
      stop("Dataset 'Hitters' requires the 'ISLR' package.", call. = FALSE)
    }
    data("Hitters", package = "ISLR")
    Hitters <- stats::na.omit(ISLR::Hitters)
    x <- stats::model.matrix(Salary ~ ., data = Hitters)[, -1, drop = FALSE]
    y <- as.numeric(Hitters[, "Salary"])
  }

  if (dataset_name == "Hitters2") {
    if (!requireNamespace("ISLR", quietly = TRUE)) {
      stop("Dataset 'Hitters2' requires the 'ISLR' package.", call. = FALSE)
    }
    data("Hitters", package = "ISLR")
    Hitters <- stats::na.omit(ISLR::Hitters)
    x <- stats::model.matrix(Salary ~ .^2, data = Hitters)[, -1, drop = FALSE]
    y <- as.numeric(Hitters[, "Salary"])
  }

  # ---- BostonHousing2 (CRAN dataset from mlbench) ----
  if (dataset_name == "BostonHousing2") {
    if (!requireNamespace("mlbench", quietly = TRUE)) {
      stop("Dataset 'BostonHousing2' requires the 'mlbench' package.", call. = FALSE)
    }
    data("BostonHousing", package = "mlbench")
    BostonHousing <- stats::na.omit(mlbench::BostonHousing)
    x <- stats::model.matrix(medv ~ .^2, data = BostonHousing)[, -1, drop = FALSE]
    y <- as.numeric(BostonHousing[, "medv"])
  }

  # ---- flare eyedata (CRAN dataset from flare) ----
  if (dataset_name == "eyedata") {
    if (!requireNamespace("flare", quietly = TRUE)) {
      stop("Dataset 'eyedata' requires the 'flare' package.", call. = FALSE)
    }
    data("eyedata", package = "flare")
    x <- flare::eyedata$x
    y <- flare::eyedata$y
  }

  # ---- final normalize (same as your current code) ----
  norm <- normalize(y, x, scale = TRUE)
  vy <- norm$vy
  mX <- norm$mX

  n <- length(vy)
  p <- ncol(mX)

  list(mX = mX, vy = vy, n = n, p = p)
}

# ---- external dataset helpers rewritten to use inst/extdata ----

generate_covid_data <- function() {
  meta_path <- .get_extdata_path("all_metadata.rds")
  dat_path  <- .get_extdata_path("all_bulk.rds")

  meta <- readRDS(meta_path)
  dat  <- readRDS(dat_path)

  inds_na <- which(is.na(meta$meta_WHO_scores))
  x <- t(dat[, -inds_na, drop = FALSE])
  y <- meta$meta_WHO_scores[-inds_na]

  y[y == "1 or 2"] <- 1.5
  y <- as.numeric(y)

  list(x = x, y = y, dataset_name = "covid")
}

generate_cookie_data <- function() {
  cookie_path <- .get_extdata_path("cookie_data.csv")
  cookie <- utils::read.csv(cookie_path)

  x <- as.matrix(cookie[, 1:700])
  y <- as.matrix(cookie[, 701:704])
  y <- y[, 2]

  list(x = x, y = y, dataset_name = "cookie")
}

generate_riboflavin_data <- function() {
  if (!requireNamespace("hdi", quietly = TRUE)) {
    stop("Dataset 'riboflavin' requires the 'hdi' package.", call. = FALSE)
  }
  data("riboflavin", package = "hdi")
  list(x = hdi::riboflavin$x, y = hdi::riboflavin$y, dataset_name = "riboflavin")
}

generate_qtl_data <- function() {
  resp_path <- .get_extdata_path("phe_simulat.csv")
  cov_path  <- .get_extdata_path("gen_simulat.csv")

  response <- utils::read.table(resp_path, header = FALSE, sep = ",")
  covariates <- utils::read.table(cov_path, header = FALSE, sep = ",")

  n <- nrow(response)
  p <- ncol(covariates)

  X <- matrix(0, n, 7381)
  vbeta <- numeric(7381)
  count <- 1

  for (i in 1:p) {
    for (j in i:p) {
      if (i == j) {
        X[, count] <- covariates[, i]
      } else {
        X[, count] <- covariates[, i] * covariates[, j]
      }

      # your original signal pattern (unchanged) :contentReference[oaicite:3]{index=3}
      vbeta[count] <- 0
      if ((i == 1)   & (j == 1))     vbeta[count] <- 4.47
      if ((i == 21)  & (j == 21))    vbeta[count] <- 3.16
      if ((i == 31)  & (j == 31))    vbeta[count] <- 2.24
      if ((i == 51)  & (j == 51))    vbeta[count] <- 1.58
      if ((i == 71)  & (j == 71))    vbeta[count] <- 1.58
      if ((i == 91)  & (j == 91))    vbeta[count] <- 1.10
      if ((i == 101) & (j == 101))   vbeta[count] <- 1.10
      if ((i == 111) & (j == 111))   vbeta[count] <- 0.77
      if ((i == 121) & (j == 121))   vbeta[count] <- 0.77
      if ((i == 1)   & (j == 11))    vbeta[count] <- 1.00
      if ((i == 2)   & (j == 119))   vbeta[count] <- 3.87
      if ((i == 10)  & (j == 91))    vbeta[count] <- 1.30
      if ((i == 15)  & (j == 75))    vbeta[count] <- 1.73
      if ((i == 20)  & (j == 46))    vbeta[count] <- 1.00
      if ((i == 21)  & (j == 22))    vbeta[count] <- 1.00
      if ((i == 26)  & (j == 91))    vbeta[count] <- 1.00
      if ((i == 41)  & (j == 61))    vbeta[count] <- 0.71
      if ((i == 56)  & (j == 91))    vbeta[count] <- 3.16
      if ((i == 65)  & (j == 85))    vbeta[count] <- 2.24
      if ((i == 86)  & (j == 96))    vbeta[count] <- 0.89
      if ((i == 101) & (j == 105))   vbeta[count] <- 1.00
      if ((i == 111) & (j == 121))   vbeta[count] <- 2.24

      count <- count + 1
    }
  }

  sigma2.true <- 20
  y_mu <- X %*% vbeta
  y <- as.vector(y_mu + stats::rnorm(n, 0, sqrt(sigma2.true)))

  X_std <- scale(X, center = TRUE, scale = TRUE)

  list(y = y, x = X_std, beta = c(0, vbeta))
}
