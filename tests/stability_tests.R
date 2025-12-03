set.seed(1)
n <- 120; p <- 8
X <- matrix(rnorm(n*p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p-5))
y <- X %*% beta + rnorm(n, sd = 1)

out <- stability(X = X, y = as.numeric(y), 
                  B = 50, 
                  resample_type = "bootstrap",
                  family = "gaussian", 
                  K = 10, 
                  eps = 1e-6, 
                  delta = 1, 
                  L = 50)

is.numeric(out$pi) # check if numeric
all(out$pi >= 0 & out$pi <= 1) # check if all values between [0, 1]
length(out$pi) == ncol(X) # check if length = ncol(x)
nrow(out$z) == out$meta$B # check if z matches dimension B x p
ncol(out$z) == ncol(X)

# ------------------

set.seed(123)
n <- 300
p <- 10

X <- as.data.frame(matrix(rnorm(n*p), n, p))
colnames(X) <- paste0("x", 1:p)

# Strong signal from x1 and x2; rest are pure noise
y <- 3*X$x1 - 2*X$x2 + rnorm(n, sd = 1)

stab <- stability(
  X, y,
  family = "gaussian",
  B = 40,
  resample_type = "bootstrap",
  K = 5,
  seed = 1
)

stab$pi # check that stab$pi[1] and stab$pi[2] are the largest values in the set, 
        # and that all other values are low due to no true predictive signal

order(stab$pi, decreasing = TRUE)

# Signal features should have clearly higher stability
stab$pi[1] > 0.3
stab$pi[2] > 0.3

# Noise features should be small
mean(stab$pi[3:10]) < 0.1
