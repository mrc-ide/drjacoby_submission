// saved as 8schools.stan
data {
  real gamma;
}
parameters {
  real mu;
}
model {
  target += -gamma*(mu^2 - 1)^2;
}
