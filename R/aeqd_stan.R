library(cmdstanr)
library(data.table)

# Import striped bass data
sb <- fread(
  "/run/media/obrien/OS/Users/darpa2/Analysis/OWF-flyways/data/joint_flyway_latitudes.csv"
) |>
  _[species == 'striped bass']

# Set up data for Stan program
stan_data <- list(
  N = nrow(sb),
  lat_rad = sb$mean_lat * pi / 180,
  lon_rad = sb$mean_lon * pi / 180
)

# Compile
mod <- cmdstan_model('Stan/aeqd_pca.stan')

fit <- mod$sample(
  data = stan_data,
  seed = 804,
  chains = 2,
  parallel_chains = 4,
  refresh = 500
)


mcmc_hist(
  fit$draws(c("lat0", "lon0")),
  transformations = function(x) {
    x * 180 / pi
  },
  facet_args = list(ncol = 1)
)


proj_center <- fit$draws(c("lon0", "lat0"), format = "df")[, 1:2] * 180 / pi

library(sf)
proj_center <- st_as_sf(k, coords = c("lon0", "lat0"), crs = 4326)
mapview::mapview(proj_center, alpha = 0.1, color = NA) +
  mapview::mapview(
    st_sf(st_sfc(
      st_point(fit$summary(c("lon0", "lat0"))$mean * 180 / pi),
      crs = 4326
    )),
    col.regions = 'red',
    alpha = 1
  )

## Fitting k
mod_k <- cmdstan_model('Stan/aeqd_k_pca.stan')

fit_k <- mod_k$sample(
  data = stan_data,
  seed = 804,
  chains = 2,
  parallel_chains = 4,
  refresh = 500
)
bayesplot::mcmc_hist(
  fit_k$draws("k_prime")
)
