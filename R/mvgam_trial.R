## adapting code from https://nicholasjclark.github.io/mvgam/articles/shared_states.html

library(mvgam)
library(mgcv)
library(data.table)


#use portal_data to start
fmod <- mvgam(
  captures ~ series - 1 + min_temp
  trend_formula = ~ s(time)
)






# flyway data?
# flyway <- fread(
#   "/run/media/obrien/OS/Users/darpa2/Analysis/OWF-flyways/data/joint_flyway_latitudes.csv"
# )
# # fill out missing obs
# flyway <- flyway[
#   data.table(expand.grid(
#     week = 1:53,
#     species = c("NARW", "striped bass", "northern gannet")
#   )),,
#   on = c('week', 'species')
# ]

# flyway[, series := as.factor(species)]
# flyway[, time := week]

# trend_map <- data.table(
#   series = unique(flyway$species),
#   trend = 1
# )

# mod <- mvgam(
#   mean_lat ~ series - 1,
#   trend_formula = ~ s(time, bs = 'cc'),

#   trend_model = GP(),
#   noncentred = T,

#   trend_map = trend_map,
#   family = gaussian(),
#   data = flyway,
#   run_model = F
# )
