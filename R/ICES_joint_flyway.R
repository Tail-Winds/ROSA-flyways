library(data.table)
library(ggplot2)
library(qs)

# Striped bass
hud_sb <- qread('data/embargo/hud_detects.qs')

setDT(hud_sb)

hud_sb[, year := year(date.local)]
hud_sb[, week := week(date.local)]

sb_agg <- hud_sb[, .(
  mean_lat = mean(lat, na.rm = T),
  mean_lon = mean(long, na.rm = T),
  species = 'striped bass'
),
by = c('week', 'year', 'transmitter')]

ggplot(data = sb_agg) +
  geom_point(aes(x = week, y = mean_lat)) +
  facet_wrap(~year)

ggplot(data = sb_agg) +
  geom_point(aes(x = week, y = mean_lon)) +
  facet_wrap(~year) +
  coord_flip()



# sb_agg <- hud_sb[, .(mean_lat = mean(lat, na.rm = T),
#                      lsd_lat = mean(lat, na.rm = T) - sd(lat, na.rm = T),
#                      usd_lat = mean(lat, na.rm = T) + sd(lat, na.rm = T),
#                      species = 'striped bass'),
#                  by = c('week', 'yr')] 

# pick 2017
# ggplot(data= sb_agg) +
#   geom_errorbar(aes(ymin = lsd, ymax = usd, x =week)) +
#   facet_wrap(~yr) 

sb_agg <- sb_agg[year == 2017, 
                 .(species, week, mean_lat, mean_lon)]


# Gannet
# gannet <- fread('extracted_figures/gannet_extracted.csv',
#                 col.names = c('date', 'latitude'))
# ## Last point had no confidence intervals
# gannet <- rbind(gannet, gannet[.N], gannet[.N])
# gannet[, date := rep(
#   seq.Date(as.Date('2012-09-01'), by = '5 days', length.out = 74),
#   each = 3)]
# 
# setorder(gannet, date, latitude)
# 
# gannet[, value := rep(c('lsd_lat', 'mean_lat', 'usd_lat'), times = 74)]
# gannet[, week := week(date)]
# gannet_agg <- gannet |>
#   dcast(week ~ value, value.var = 'latitude', fun.aggregate = mean) |> 
#   _[, species := 'gannet']




gannet <- fread('data/northern gannet (diving bird study).csv')
region <- fread('data/spiegel_A-1.csv')
gannet <- gannet[region, on = c(`tag-local-identifier` = 'tag_id',
                                `individual-local-identifier` = 'band_number')]
gannet <- gannet[region == 'St. Lawrence']
gannet <- gannet[`argos:lc` %in% c('G', 3, 2, 1)]
gannet[, let(week = week(timestamp),
             year = year(timestamp))]

gannet_agg <- gannet[, .(
  mean_lat = mean(`location-lat`),
  mean_lon = mean(`location-long`),
  species = 'northern gannet'
),
by = c('week', 'year', 'individual-local-identifier')
]

ggplot(data = gannet_agg) +
  geom_point(aes(x = week, y = mean_lat)) +
  facet_wrap(~year)

ggplot(data = gannet_agg) +
  geom_point(aes(x = week, y = mean_lon)) +
  facet_wrap(~year) +
  coord_flip()

gannet_agg <- gannet_agg[year == 2015, 
                         .(species, week, mean_lat, mean_lon)]


## NARW
narw <- fread('data/embargo/davis-et-al-2017_narw-daily-detection-data.csv')
narw <- narw[NARW_PRESENCE == 1]
#correct negative longitude
narw[, RECORDER_LONGITUDE := fifelse(RECORDER_LONGITUDE > 0,
                                     -RECORDER_LONGITUDE,
                                     RECORDER_LONGITUDE)]
narw[, week := week(ANALYSIS_END_ISO)]
narw[, year := year(ANALYSIS_END_ISO)]


ggplot(data = narw) +
  geom_point(aes(x = week, y = RECORDER_LATITUDE)) +
  facet_wrap(~yr)

ggplot(data = narw) +
  geom_point(aes(x = week, y = RECORDER_LONGITUDE)) +
  facet_wrap(~yr) +
  coord_flip()

# narw_agg <- narw[, .(mean_lat = mean(RECORDER_LATITUDE, na.rm = T),
#                      lsd_lat = mean(RECORDER_LATITUDE, na.rm = T) -
#                        sd(RECORDER_LATITUDE, na.rm = T),
#                      usd_lat = mean(RECORDER_LATITUDE, na.rm = T) +
#                        sd(RECORDER_LATITUDE, na.rm = T),
#                      species = 'NARW'),
#                  by = c('week', 'yr')] 

# There are no individual records, so just take the unique locationXweekXyear combos
narw_agg <- unique(narw, by = c('week', 'year',
                                'RECORDER_LATITUDE','RECORDER_LONGITUDE')) |>
  _[, let(mean_lat = RECORDER_LATITUDE,
          mean_lon = RECORDER_LONGITUDE,
          species = 'NARW')]

ggplot(data= narw_agg, aes(x =week)) +
    geom_point(aes(y = mean_lat), size = 0.5) +
    facet_wrap(~year)
ggplot(data= narw_agg, aes(x =week)) +
  geom_point(aes(y = mean_lon), size = 0.5) +
  facet_wrap(~year) +
  coord_flip()

# pick 2013
# ggplot(data= narw_agg, aes(x =week)) +
#   geom_errorbar(aes(ymin = lsd, ymax = usd)) +
#   geom_point(aes(y = mean), size = 0.5) +
#   facet_wrap(~yr)

# narw_agg <- narw_agg[yr == 2013, -'yr']
narw_agg <- narw_agg[year == 2013, 
                     .(species, week, mean_lat, mean_lon)]
#correct missing values
# narw_agg <- narw_agg[data.table(week = 1:52,
#                                 species = 'NARW'), on = c('week', 'species')]

# Merge
agg <- rbind(sb_agg, gannet_agg, narw_agg)

## Export
# fwrite(agg, 'data/joint_flyway_latitudes.csv')

agg[, week_date := as.Date('2024-01-01') + (week-1) * 7]

# Same plot, different colors
lat <- ggplot(data = agg, aes(x = week_date, color = species)) +
  # geom_errorbar(aes(ymin = lsd, ymax = usd)) +
  geom_point(aes(y = mean_lat)) +
  scale_x_date(date_labels = '%b') +
  labs(x = '', y = 'Latitude') +
  theme_minimal()
lon <- ggplot(data = agg, aes(x = week_date, color = species)) +
  # geom_errorbar(aes(ymin = lsd, ymax = usd)) +
  geom_point(aes(y = mean_lon)) +
  scale_x_date(date_labels = '%b') +
  labs(x = '', y = 'Longitude') +
  theme_minimal() +
  coord_flip()

ggplot() +
  geom_density_2d_filled(data = agg, aes(x = week, y = mean_lat)) +
  facet_wrap(~species, scales = 'free')

# Facet by species
ggplot(data = agg, aes(x = week_date)) +
  geom_errorbar(aes(ymin = lsd, ymax = usd)) +
  geom_point(aes(y = mean)) +
  facet_wrap(~species, ncol = 1, scales = 'free_y') +
  scale_x_date(date_labels = '%b') +
  labs(x = '', y = 'Latitude') +
  theme_minimal()


############# Model ----
library(mgcv)
agg[, species := as.factor(species)]
m1 <- gam(mean_lat ~  
            0 + species +
            s(week, by = species, k = 15, bs = 'cc'),
          knots = list(week = c(1,53)),
          data = agg, method = 'REML')

pred_dat <- data.table(
  week = rep(seq(1,53, length.out = 100), times = 3),
  species = rep(c('northern gannet', 'striped bass', 'NARW'),
                each = 100)
)
pd <- predict(m1, newdata = pred_dat, se.fit = TRUE)

pred_dat[, ':='(pred = pd$fit,
                lci = pd$fit - 2*pd$se.fit,
                uci = pd$fit + 2*pd$se.fit)]



smooth_plot <-
  ggplot(data = pred_dat) +
  geom_ribbon(aes(x = week, ymin = lci, ymax = uci),
              fill = 'lightgray') +
  geom_line(aes(x = week, y = pred)) +
  # geom_errorbar(data = agg, aes(x = week, ymin = lsd, ymax = usd),
  #               color = 'orange') +
  geom_point(data = agg, aes(x = week, y = mean_lat),
             color = 'orange', alpha = 0.5, size = 0.5) +
  facet_wrap(~species, ncol = 1, scales = 'free_y') +
  labs(x = 'Week of year', y = 'Latitude') +
  theme_minimal() +
  theme(
    plot.margin = unit(rep(0, 4), 'mm')
  )




library(gratia)

deriv_plot <-
  derivatives(m1, data = pred_dat, 
              eps = 5) |> 
  ggplot() +
  geom_ribbon(aes(x = week, ymin = .lower_ci, ymax = .upper_ci),
              fill = 'lightgray') +
  geom_line(aes(x = week, y = .derivative)) +
  facet_wrap(~species, ncol = 1, scales = 'free_y')+
  labs(x = 'Week of year', y = 'Latitudinal velocity (degrees/week)') +
  theme_minimal() +
  theme(
    plot.margin = unit(rep(0, 4), 'mm')
  )


library(patchwork)
library(svglite)
svglite(
  'ms_figures/figure7.svg',
  width = 170 / 25.4,
  height = 170 / 25.4 / (650/575)
)
smooth_plot + deriv_plot
dev.off()


m2 <- gam(mean_lon ~  
            0 + species +
            s(week, by = species, k = 15, bs = 'cc'),
          knots = list(week = c(1,53)),
          data = agg, method = 'REML')


m3 <- gam(
  list(
    mean_lon ~  
      0 + species +
      s(week, by = species, k = 15, bs = 'cc'),
    mean_lat ~  
      0 + species +
      s(week, by = species, k = 15, bs = 'cc')
  ),
  family = mvn(),
  data = agg
  )


pred_dat <- data.table(
  week = rep(1:53, times = 3),
  species = rep(c('northern gannet', 'striped bass', 'NARW'),
                each = 53)
)
pd <- predict(m3, newdata = pred_dat, se.fit = TRUE)

pred_dat[, ':='(pred_lon = pd$fit[,1],
                lci_lon = pd$fit[,1] - 2*pd$se.fit[,1],
                uci_lon = pd$fit[,1] + 2*pd$se.fit[,1],
                pred_lat = pd$fit[,2],
                lci_lat = pd$fit[,2] - 2*pd$se.fit[,2],
                uci_lat = pd$fit[,2] + 2*pd$se.fit[,2])]

ggplot(data = pred_dat) +
  geom_ribbon(aes(x = week, ymin = lci_lon, ymax = uci_lon),
              fill = 'lightgray') +
  geom_line(aes(x = week, y = pred_lon)) +
  # geom_errorbar(data = agg, aes(x = week, ymin = lsd, ymax = usd),
  #               color = 'orange') +
  geom_point(data = agg, aes(x = week, y = mean_lon),
             color = 'orange', alpha = 0.5) +
  facet_wrap(~species, ncol = 3, scales = 'free_y') +
  labs(x = 'Week of year', y = 'Latitude') +
  theme_minimal() +
  coord_flip()


################# OMERC ----
# translate to a line between GA/FL border and Cape Sable, NS
# (-81.45, 30.71) # GA/FL
# (-65.62, 43.39) # Cape Sable, NS

# Start with all points at constant longitude (-81.45)
# agg[, lon := -81.45]
# agg[, lon := -74.4]

library(sf)

land <- '/vsizip/data/spatial/ne_10m_land.zip/ne_10m_land.shp' |> 
  sf::read_sf() |>
  st_intersection(
    st_bbox(c(xmin = -90, xmax = -50,
              ymin = 20, ymax = 60), crs = 4326) |> 
      st_as_sfc()
  )



# agg[, lon_adj := (tan(atan(15.83/12.68)) * (mean - 30.71)) + -81.45]
# agg_sf <- st_as_sf(agg[!is.na(mean)], coords = c('lon_adj', 'mean'), crs = 4326)
agg_sf <- st_as_sf(agg, coords = c('mean_lon', 'mean_lat'), crs = 4326)
proj_string <- '+proj=omerc +lonc=-81.45 +lat_0=30.71 +gamma=-40'
agg_omerc <- agg_sf |> 
  st_transform(proj_string)
st_geometry(agg_omerc) <- st_geometry(agg_omerc) * matrix(c(0,1))
st_crs(agg_omerc) <- st_crs(proj_string)

land_omerc <- land |> 
  st_transform(proj_string)

ggplot() +
  geom_sf(data = land) + 
  geom_sf(data = agg_sf) +
  facet_wrap(~species)
+
  coord_sf(crs = st_crs(land_omerc))

# l <- land |> 
#   st_make_valid() |> 
#   st_crop(ymin = 30, ymax = 60, xmin = -80, xmax = -50)

## Convert mean lat to linestring spanning the coast (-90 to -50?)
## intersect with land

agg_sf <- copy(agg)

agg_sf[, geom :=
         list(
           list(
             st_linestring(
               matrix(
                 c(-90, -50,rep(as.numeric(.SD[,'mean']), 2)),
                 ncol = 2
               )
             )
           )
         ),
       by = .I]  
agg_sf[, geom := st_sfc(geom, crs = 4326)]

agg_sf <- agg_sf |> 
  st_as_sf()

k <- st_intersection(
  agg_sf,
  l[l$min_zoom == 0,]
)
|> 
  st_nearest_points(
    st_linestring(matrix(c( -45, -45, 0, 90), ncol=2)) |> 
      st_sfc(crs = 4326)
  )

plot(k)

st_nearest_points(agg_sf, l) |>
  plot()

|> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-68 +alpha=-40')

crd <- st_coordinates(agg_sf) |> 
  data.frame()
setDT(crd)

crd <- cbind(agg, crd)

crd[, spd := (Y - shift(Y, 1)) / 7 / 1000, by = species]

ggplot(crd) +
  geom_line(aes(x = week, y = spd, color = species))













############ Speeds ----
setorder(agg, species, week)

# Pointwise lat/week
agg[, spd := (shift(mean, -1) - mean) / 7, by = species]
ggplot(agg) +
  geom_line(aes(x = week, y = spd)) +
  facet_wrap(~species, ncol = 1, scale = 'free_y') +
  scale_x_date(date_labels = '%b') +
  labs(x = NULL, y = 'Latitudinal velocity') +
  theme_minimal()

spd_mo <- copy(agg)
spd_mo[, mon := month(week)]
spd_mo <- spd_mo[, .(mean = mean(mean, na.rm = T)),
                 by = c('mon', 'species')]

spd_mo[, spd := (shift(mean, -1) - mean), by = species]



agg_na <- copy(agg)
agg_na[,lat_before := nafill(mean, "locf")]
# Bring back the next non-missing dist
agg_na[,lat_after := nafill(mean, "nocb")]
# rleid will create groups based on run-lengths of values within the data.
# This means 4 NA's in a row will be grouped together, for example.
# We then count the missings and add 1, because we want the 
# last NA before the next non-missing to be less than the non-missing value.
agg_na[, rle := rleid(mean)][,missings := max(.N +  1 , 2), by = rle] 
agg_na[is.na(mean), ':='(mean = lat_before + .SD[,.I] *
                           (lat_after - lat_before)/(missings)), by = rle]
agg_na[,':='(lat_before = NULL,
             lat_after = NULL,
             rle = NULL,
             missings = NULL,
             type = 'interp')]
agg[,type := '']

agg <- rbind(agg, agg_na)
setorder(agg, species, week)

agg[, spd := (shift(mean, -1) - mean) / 7, by = c('species', 'type')]

ggplot() +
  geom_line(data = agg[type == 'interp'],
            aes(x = week, y = spd),
            linetype = 'dashed') +
  geom_line(data = agg[type == ''],
            aes(x = week, y = spd)) +
  facet_wrap(~species, ncol = 1, scale = 'free_y') +
  scale_x_date(date_labels = '%b') +
  labs(x = NULL, y = 'Latitudinal velocity') +
  theme_minimal()