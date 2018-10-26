

# Credits to:

# Grzegorz Sionkowski:
# Bayesian Optimization - https://www.kaggle.com/sionek/bayesian-optimization

# Heng CherKeng
# Track Extension - https://www.kaggle.com/c/trackml-particle-identification/discussion/58194

# Improvement ideas:
#   - Implement an efficient merging function.
#   - Try several architectures of merging / extending tracks (based on size and quality criterion).
#   - Add a TrackQuality function to sort and merge more efficiently.
#   - Use a garbage track (0) in trackML.
#   - Compute more precise spatial position using cells data.



rm(list=ls())
cat('\014')



#------------------------------------------------------------------------------------------------
#-------------------- 0. MODULES ----------------------------------------------------------------


# Data hangling
library(data.table)
library(bit64)

# Clustering
library(dbscan)

# Nearest neighbors search
library(RANN)

# Parallel environment
library(doParallel)

# Hyperparameter tuning
library(rBayesianOptimization)

path_to_train = '/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/train_100_events/'
path_to_test = '/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/test/'
path_to_submit = '/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/'


#------------------------------------------------------------------------------------------------
#-------------------- 1. SCORE ------------------------------------------------------------------


score <- function(predictions, expected) {
  
  data = merge(predictions, expected[, .(hit_id, particle_id, weight)])
  
  data[, Np := .N, by = particle_id]                         # Np  = Hits per Particle
  data[, Nt := .N, by = track_id]                            # Nt  = Hits per Track
  data[, Ntp := .N, by = list(track_id, particle_id)]        # Ntp = Hits per Particle per Track
  
  data[, r1 := Ntp/Nt]
  data[, r2 := Ntp/Np]
  
  sum(data[r1>.5 & r2>.5, weight])
}


#------------------------------------------------------------------------------------------------
#-------------------- 2. TRACKING FUNCTION ------------------------------------------------------


trackML <- function(data, w1, w2, w3, coef, eps, iterNumber) {
  
  # Hough transformation
  data[, track := hit_id]
  data[, N1    := 1L]
  data[, r     := sqrt(x*x + y*y + z*z)]
  data[, rt    := sqrt(x*x + y*y)]
  data[, a0    := atan2(y, x)]             # (z, r) ?
  data[, z1    := z/rt]
  data[, z2    := z/r]
  
  mm = 1
  
  for (i in 0:iterNumber) {
    
    mm = mm * (-1)
    
    data[, a1    := a0+mm*(rt+coef*rt^2)/1000*(i/2)/180*pi]
    data[, sina1 := sin(a1)]
    data[, cosa1 := cos(a1)]
    
    data_scaled = scale(data[, .(sina1,cosa1,z1,z2)])
    cx = c(w1,w1,w2,w3)
    
    for (j in 1:ncol(data_scaled)) { data_scaled[,j] = data_scaled[,j] * cx[j] }
    
    # clustering
    cluster = dbscan(data_scaled, eps, minPts=1)
    
    data[, new_track := cluster$cluster]
    data[, N2 := .N, by=new_track]
    
    max_track = max(data$track)
    
    data[, track := ifelse( (N2 > N1) & (N2 < 20), new_track + max_track, track)]
    data[, track := as.integer(as.factor(track))]
    data[, N1 := .N, by=track]
  }
  
  return(data$track)
}


eventNumber = 1018

hits        = fread(paste0(path_to_train,'event00000',eventNumber,'-hits.csv'))

# Clustering
track = trackML(hits, 1, 1, 1, 0.003, 0.000004, 50L)

predictions = data.table(event_id = eventNumber,
                         hit_id   = hits$hit_id,
                         track_id = track)

limit = 0.04
numNeighbours = 15

### TEST END

# Merging
setkey(predictions, hit_id)
setkey(hits, hit_id)

data = hits[predictions]

# Hough transform
data[, `:=` (d       = sqrt(x*x + y*y + z*z),
             r       = sqrt(x*x + y*y))]

data[, arctan2 := atan2(z, r)]          # (x, y) ?

# Checking all angles
angle = 1
  
inter = data[(arctan2 > pi*(angle-1.5)/180) & (arctan2 < pi*(angle+1.5)/180)]

minNumNeighbours = nrow(inter)
  
hit_ids = inter[, hit_id]

inter[, `:=` (r = r/1000,
              a = atan2(y, x))]

inter[, `:=` (c = cos(a),
              s = sin(a))]

# # Exact search
# tree = nn(inter[, .(c, s, r)], p=min(minNumNeighbours, numNeighbours))

# Approximate search
tree = nn2(inter[, .(c, s, r)], k=min(minNumNeighbours, numNeighbours), treetype='kd')

tracks = inter[, unique(track_id)]
numTracks = length(tracks)
minLength = 3

# Extension loop
for (track in tracks) { 

idx = inter[track_id == track, which = TRUE]
  
if (sign(angle) == 1) {sign = 1} else {sign = -1}
idx = idx[order(inter[idx, sign * z])]

## Starting & ending direction
da0 = inter[idx[2], a] - inter[idx[1], a]
dr0 = inter[idx[2], r] - inter[idx[1], r]
direction0 = atan2(dr0, da0)

da1 = inter[idx[length(idx)], a] - inter[idx[length(idx)-1], a]
dr1 = inter[idx[length(idx)], r] - inter[idx[length(idx)-1], r]
direction1 = atan2(dr1, da1)

## Extend starting point
idxNeighbours = tree$nn.idx[1, ]
direction = atan2(inter[idx[1], r] - inter[idxNeighbours, r], inter[idx[1], a] - inter[idxNeighbours, a])
diff = 1 - cos(direction - direction0)

idxNeighbours = idxNeighbours[(inter[idx[1], r] - inter[idxNeighbours, r] > 0.01) & (diff < 1 - cos(limit))]
data[hit_id %in% hit_ids[idxNeighbours], track_id := track]

## Extend ending point
idxNeighbours = tree$nn.idx[length(idx), ]
direction = atan2(inter[idxNeighbours, r] - inter[idx[length(idx)], r], inter[idxNeighbours, a] - inter[idx[length(idx)], a])
diff = 1 - cos(direction - direction1)

idxNeighbours = idxNeighbours[(inter[idxNeighbours, r] - inter[idx[length(idx)], r] > 0.01) & (diff < 1 - cos(limit))]
data[hit_id %in% hit_ids[idxNeighbours], track_id := track]

}

debug = data[, .(event_id, hit_id, track_id)]
