


rm(list=ls())



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

names = list.files(path_to_test, pattern = 'hits')
predictions = fread('/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/test_submission_BO.csv')
hits = fread(paste0(path_to_test,'/', names[1]))

# angle = 10
# numNeighbours = 18
# limit = 0.04
# track = 6220

#------------------------------------------------------------------------------------------------
#-------------------- 1. EXTENSION --------------------------------------------------------------


extend <- function(predictions, hits, limit, numNeighbours) {
  
  # Merging
  setkey(predictions, hit_id)
  setkey(hits, hit_id)
  
  data = hits[predictions]
  trackID = data[, track_id]
  hit_id = data[, hit_id]
  
  rm(predictions, hits)
  
  # Hough transform
  data[, `:=` (d       = sqrt(x*x + y*y + z*z),
               r       = sqrt(x*x + y*y))]
  
  data[, arctan2 := atan2(z, r)]          # (x, y) ?
  
  # Checking all angles
  for (angle in -89:89) {
    
    if (sign(angle) == 1) {sign = 1} else {sign = -1}
    
    inter = data[(arctan2 > pi*(angle-1.5)/180) & (arctan2 < pi*(angle+1.5)/180)]
    trackInter = inter[, track_id]
    
    minNumNeighbours = inter[, .N]
    condition1 = (minNumNeighbours >= 3)
    
    if (condition1) {
      
      hit_ids = inter[, hit_id]
      z = inter[, z]
      
      r = inter[, r/1000]
      a = inter[, atan2(y, x)]
      c = cos(a)
      s = sin(a)
      
      # Approximate search
      tree = nn2(cbind(c, s, r), k=min(minNumNeighbours, numNeighbours)+1, treetype='kd')$nn.idx[, -c(1)]

      tracks = unique(trackInter)
      minLength = 3
      condition = inter[, .(.N>=3), by = track_id][, V1]
      
      rm(inter)

      # Extension loop
      for (track in tracks[condition]) {
        
        idx = which(trackInter == track)
        idx = idx[order(sign * z[idx])]
        len_idx = length(idx)
        
        ## Starting & ending points
        first = idx[1]
        last = idx[len_idx]
        
        ## Starting & ending direction
        da_first = a[idx[2]] - a[first]
        dr_first = r[idx[2]] - r[first]
        direction0 = atan2(dr_first, da_first)
        
        da_last = a[last] - a[idx[len_idx-1]]
        dr_last = r[last] - r[idx[len_idx-1]]
        direction1 = atan2(dr_last, da_last)
        
        ## Extend starting point
        idxNeighbours = tree[first,]
        direction = atan2(r[first] - r[idxNeighbours], a[first] - a[idxNeighbours])
        diff = 1 - cos(direction - direction0)
        
        idxNeighbours = idxNeighbours[(r[idx[1]] - r[idxNeighbours] > 0.01) & (diff < 1 - cos(limit))]
        trackID[hit_id %in% hit_ids[idxNeighbours]] = track
        
        ## Extend ending point
        idxNeighbours = tree[idx[len_idx],]
        direction = atan2(r[idxNeighbours] - r[last], a[idxNeighbours] - a[last])
        diff = 1 - cos(direction - direction1)
        
        idxNeighbours = idxNeighbours[(r[idxNeighbours] - r[last] > 0.01) & (diff < 1 - cos(limit))]
        trackID[hit_id %in% hit_ids[idxNeighbours]] = track
      }
    }
  }
  data[, track_id := trackID]
  return(data[, .(event_id, hit_id, track_id)])
}

###### TESTING PERFORMANCE ######

library(compiler)
extendCompiled = cmpfun(extend)

bigtime = Sys.time()

predictions_fast = extend(predictions, hits, 0.04, 18)
# predictions_fast = extendCompiled(predictions, hits, 0.04, 18)
# predictions_fast2 = extend2(predictions, hits, 0.04, 18)

cat("Computing extension time:", round(difftime(Sys.time(), bigtime,units = c("sec")), digits = 2), "sec")


py = fread('/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/pythonExtension.csv')

check = data.table(pred=predictions[, track_id],
                   pred_fast=predictions_fast[, track_id],
                   # pred_fast2=predictions_fast2[, track_id],
                   pred_python=py[, track_id])

check[, norm_fast := (pred == pred_fast)]
check[, norm_python := (pred_python == pred)]

check[, fast_python := (pred_python == pred_fast)]
# check[, pred_fast2 := (pred_fast2 == pred_fast)]



