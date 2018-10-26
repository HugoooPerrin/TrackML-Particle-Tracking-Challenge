

# Credits to:

# Grzegorz Sionkowski:
# Bayesian Optimization - https://www.kaggle.com/sionek/bayesian-optimization

# Heng CherKeng:
# Track Extension - https://www.kaggle.com/c/trackml-particle-identification/discussion/58194

# Improvement ideas:
#   - Implement an efficient merging function.
#   - Try several architectures of merging / extending tracks (based on size and quality criterion).
#   - Add a TrackQuality function to extend and merge more efficiently.
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

  data[, Np  := .N, by = particle_id]                        # Np  = Hits per Particle
  data[, Nt  := .N, by = track_id]                           # Nt  = Hits per Track
  data[, Ntp := .N, by = list(track_id, particle_id)]        # Ntp = Hits per Particle per Track

  data[, r1 := Ntp/Nt]
  data[, r2 := Ntp/Np]

  sum(data[r1 > .5 & r2 > .5, weight])
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


#------------------------------------------------------------------------------------------------
#-------------------- 3. TRACK EXTENSION --------------------------------------------------------


trackExtension <- function(predictions, hits, limit, numNeighbours, iter) {

  extend <- function(predictions, hits, limit, numNeighbours) {
    
    # Merging
    setkey(predictions, hit_id)
    setkey(hits, hit_id)
    
    data = hits[predictions]
    
    # Extending long tracks first
    data[, count := .N, by = track_id]
    data = data[order(-count)]
    
    trackID = data[, track_id]
    hit_id = data[, hit_id]
    
    rm(predictions, hits)
    
    # Hough transform
    data[, `:=` (d = sqrt(x*x + y*y + z*z),
                 r = sqrt(x*x + y*y))]
    
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
        
        # When it comes to indexing R is optimized for vectors/arrays, and not dataframes/datatables
        r = inter[, r/1000]
        a = inter[, atan2(y, x)]
        c = cos(a)
        s = sin(a)
        
        # Approximate search
        tree = nn2(cbind(c, s, r), k=min(minNumNeighbours, numNeighbours)+1, treetype='kd')$nn.idx[, -c(1)]
        
        tracks = unique(trackInter)
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
          idxNeighbours = tree[last,]
          direction = atan2(r[idxNeighbours] - r[last], a[idxNeighbours] - a[last])
          diff = 1 - cos(direction - direction1)
          
          idxNeighbours = idxNeighbours[(r[idxNeighbours] - r[last] > 0.01) & (diff < 1 - cos(limit))]
          trackID[hit_id %in% hit_ids[idxNeighbours]] = track
        }
      }
    }
    # Back into datatable !
    data[, track_id := trackID]
    return(data[, .(event_id, hit_id, track_id)])
  }
  
  for (i in 1:iter) { predictions = extend(predictions, hits, limit, numNeighbours) }
  
  return(predictions)
  
}


#------------------------------------------------------------------------------------------------
#-------------------- 3. TRACK MERGING ----------------------------------------------------------


# Not working ! 

trackMerging <- function(predictions, hits, limit, numNeighbours, iter) {
  
  merge <- function(predictions, hits, limit, numNeighbours) {
    
    # Merging
    setkey(predictions, hit_id)
    setkey(hits, hit_id)
    
    data = hits[predictions]
    
    # Merging long tracks first (only ?)
    data[, count := .N, by = track_id]
    data = data[order(-count)]
    
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

      minNumNeighbours = inter[, .N]
      condition1 = (minNumNeighbours >= 3)
      
      if (condition1) {
        
        # Computing merging metrics
        inter = inter[order(-z*sign)]
        inter[, isFirst := !duplicated(track_id)]
        inter = inter[order(-z*sign)]
        inter[, isLast := !duplicated(track_id)]
        
        # Careful with the order of index
        hit_ids = inter[, hit_id]
        track_ids = inter[, track_id]
        isFirst = inter[, isFirst]
        isLast = inter[, isLast]
          
        z = inter[, z]
        
        r = inter[, r/1000]
        a = inter[, atan2(y, x)]
        c = cos(a)
        s = sin(a)
        
        # Approximate search
        tree = nn2(cbind(c, s, r), k=min(minNumNeighbours, numNeighbours)+1, treetype='kd')$nn.idx[, -c(1)]
        
        tracks = unique(track_ids)
        condition = inter[, .(.N>=3 & .N<10), by = track_id][, V1]
        
        # merging loop
        for (track in tracks[condition]) {
          
          idx = which(track_ids == track)
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
          
          ## Merge starting point
          idxNeighbours = tree[first,]
          direction = atan2(r[first] - r[idxNeighbours], a[first] - a[idxNeighbours])
          diff = 1 - cos(direction - direction0)
          
          idxNeighbours = idxNeighbours[(r[idx[1]] - r[idxNeighbours] > 0.01) & (diff < 1 - cos(limit))]
          idxNeighbours = idxNeighbours[isLast[idxNeighbours]][1]
          trackID[trackID %in% track_ids[idxNeighbours]] = track
          
          ## Extend ending point
          idxNeighbours = tree[last,]
          direction = atan2(r[idxNeighbours] - r[last], a[idxNeighbours] - a[last])
          diff = 1 - cos(direction - direction1)
          
          idxNeighbours = idxNeighbours[(r[idxNeighbours] - r[last] > 0.01) & (diff < 1 - cos(limit))]
          idxNeighbours = idxNeighbours[isFirst[idxNeighbours]][1]
          trackID[trackID %in% track_ids[idxNeighbours]] = track
        }
      }
    }
    data[, track_id := trackID]
    return(data[, .(event_id, hit_id, track_id)])
  }
  
  for (i in 1:iter) { predictions = merge(predictions, hits, limit, numNeighbours) }
  
  return(predictions)
}


#------------------------------------------------------------------------------------------------
#-------------------- 4. BAYESIAN OPTIMIZATION FOR HYPERPARAMETER TUNING ------------------------


optimize = FALSE

if (optimize) {
  
    print("Optimizing")

    bigtime = Sys.time()

    registerDoParallel(cores=8)
    print("All cores initialized")

    Opt = foreach(eventNumber = 1009:1024, .combine = rbind)  %dopar%  {                                  # Optimization using 16 events
      
      hits        = fread(paste0(path_to_train,'event00000',eventNumber,'-hits.csv'))
      expected    = fread(paste0(path_to_train,'event00000',eventNumber,'-truth.csv'), stringsAsFactors = T)

      BayesianOptFun <- function(w1, w2, w3, coef, eps, iterNumber) {
        
        track = trackML(hits, w1, w2, w3, coef, eps, iterNumber)
        predictions = data.table(event_id = eventNumber,
                                 hit_id   = hits$hit_id,
                                 track_id = track)

        score = score(predictions, expected)
        return(list(Score=score, Pred=0))
      }
      
      tempOpt = BayesianOptimization(BayesianOptFun,
                                     
                                     bounds = list(w1            = c(1.3, 2),
                                                   w2            = c(0.45, 0.65),
                                                   w3            = c(0.2, 0.4), 
                                                   coef          = c(0.000004, 0.00001),
                                                   eps           = c(0.003, 0.005), 
                                                   iterNumber    = c(195L, 245L)),
                                     
                                     init_points = 4, n_iter = 40,
                                     acq = "ucb", kappa = 2.576,
                                     verbose = FALSE)
      
      subOpt = data.table(w1            = tempOpt$Best_Par[[1]],
                          w2            = tempOpt$Best_Par[[2]],
                          w3            = tempOpt$Best_Par[[3]],
                          coef          = tempOpt$Best_Par[[4]],
                          eps           = tempOpt$Best_Par[[5]],
                          iterNumber    = tempOpt$Best_Par[[6]])
      
      return(subOpt)
    }

    stopImplicitCluster()

    fwrite(Opt, paste0(path_to_submit, 'BestParameters2.csv'))
    cat("Bayesian Optimization done in:", round(difftime(Sys.time(), bigtime,units = c("min")), digits = 2), "mins")   # 229.41 mins (1st) / X mins (2nd)

} else {
    
    print("Loading previous parameters")
    Opt = fread(paste0(path_to_submit, 'BestParameters2.csv'))
}


#------------------------------------------------------------------------------------------------
#-------------------- 5. MERGING/EXTENSION ARCHITECTURE TESTING ---------------------------------


# Opening
eventNumber = 1052
hits = fread(paste0(path_to_train,'event00000',eventNumber,'-hits.csv'))
truth = fread(paste0(path_to_train,'event00000',eventNumber,'-truth.csv'), stringsAsFactors = T)

## Optimized parameters
w1             = Opt[, mean(w1)]
w2             = Opt[, mean(w2)]
w3             = Opt[, mean(w3)]
coef           = Opt[, mean(coef)]
eps            = Opt[, mean(eps)]
iterNumber     = Opt[, round(mean(iterNumber), digits = 0)]
limitE         = 0.04
numNeighboursE = 17
iterExtend     = 8
limitM         = 0.04
numNeighboursM = 17
iterMerge      = 1

# Clustering
time=Sys.time()
track = trackML(hits, w1, w2, w3, coef, eps, iterNumber)
sub = data.table(event_id = eventNumber,
                 hit_id   = hits$hit_id,
                 track_id = track)
time2=Sys.time()
cat('Clustering score', score(sub, truth),'in',round(difftime(time2, time, units=c("sec")), digits=0), 'secs\n')

# Extension
time=Sys.time()
sub = trackExtension(sub, hits, limitE, numNeighboursE, iterExtend)
time2=Sys.time()
cat('Extension score', score(sub, truth),'in',round(difftime(time2, time, units=c("sec")), digits=0), 'secs\n')

# Merging
time=Sys.time()
sub = trackMerging(sub, hits, limitM, numNeighboursM, iterMerge)
time2=Sys.time()
cat('Merging score', score(sub, truth),'in',round(difftime(time2, time, units=c("sec")), digits=0), 'secs\n')


#------------------------------------------------------------------------------------------------
#-------------------- 6. COMPUTING PREDICTIONS --------------------------------------------------

predict = FALSE

if (predict) {
  
  bigtime = Sys.time()
  
  names = list.files(path_to_test, pattern = 'hits')
  
  ## Optimized parameters
  w1            = Opt[, mean(w1)]
  w2            = Opt[, mean(w2)]
  w3            = Opt[, mean(w3)]
  coef          = Opt[, mean(coef)]
  eps           = Opt[, mean(eps)]
  iterNumber    = Opt[, round(mean(iterNumber), digits = 0)]
  limit         = 0.04
  numNeighbours = 17
  iterExtend    = 8
  
  registerDoParallel(cores=8)
  print("All cores initialized")
  
  ## Predicting
  sub = foreach(eventNumber = 0:124, .combine = rbind)  %dopar%  {
    
    hits = fread(paste0(path_to_test, names[eventNumber+1]))
    
    track = trackML(hits, w1, w2, w3, coef, eps, iterNumber)
    
    subEvent = data.table(event_id = eventNumber,
                          hit_id   = hits$hit_id,
                          track_id = track)
    
    subEvent = trackExtension(subEvent, hits, limit, numNeighbours, iterExtend)
    
    rm(hits)
    
    return(subEvent)
  }
  
  stopImplicitCluster()
  
  ## Saving
  sub = sub[order(event_id, hit_id)]
  fwrite(sub, paste0(path_to_submit, '4th_submission.csv'))
  
  cat("Computing final predictions done in:", round(difftime(Sys.time(), bigtime,units = c("min")), digits = 2), "mins")  
  # 28 mins (1st) / 174 mins (2nd) / 275 mins (3rd) / X (4th)
}
