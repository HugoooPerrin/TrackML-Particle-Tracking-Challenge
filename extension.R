

library(data.table)
library(RANN)


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
    data[, track_id := trackID]
    return(data[, .(event_id, hit_id, track_id)])
  }
  
  for (i in 1:iter) { predictions = extend(predictions, hits, limit, numNeighbours) }
  
  return(predictions)
  
}