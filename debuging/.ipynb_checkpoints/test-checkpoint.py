
from sklearn.neighbors import KDTree
import pandas as pd
import numpy as np

submission = pd.read_csv('/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/test_submission_BO.csv')
hits = pd.read_csv('/home/hugoperrin/Bureau/Datasets/TrackML Particle Tracking Challenge/test/event000000000-hits.csv')

df = submission.merge(hits,  on=['hit_id'], how='left')
df = df.assign(d = np.sqrt( df.x**2 + df.y**2 + df.z**2 ))
df = df.assign(r = np.sqrt( df.x**2 + df.y**2))
df = df.assign(arctan2 = np.arctan2(df.z, df.r))

angle = 10
num_neighbours = 18
limit = 0.04

#df1 = df.loc[(df.arctan2>(angle-0.5)/180*np.pi) & (df.arctan2<(angle+0.5)/180*np.pi)]
df1 = df.loc[(df.arctan2>(angle-1.5)/180*np.pi) & (df.arctan2<(angle+1.5)/180*np.pi)]

min_num_neighbours = len(df1)

hit_ids = df1.hit_id.values
x,y,z = df1[['x', 'y', 'z']].values.T
r  = (x**2 + y**2)**0.5
r  = r/1000
a  = np.arctan2(y,x)
c = np.cos(a)
s = np.sin(a)
#tree = KDTree(np.column_stack([a,r]), metric='euclidean')
tree = KDTree(np.column_stack([c, s, r]), metric='euclidean')


track_ids = list(df1.track_id.unique())
num_track_ids = len(track_ids)
min_length=3

for track in track_ids:
    p = track
    idx = np.where(df1.track_id==p)[0]
    if len(idx) >= 3: break

idx
p

if angle>0:
    idx = idx[np.argsort( z[idx])]
else:
    idx = idx[np.argsort(-z[idx])]


## start and end points  ##
idx0,idx1 = idx[0],idx[-1]
a0 = a[idx0]
a1 = a[idx1]
r0 = r[idx0]
r1 = r[idx1]
c0 = c[idx0]
c1 = c[idx1]
s0 = s[idx0]
s1 = s[idx1] 

da0 = a[idx[1]] - a[idx[0]]  #direction
dr0 = r[idx[1]] - r[idx[0]]
direction0 = np.arctan2(dr0,da0)

da1 = a[idx[-1]] - a[idx[-2]]
dr1 = r[idx[-1]] - r[idx[-2]]
direction1 = np.arctan2(dr1,da1)

## extend start point
ns = tree.query([[c0, s0, r0]], k=min(num_neighbours, min_num_neighbours), return_distance=False)
ns = np.concatenate(ns)

direction = np.arctan2(r0 - r[ns], a0 - a[ns])
diff = 1 - np.cos(direction - direction0)
ns = ns[(r0 - r[ns] > 0.01) & (diff < (1 - np.cos(limit)))]
for n in ns: df.loc[df.hit_id == hit_ids[n], 'track_id'] = p

## extend end point
ns = tree.query([[c1, s1, r1]], k=min(num_neighbours, min_num_neighbours), return_distance=False)
ns = np.concatenate(ns)

direction = np.arctan2(r[ns] - r1, a[ns] - a1)
diff = 1 - np.cos(direction - direction1)
ns = ns[(r[ns] - r1 > 0.01) & (diff < (1 - np.cos(limit)))]
for n in ns:  df.loc[df.hit_id == hit_ids[n], 'track_id'] = p

#print ('\r')
df[['event_id', 'hit_id', 'track_id']]