# rinexsatposition
Vectorized Python implementation of Keplerian orbit propagator that gives satellite positions from RINEX navigation message files.

Has not been checked extensively for accuracy, but is accurate 'enough' for GPS satellites.

Not recommended for highly precise use.

RINEX nav files are read with Georinex. Requires Pandas and Numpy. Use pymap3d to convert to lat/lon/alt. 

In UTC, the input is:

``test=rinexnav2satxyz('brdc0060.18n.Z','G05',np.array(['2018-01-06T02:30:43.000','2018-01-06T10:00:00.000']))``

and get coordinates with 
``pymap3d.ecef.ecef2geodetic(test[:,0],test[:,1],test[:,2])``

Note that the output will be a lists of lats, lons, and altitudes corresponding to the given timestamps.
