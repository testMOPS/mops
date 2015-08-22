# mops
Code related to MOPS test

In our multi-pronged attack to understand the capabilities of LSST in detecting NEOs, a major angle of attack is related to running MOPS on a large set of moving objects (full density solar system model). The primary goal here is to determine: can MOPS run and produce tracks from a full-density + noise input? can we unambiguously separate real and false tracks (or what is the contamination rate)? Requirements:
- complete model of the solar system - this can be lower fidelity than the NEO model used in the cadence study, but must be much larger
- create detections from this model (for all objects) over part of an opsim run
  + start with one day, as a verification step
  + increase to one month (pick a month with a high number of potential tracks, as this is hardest for MOPS)
  + increase to one year
  + doesn't need to include lightcurves, camera footprint (at least at first)
- evaluate the output from this simulation
  + compute resources required to run MOPS
  + true/false ratio of orbits
  + found/missed (but could potentially have been found) ratio of objects (per discovery opportunity?)
  + found/missed ratio as a function of orbital parameters
