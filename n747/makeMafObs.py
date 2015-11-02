from lsst.sims.maf.objUtils import runMoObs

runMoObs('neo_S0.des', 'neo_S0_allObs.txt', 'enigma_1189_sqlite.db', sqlconstraint='night=747', tstep=1.5/24.)


