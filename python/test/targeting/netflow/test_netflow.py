import os
import numpy as np
import matplotlib.pyplot as plt

from test_base import TestBase

import pfs.ga.targeting.netflow as nf
from pfs.ga.targeting.netflow import ScienceTarget, CalibTarget, Telescope
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.Bench import Bench
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID

class NetflowTest(TestBase):
    def test_netflow(self):
        obs = self.load_test_observation()
        cmd, _ = self.get_test_cmd()
        sel = self.get_test_color_selection(cmd)

        mask = sel.apply(obs)
        
        targets = []
        ids = obs.get_id(mask=mask)
        coords = obs.get_coords(mask=mask)
        for id, ra, dec in zip(ids, *coords):
            exp_time = 1200
            priority = 1
            prefix = 'sci'
            s = ScienceTarget(id, ra, dec, exp_time, priority, prefix)

            targets.append(s)

        ra = coords[0].mean()
        dec = coords[1].mean()
        targets.append(CalibTarget(id, ra, dec, 'cal'))
        targets.append(CalibTarget(id, ra, dec, 'sky'))

        bench = Bench(layout="full")

        nvisit = 1

        telescopes = []
        forbiddenPairs = []
        for _ in range(nvisit):
            tra = ra + np.random.normal() * 1e-2
            tdec = dec + np.random.normal() * 1e-2
            posang = 0.
            obs_time = "2016-04-03T08:00:00Z"
            telescope = Telescope(tra, tdec, posang, obs_time)
            telescopes.append(telescope)

            forbiddenPairs.append([])

        # get focal plane positions for all targets and all visits
        fp_pos = [ t.get_fp_positions(targets) for t in telescopes ]

        # create the dictionary containing the costs and constraints for all classes
        # of targets
        classdict = {}
        classdict["sci_P1"] = {"nonObservationCost": 100,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P2"] = {"nonObservationCost": 90,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P3"] = {"nonObservationCost": 80,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P4"] = {"nonObservationCost": 70,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P5"] = {"nonObservationCost": 60,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P6"] = {"nonObservationCost": 50,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sci_P7"] = {"nonObservationCost": 40,
                            "partialObservationCost": 1e9, "calib": False}
        classdict["sky"] = {"numRequired": 1,
                            "nonObservationCost": 1e6, "calib": True}
        classdict["cal"] = {"numRequired": 1,
                            "nonObservationCost": 1e6, "calib": True}
        
        # optional: slightly increase the cost for later observations,
        # to observe as early as possible
        vis_cost = [i*10. for i in range(nvisit)]


        # optional: penalize assignments where the cobra has to move far out
        def cobraMoveCost(dist):
            return 0.1*dist


        # duration of one observation in seconds
        t_obs = 300.

        gurobiOptions = dict(seed=0, presolve=1, method=4, degenmoves=0,
                            heuristics=0.8, mipfocus=0, mipgap=1.0e-04)


        # let's pretend that most targets have already been completely observed,
        # and that the rest has been partially observed
        alreadyObserved={}
        for t in targets:
            alreadyObserved[t.ID] = 3
        for t in targets[::10]:
            alreadyObserved[t.ID] = 1

        forbiddenPairs = []
        for i in range(nvisit):
            forbiddenPairs.append([])

        done = False
        while not done:
            # compute observation strategy
            problem = nf.buildProblem(bench, targets, fp_pos, classdict, t_obs,
                                vis_cost, cobraMoveCost=cobraMoveCost,
                                collision_distance=2., elbow_collisions=True,
                                gurobi=False, gurobiOptions=gurobiOptions,
                                alreadyObserved=alreadyObserved,
                                forbiddenPairs=forbiddenPairs)

            print("solving the problem")
            problem.solve()

            # extract solution
            res = [{} for _ in range(nvisit)]
            for k1, v1 in problem._vardict.items():
                if k1.startswith("Tv_Cv_"):
                    visited = problem.value(v1) > 0
                    if visited:
                        _, _, tidx, cidx, ivis = k1.split("_")
                        res[int(ivis)][int(tidx)] = int(cidx)

            print("Checking for trajectory collisions")
            ncoll = 0
            for ivis, (vis, tp) in enumerate(zip(res, fp_pos)):
                selectedTargets = np.full(len(bench.cobras.centers), NULL_TARGET_POSITION)
                ids = np.full(len(bench.cobras.centers), NULL_TARGET_ID)
                for tidx, cidx in vis.items():
                    selectedTargets[cidx] = tp[tidx]
                    ids[cidx] = ""
                for i in range(selectedTargets.size):
                    if selectedTargets[i] != NULL_TARGET_POSITION:
                        dist = np.abs(selectedTargets[i]-bench.cobras.centers[i])

                simulator = CollisionSimulator(bench, TargetGroup(selectedTargets, ids))
                simulator.run()
                if np.any(simulator.endPointCollisions):
                    print("ERROR: detected end point collision, which should be impossible")
                coll_tidx = []
                for tidx, cidx in vis.items():
                    if simulator.collisions[cidx]:
                        coll_tidx.append(tidx)
                ncoll += len(coll_tidx)
                for i1 in range(0,len(coll_tidx)):
                    for i2 in range(i1+1,len(coll_tidx)):
                        if np.abs(tp[coll_tidx[i1]]-tp[coll_tidx[i2]])<10:
                            forbiddenPairs[ivis].append((coll_tidx[i1],coll_tidx[i2]))

            print("trajectory collisions found:", ncoll)
            done = ncoll == 0
        