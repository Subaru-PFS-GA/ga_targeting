import numpy as np

from ics.cobraCharmer.cobraCoach import engineer
from ics.cobraCharmer.pfi import PFI
from ics.cobraCharmer.pfiDesign import PFIDesign

from .setup_logger import logger

# TODO: this should be removed

class CollisionSimulator():
    """
    Simulate trajectory collisions for cobras. This is a superset of the features
    implemented in ics.cobraOps.CollisionSimulator and ics.cobraOps.CollisionSimulator2.
    """

    def __init__(self, 
                 instrument,
                 targets,
                 trajectory_steps=None,
                 trajectory_step_width=None):
        
        """
        Initializes the collisions simulator.

        Targets must be ordered by cobra ID.
        """
        
        self.__instrument = instrument
        self.__targets = targets
        self.__trajectory_steps = trajectory_steps if trajectory_steps is not None else 200
        self.__trajectory_step_width = trajectory_step_width if trajectory_step_width is not None else 50

        # Check which cobras are assigned to a target
        self.__assigned_cobras = self.__targets.notNull.copy()
        self.__good_cobras = ~self.__instrument.bench.cobras.hasProblem

        # Define some internal variables that will filled by the run method
        self.__final_fiber_positions = None
        self.__trajectories = None
        self.__fiber_positions = None
        self.__elbow_positions = None
        self.__steps_count = None
        self.__association_trajectory_collisions = None
        self.__association_endpoint_collisions = None
        self.__trajectory_collisions = None
        self.__endpoint_collisions = None

    def __get_trajectory_collisions(self):
        return self.__trajectory_collisions
    
    trajectory_collisions = property(__get_trajectory_collisions)

    def __get_endpoint_collisions(self):
        return self.__endpoint_collisions
    
    endpoint_collisions = property(__get_endpoint_collisions)

    def run(self, solveCollisions=True):
        """Runs the collisions simulator.

        Parameters
        ----------
        solveCollisions: bool, optional
            If True, the simulator will try to solve trajectory collisions,
            changing the movement directions and strategies of the affected
            cobras. Default is True.

        """
        # Calculate the final fiber positions
        self.__calculate_final_fiber_positions()
    
        # Optimize the unassigned cobra positions to minimize their possible
        # collisions with other cobras
        # self.__optimize_unassigned_cobra_positions()
        self.__send_unassigned_cobras_to_home()

        # # Define the theta and phi movement directions
        # self.defineMovementDirections()

        # # Define the theta and phi movement strategies
        # self.defineMovementStrategies()

        # Calculate the cobra trajectories
        self.__calculate_trajectories()

        # # Detect cobra collisions during the trajectory
        self.__detect_trajectory_collisions()

        # # Solve trajectory collisions if requested
        # if solveCollisions:
        #     self.solveTrajectoryCollisions(True)
        #     self.solveTrajectoryCollisions(False)
        #     self.solveTrajectoryCollisions(True)
        #     self.solveTrajectoryCollisions(False)

    def __calculate_final_fiber_positions(self):
        """
        Calculates the cobras final fiber positions.
        """

        # Make sure that broken cobras are in the home position

        # Set the final fiber positions to their associated target positions,
        # leaving unassigned cobras at their home positive positions
        self.__final_fiber_positions = self.__instrument.bench.cobras.home0.copy()
        self.__final_fiber_positions[self.__assigned_cobras] = self.__targets.positions[self.__assigned_cobras]

    def __send_unassigned_cobras_to_home(self):
        """
        Sends the unassigned cobras to their home positions.
        """

        # Set the final fiber positions to their home positions
        self.__final_fiber_positions[~self.__assigned_cobras] = self.__instrument.bench.cobras.home0[~self.__assigned_cobras]

    def __optimize_unassigned_cobra_positions(self):
        """
        Finds the unassigned cobras final fiber positions that minimize
        their collisions with other cobras.
        """
        # Calculate the cobras elbow positions at their current final fiber
        # positions
        elbow_positions = self.__instrument.bench.cobras.calculateElbowPositions(self.__final_fiber_positions)

        # Get indexes of unassigned cobras that can be moved
        broken_cobras = (self.__instrument.calib_model.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0
        (unassigned_cobra_idx,) = np.where(~self.__assigned_cobras & ~broken_cobras)

        # Find the optimal position for each unassigned cobra
        for cidx in unassigned_cobra_idx:
            # Get the cobra nearest neighbors
            ncidx = self.__instrument.bench.getCobraNeighbors(cidx)

            # Select only those that are assigned to a target
            ncidx = ncidx[self.__assigned_cobras[ncidx]]

            # Jump to the next cobra if all the neighbors are unassigned
            if len(ncidx) == 0:
                continue

            # Calculate all the possible cobra elbow rotations of the unassigned cobra
            rotation_angles = np.arange(0, 2 * np.pi, 0.1 * np.pi)
            center = self.__instrument.bench.cobras.centers[cidx]
            rotated_elbow_positions = center + (elbow_positions[cidx] - center) * np.exp(1j * rotation_angles)

            # Obtain the angle that maximizes the closer distance to a neighbor
            distances = np.abs(self.__final_fiber_positions[ncidx] - rotated_elbow_positions[:, np.newaxis])
            min_distances = np.min(distances, axis=1)
            optimal_angle = rotation_angles[np.argmax(min_distances)]

            # Update the cobra final fiber position
            self.__final_fiber_positions[cidx] = (self.__final_fiber_positions[cidx] - center) * np.exp(1j * optimal_angle) + center

    def __calculate_trajectories(self, time_step=10, max_steps=2000):
        """Calculates the cobra trajectories.

        Parameters
        ----------
        timeStep: int
            The trajectories time step resolution in steps.
        maxSteps: int
            The trajectories maximum number of steps.

        """

        # Calculate the theta and phi angle ranges based on the calibration product
        # cidx = np.where(self.__good_cobras)
        cidx = np.arange(self.__instrument.bench.cobras.nCobras)

        calib_model = self.__instrument.cobra_coach.pfi.calibModel

        # Calculate the final theta and phi angles for the good cobras
        theta, phi, flags = self.__instrument.cobra_coach.pfi.positionsToAngles(
            self.__instrument.cobra_coach.allCobras[cidx],
            self.__final_fiber_positions[cidx])

        # TODO: positionsToAngles returns two solutions, we should select the best one

        angles_ok = (flags & PFI.SOLUTION_OK) != 0
        if ~np.all(angles_ok[:, 0]):
            raise ValueError("Some cobras have no solution for the final fiber positions")

        # Select the first angles solution
        theta = theta[:, 0]
        phi = phi[:, 0]

        # Initialize the engineer module
        engineer.setCobraCoach(self.__instrument.cobra_coach)
        engineer.setConstantOntimeMode(maxSteps=max_steps)

        # Calculate the cobra trajectories
        # NOTE: this function will report out-of-range error for broken cobras.
        #       this seems to be normal behavior, but we should check it.
        #       cobra status is in self.__instrument.calib_model.status which provides a
        #       status for each cobra.
        self.__trajectories, _ = engineer.createTrajectory(
            cidx, theta, phi,
            tries=8, twoSteps=True, threshold=20.0, timeStep=time_step)

        # Calculate the fiber and elbow positions along the cobra trajectories
        self.__fiber_positions = self.__trajectories.calculateFiberPositions(self.__instrument.cobra_coach)
        self.__elbow_positions = self.__trajectories.calculateElbowPositions(self.__instrument.cobra_coach)
        self.__steps_count = self.__fiber_positions.shape[1]

    def __detect_trajectory_collisions(self):
        """
        Detects collisions in the cobra trajectories.
        """

        # Extract some useful information from the bench instance
        cobra_associations = self.__instrument.bench.cobraAssociations
        link_radius = self.__instrument.bench.cobras.linkRadius

        # Calculate the distances between the cobras links for each step in the
        # trajectory
        start_points1 = self.__fiber_positions[cobra_associations[0]].ravel()
        end_points1 = self.__elbow_positions[cobra_associations[0]].ravel()
        start_points2 = self.__fiber_positions[cobra_associations[1]].ravel()
        end_points2 = self.__elbow_positions[cobra_associations[1]].ravel()
        distances = self.__instrument.bench.distancesBetweenLineSegments(
            start_points1, end_points1, start_points2, end_points2)

        # Reshape the distances array
        distances = distances.reshape((len(cobra_associations[0]), self.__steps_count))

        # Detect trajectory collisions between cobra associations
        minimum_separation = link_radius[cobra_associations[0]] + link_radius[cobra_associations[1]]
        trajectory_collisions = distances < minimum_separation[:, np.newaxis]

        # Check which cobra associations are affected by collisions
        self.__association_trajectory_collisions = np.any(trajectory_collisions, axis=1)
        self.__association_endpoint_collisions = trajectory_collisions[:, -1]

        self.__trajectory_collisions = self.__instrument.bench.cobraAssociations[:, self.__association_trajectory_collisions]
        self.__endpoint_collisions = self.__instrument.bench.cobraAssociations[:, self.__association_endpoint_collisions]
        