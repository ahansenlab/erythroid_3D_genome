# The following code was adapted from the simulation code used in JH Yang, HB Brandao, and AS Hansen, DNA double-strand break end synapsis by DNA loop extrusion. Nature Communications 14, 1913 (2023).
# Further adaptation was done from Jusuf et al., NSMB (2025) to include targeted loading of cohesin

# In this condition:
# [Yes] Loop Extrusion
# [Yes] Enhancer-promoter attraction (3 kBT attractive energy)
# [Yes] targeted loading at specified coordinates

from __future__ import absolute_import, division, print_function

import numpy as np
import pyximport; pyximport.install()
pyximport.install(setup_args={'include_dirs': np.get_include()})

from DSB_smcTranslocator_v2 import smcTranslocatorDirectional

from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
import simtk.unit
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file

import copy
import time
import logging
logging.basicConfig(level=logging.INFO)
import sys
import os
from itertools import product

import matplotlib.pyplot as plt

from datetime import datetime
class DSB_Simulation(Simulation):
    def do_block(
            self,
            steps=None,
            check_functions=[],
            get_velocities=False,
            save=True,
            save_extras={},
            positions_to_sample=None,
    ):
        """performs one block of simulations, doing steps timesteps,
        or steps_per_block if not specified.

        Parameters
        ----------

        steps : int or None
            Number of timesteps to perform.
        increment : bool, optional
            If true, will not increment self.block and self.steps counters
        """

        if not self.forces_applied:
            if self.verbose:
                logging.info("applying forces")
                sys.stdout.flush()
            self._apply_forces()
            self.forces_applied = True

        a = time.time()
        self.integrator.step(steps)  # integrate!

        self.state = self.context.getState(
            getPositions=True, getVelocities=get_velocities, getEnergy=True
        )

        b = time.time()
        coords = self.state.getPositions(asNumpy=True)
        newcoords = coords / simtk.unit.nanometer
        newcoords = np.array(newcoords, dtype=np.float32)
        if self.kwargs["save_decimals"] is not False:
            newcoords = np.round(newcoords, self.kwargs["save_decimals"])

        self.time = self.state.getTime() / simtk.unit.picosecond

        # calculate energies in KT/particle
        eK = self.state.getKineticEnergy() / self.N / self.kT
        eP = self.state.getPotentialEnergy() / self.N / self.kT
        curtime = self.state.getTime() / simtk.unit.picosecond

        msg = "block %4s " % int(self.block)
        msg += "pos[1]=[%.1lf %.1lf %.1lf] " % tuple(newcoords[0])

        check_fail = False
        for check_function in check_functions:
            if not check_function(newcoords):
                check_fail = True

        if np.isnan(newcoords).any():
            raise IntegrationFailError("Coordinates are NANs")
        if eK > self.eK_critical:
            print(eK)
            print(self.eK_critical)
            raise EKExceedsError("Ek={1} exceeds {0}".format(self.eK_critical, eK))
        if (np.isnan(eK)) or (np.isnan(eP)):
            raise IntegrationFailError("Energy is NAN)")
        if check_fail:
            raise IntegrationFailError("Custom checks failed")

        dif = np.sqrt(np.mean(np.sum((newcoords - self.get_data()) ** 2, axis=1)))
        msg += "dr=%.2lf " % (dif,)
        self.data = coords
        msg += "t=%2.1lfps " % (self.state.getTime() / simtk.unit.picosecond)
        msg += "kin=%.2lf pot=%.2lf " % (eK, eP)
        msg += "Rg=%.3lf " % self.RG()
        msg += "SPS=%.0lf " % (steps / (float(b - a)))

        if (
                self.integrator_type.lower() == "variablelangevin"
                or self.integrator_type.lower() == "variableverlet"
        ):
            dt = self.integrator.getStepSize()
            msg += "dt=%.1lffs " % (dt / simtk.unit.femtosecond)
            mass = self.system.getParticleMass(1)
            dx = simtk.unit.sqrt(2.0 * eK * self.kT / mass) * dt
            msg += "dx=%.2lfpm " % (dx / simtk.unit.nanometer * 1000.0)

        logging.info(msg)

        if positions_to_sample is None:
            result = {
                "pos": newcoords,
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }
        else:

            result = {
                "pos": [newcoords[pos] for pos in positions_to_sample],
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }

        if get_velocities:
            result["vel"] = self.state.getVelocities() / (
                    simtk.unit.nanometer / simtk.unit.picosecond
            )
        result.update(save_extras)
        if save:
            for reporter in self.reporters:
                reporter.report("data", result)

        self.block += 1
        self.step += steps

        return result



def run_simulation(
        smc_params, stallLeftArray, stallRightArray, boundary_coords, loading_loc_biases,
        PAUSEPROB, L, N_monomers, translocator_initialization_steps,
        smcStepsPerBlock, smcBondDist, smcBondWiggleDist,
        steps_per_block, volume_density, block,
        chain, extra_bonds, save_folder, save_every_x_blocks,
        total_saved_blocks, restartBondUpdaterEveryBlocks, GPU_choice=0,
        overwrite=False, density=0.3, positions_to_sample=None, monomer_types=None,
        interaction_matrix=None, colrate=0.3, errtol=0.01, trunc=3,
        initial_conformation=None, block_to_save_all=None, save_length=100
):
    print(f"Starting simulation. Doing {total_saved_blocks} steps, saving every {save_every_x_blocks} block(s).")
    print(f"A total of {total_saved_blocks} steps will be performed")

    # Assertions to ensure math aligns for bond updater
    assert restartBondUpdaterEveryBlocks % save_every_x_blocks == 0
    assert (total_saved_blocks * save_every_x_blocks) % restartBondUpdaterEveryBlocks == 0

    BondUpdaterInitsTotal = (total_saved_blocks) * save_every_x_blocks // restartBondUpdaterEveryBlocks
    print("BondUpdater will be initialized {0} times".format(BondUpdaterInitsTotal))

    SMCTran = initModel(
        smc_params=smc_params,
        stallLeftArray=stallLeftArray,
        stallRightArray=stallRightArray,
        boundary_coords=boundary_coords,
        loading_loc_biases=loading_loc_biases,
        L=L,
        PAUSEPROB=PAUSEPROB
    )

    print("Equilibrating...")
    SMCTran.steps(translocator_initialization_steps)  # Equilibrate SMC dynamics

    BondUpdater = simulationBondUpdater(SMCTran)

    if initial_conformation is None:
        box = (N_monomers / volume_density) ** 0.33
        data = grow_cubic(N_monomers, int(box) - 2, method='linear')
    else:
        data = initial_conformation

    reporter = HDF5Reporter(folder=save_folder, max_data_length=save_length, overwrite=overwrite, blocks_only=True)

    print("Beginning Saved Steps...")
    for BondUpdaterCount in range(BondUpdaterInitsTotal):

        # Initialize OpenMM Simulation
        sim = DSB_Simulation(
            platform="cuda",
            integrator="variableLangevin",
            error_tol=errtol,
            GPU="{}".format(GPU_choice),
            collision_rate=colrate,
            N=len(data),
            reporters=[reporter],
            PBCbox=False,
            precision="mixed"
        )

        sim.set_data(data)  # Load polymer, center of mass at zero
        sim.add_force(forces.spherical_confinement(sim, density=density, k=1))

        # Add forces, passing the dynamically generated 3D arrays
        sim.add_force(
            forcekits.polymer_chains(
                sim,
                chains=chain,
                bond_force_func=forces.harmonic_bonds,
                bond_force_kwargs={
                    'bondLength': 1.0,
                    'bondWiggleDistance': 0.1,
                    'override_checks': True,
                },
                angle_force_func=forces.angle_force,
                angle_force_kwargs={
                    'k': 0.05,
                    'override_checks': True,
                },
                nonbonded_force_func=forces.heteropolymer_SSW,
                nonbonded_force_kwargs={
                    "interactionMatrix": interaction_matrix,  # Tiled interaction matrix
                    "monomerTypes": monomer_types,  # Tiled sticky elements
                    "extraHardParticlesIdxs": [],
                    "attractionEnergy": 0,
                    "attractionRadius": 2,
                    "repulsionRadius": 1.05,
                    "repulsionEnergy": trunc,
                },
                except_bonds=True,
                override_checks=True,
            )
        )

        sim.step = block

        # Configure SMC Bond parameters
        kbond = sim.kbondScalingFactor / (smcBondWiggleDist ** 2)
        bondDist = smcBondDist * sim.length_scale
        activeParams = {"length": bondDist, "k": kbond}
        inactiveParams = {"length": bondDist, "k": 0}
        BondUpdater.setParams(activeParams, inactiveParams)

        # Apply current bonds from 1D translocator to the 3D simulation
        BondUpdater.setup(
            BondUpdaterCount,
            bondForce=sim.force_dict['harmonic_bonds'],
            smcStepsPerBlock=smcStepsPerBlock,
            blocks=restartBondUpdaterEveryBlocks
        )
        print("Restarting BondUpdater")

        # Minimize energy after restructuring bonds
        sim.local_energy_minimization()

        for b in range(restartBondUpdaterEveryBlocks):
            # BondUpdater updates bonds at each time step.
            curBonds, pastBonds = BondUpdater.step(sim.context)
            if (b % save_every_x_blocks == 0):
                # save SMC positions and monomer positions
                total_blocks_so_far = BondUpdaterCount * restartBondUpdaterEveryBlocks + b
                if total_blocks_so_far in block_to_save_all:
                    sim.do_block(steps=steps_per_block, save_extras={"SMCs": curBonds, "SMC_step": sim.step},
                                 positions_to_sample=[i for i in range(N_monomers)])
                else:
                    sim.do_block(steps=steps_per_block, save_extras={"SMCs": curBonds, "SMC_step": sim.step},
                                 positions_to_sample=positions_to_sample)
            else:
                sim.integrator.step(steps_per_block)  # do steps without getting the positions from the GPU (faster)

                # Extract data to pass to the next BondUpdater chunk
        data = sim.get_data()
        block = sim.step
        del sim
        time.sleep(0.2)

        # Dump final data and close out
    reporter.dump_data()
    with open(os.path.join(save_folder, 'sim_done.txt'), "w+") as done_file:
        pass

def initModel(smc_params, stallLeftArray, stallRightArray, boundary_coords, loading_loc_biases, L, PAUSEPROB):
    # unchanging
    BELT_ON = 0
    BELT_OFF = 1
    switchRate = 0
    SWITCH_PROB = switchRate
    PUSH = 0
    PAIRED = 0
    SLIDE = 1
    SLIDE_PAUSEPROB = 0.99
    loop_prefactor = 1.5
    FULL_LOOP_ENTROPY = 1
    FRACTION_ONESIDED = 0

    # sweep
    processivity = smc_params['processivity']
    separations = smc_params['separations']
    longlived_fraction = smc_params.get('longlived_fraction', 0)
    longlived_boost = smc_params.get('longlived_boost_factor', 1)
    lifetime_extension = smc_params.get('ctcf_boost_factor', 1)

    # birth array (loading targeted)
    if loading_loc_biases is not None:
        # for loader_list in smc_params['loading_loc_biases']:
        #     loader_loc, len_bias, bias = loader_list[0], loader_list[1], loader_list[2]
        #
        #     start = max(0, loader_loc - len_bias)
        #     end = min(L, loader_loc + len_bias)
        #     birthArray[start:end] = bias
        birthArray = loading_loc_biases.copy()
        birthArray = birthArray / sum(birthArray)  # Normalize
    else:
        birthArray = np.ones(L)/L

    # Base rate (times 0.5 to account for two-sided extrusion)
    base_death_rate = 1. / (0.5 * processivity / (1 - PAUSEPROB))
    deathArray = np.zeros(L, dtype=np.double) + base_death_rate

    if len(boundary_coords) > 0:
        deathArray[boundary_coords] = base_death_rate / lifetime_extension

    if longlived_fraction == 0:
        smcNum = int(L // separations)
        longlived_smcNum = 0
        longlived_deathArray = np.zeros(L, dtype=np.double) + base_death_rate
    else:
        normal_sep = separations / (1 - longlived_fraction)
        smcNum = int(L // normal_sep)

        longlived_sep = separations / longlived_fraction
        longlived_smcNum = int(L // longlived_sep)

        longlived_base_rate = base_death_rate / longlived_boost
        longlived_deathArray = np.zeros(L, dtype=np.double) + longlived_base_rate

        if len(boundary_coords) > 0:
            longlived_deathArray[boundary_coords] = longlived_base_rate / lifetime_extension

    # --- 5. Extruder Dynamics Frequency Arrays ---
    SWITCH = np.ones(L, dtype=np.double) * SWITCH_PROB
    pauseArray = PAUSEPROB * np.ones(L, dtype=np.double)
    slidePauseArray = np.zeros(L, dtype=np.double) + SLIDE_PAUSEPROB
    oneSidedArray = np.zeros(smcNum, dtype=np.int64)
    longlived_oneSidedArray = np.zeros(longlived_smcNum, dtype=np.int64)
    belt_on_array = np.zeros(smcNum, dtype=np.double) + BELT_ON
    belt_off_array = np.zeros(smcNum, dtype=np.double) + BELT_OFF

    # Calculate directional sliding pause probabilities
    spf = slidePauseArray * (1. - (1. - SLIDE_PAUSEPROB) * np.exp(-1. * loop_prefactor))
    spb = slidePauseArray * (1. - (1. - SLIDE_PAUSEPROB) * np.exp(loop_prefactor))

    # --- 6. Initialize Translocator ---
    transloc = smcTranslocatorDirectional(
        birthArray,
        deathArray,
        longlived_deathArray,
        stallLeftArray,  # The fully tiled arrays
        stallRightArray,  # The fully tiled arrays
        pauseArray,
        smcNum,
        longlived_smcNum,
        oneSidedArray,
        longlived_oneSidedArray,
        FRACTION_ONESIDED,
        slide=SLIDE,
        slidepauseForward=spf,
        slidepauseBackward=spb,
        switch=SWITCH,
        pushing=PUSH,
        belt_on=belt_on_array,
        belt_off=belt_off_array,
        SLIDE_PAUSEPROB=SLIDE_PAUSEPROB
    )

    return transloc

class simulationBondUpdater(object):
    """
    This class precomputes simulation bonds for faster dynamic allocation.
    """

    def __init__(self, smcTransObject):  # , plectonemeObject):
        """
        :param smcTransObject: smc translocator object to work with
        """
        self.smcObject = smcTransObject
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict

    def setup(self, BondUpdaterCount, bondForce, smcStepsPerBlock, blocks=100):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """

        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce

        # precalculating all bonds
        allBonds = []

        left, right = self.smcObject.getSMCs()
        # add SMC bonds
        bonds = [(int(i), int(j)) for i, j in zip(left, right)]

        left, right = self.smcObject.getlonglivedSMCs()
        # add longlivedSMC bonds
        bonds += [(int(i), int(j)) for i, j in zip(left, right)]
        self.curBonds = bonds.copy()

        allBonds.append(bonds)
        for dummy in range(blocks):
            self.smcObject.steps(smcStepsPerBlock)
            left, right = self.smcObject.getSMCs()
            # add SMC bonds
            bonds = [(int(i), int(j)) for i, j in zip(left, right)]

            left, right = self.smcObject.getlonglivedSMCs()
            # add longlivedSMC bonds
            bonds += [(int(i), int(j)) for i, j in zip(left, right)]

            allBonds.append(bonds)

        self.allBonds = allBonds
        self.uniqueBonds = list(set(sum(allBonds, [])))

        allBonds.pop(0)

        # adding forces and getting bond indices
        self.bondInds = []

        for bond in self.uniqueBonds:
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict
            ind = bondForce.addBond(bond[0], bond[1], **paramset)  # changed from addBond
            self.bondInds.append(ind)
        self.bondToInd = {i: j for i, j in zip(self.uniqueBonds, self.bondInds)}
        return self.curBonds, []

    def step(self, context, verbose=False):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds]
        bondsAdd = [i for i in self.curBonds if i not in pastBonds]
        bondsStay = [i for i in pastBonds if i in self.curBonds]
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove)
        for bond, isAdd in zip(bondsToChange, bondsIsAdd):
            ind = self.bondToInd[bond]
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds


def generate_param_grid(param_dict):
    """Yields a dictionary for each combination in the sweep space."""
    keys, values = zip(*param_dict.items())
    for bundle in product(*values):
        yield dict(zip(keys, bundle))


def build_tiled_boundaries(ctcf_configs, chrom_size, region_size):
    """Tiles multiple CTCF configurations across the chromosome."""
    num_regions = chrom_size // region_size
    boundaryStrengthsL = np.zeros(chrom_size, dtype=np.double)
    boundaryStrengthsR = np.zeros(chrom_size, dtype=np.double)
    boundary_coordinates = []

    for i in range(num_regions):
        region_start = i * region_size
        config = ctcf_configs[i % len(ctcf_configs)]

        abs_sites_L = region_start + config['sites_L']
        abs_sites_R = region_start + config['sites_R']

        boundaryStrengthsL[abs_sites_L] = config['probs_L']
        boundaryStrengthsR[abs_sites_R] = config['probs_R']

        boundary_coordinates.extend(abs_sites_L)
        boundary_coordinates.extend(abs_sites_R)

    return boundaryStrengthsL, boundaryStrengthsR, np.array(boundary_coordinates, dtype=int)


def build_tiled_monomer_types(sticky_configs, chrom_size, region_size):
    """Tiles sticky elements (e.g., enhancers/promoters) across the chromosome."""
    num_regions = chrom_size // region_size
    monomer_types = np.zeros(chrom_size, dtype=int)

    for i in range(num_regions):
        region_start = i * region_size
        config = sticky_configs[i % len(sticky_configs)]
        abs_sites = region_start + config['sites']
        monomer_types[abs_sites] = config['types']

    return monomer_types

def build_tiled_mm(mm_configs, chrom_size, region_size, loading_bias, loading_width):
    """Tiles the targeted loading elements across the chromosome."""
    num_regions = chrom_size // region_size
    mmLoc = np.ones(chrom_size, dtype=np.double)

    for i in range(num_regions):
        region_start = i * region_size
        config = mm_configs[i % len(mm_configs)]
        abs_sites = region_start + config['sites']

        for site in abs_sites:
            mmLoc[(site-loading_width//2):(site+loading_width//2)] = loading_bias

    return mmLoc # this fulfills the need for birthArray

def build_tiled_mm_permm(mm_configs, chrom_size, region_size, cohesins_per_mm, loading_width):
    """Tiles the targeted loading elements across the chromosome."""
    num_regions = chrom_size // region_size
    mmLoc = np.ones(chrom_size, dtype=np.double)

    loading_bias = cohesins_per_mm / loading_width

    for i in range(num_regions):
        region_start = i * region_size
        config = mm_configs[i % len(mm_configs)]
        abs_sites = region_start + config['sites']

        for site in abs_sites:
            mmLoc[(site-loading_width//2):(site+loading_width//2)] = loading_bias

    return mmLoc, loading_bias # this fulfills the need for birthArray, but is not normalized

def test_sweep(test_idx, region_size, CTCF_L, CTCF_R, sticky_elts, birthArray):
    """
    Plots a 1D 'dummy map' of the simulation parameters for a specific region.
    """
    test_region_start = region_size * test_idx
    test_region_end = test_region_start + region_size

    # 1. Slice out the specific region from the full chromosome arrays
    x = np.arange(test_region_start, test_region_end)
    ctcf_l_slice = CTCF_L[test_region_start:test_region_end]
    ctcf_r_slice = CTCF_R[test_region_start:test_region_end]
    sticky_slice = sticky_elts[test_region_start:test_region_end]
    birth_slice = birthArray[test_region_start:test_region_end]
    birth_slice = birth_slice / sum(birth_slice)  # Normalize
    print(sum(birth_slice))

    # 2. Set up a stacked plot sharing the same X-axis
    fig, axes = plt.subplots(4, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f'Region {test_idx} (Monomers {test_region_start} - {test_region_end})', fontsize=14)

    # Track 1: Left-pointing CTCF (Stall Right for translocator)
    axes[0].stem(x, ctcf_l_slice, linefmt='b-', markerfmt='bo', basefmt='k-', use_line_collection=True)
    axes[0].set_ylabel('CTCF Left\nProb')
    axes[0].set_ylim(0, 1.1)
    axes[0].grid(axis='x', linestyle='--', alpha=0.5)

    # Track 2: Right-pointing CTCF (Stall Left for translocator)
    axes[1].stem(x, ctcf_r_slice, linefmt='r-', markerfmt='ro', basefmt='k-', use_line_collection=True)
    axes[1].set_ylabel('CTCF Right\nProb')
    axes[1].set_ylim(0, 1.1)
    axes[1].grid(axis='x', linestyle='--', alpha=0.5)

    # Track 3: Sticky Elements (Enhancers/Promoters)
    # Using fill_between creates a nice blocky look for continuous regions
    axes[2].fill_between(x, 0, sticky_slice, color='green', step="mid", alpha=0.7)
    axes[2].set_ylabel('Sticky\nElements')
    axes[2].set_yticks([0, 1])
    axes[2].set_ylim(0, 1.2)
    axes[2].grid(axis='x', linestyle='--', alpha=0.5)

    # Track 4: Targeted Loading Bias (birthArray)
    axes[3].plot(x, birth_slice, color='purple', drawstyle='steps-mid')
    axes[3].fill_between(x, 0, birth_slice, color='purple', alpha=0.3, step="mid")
    axes[3].set_ylabel('Loading Bias\n(birthArray)')

    # We let matplotlib auto-scale the Y-axis here because if birthArray
    # was already normalized (birthArray / sum), the values will be tiny.
    axes[3].set_ylim(bottom=0)
    axes[3].grid(axis='x', linestyle='--', alpha=0.5)

    # Clean up the bottom axis
    axes[3].set_xlabel('Genomic Coordinate (Monomer Index)', fontsize=12)

    plt.tight_layout()
    plt.show()

save_base_folder = f"/mnt/md1/varshini/Blood/sim_data_sweep_cohesins_single/{datetime.today()}"
chrom_size = 70000
region_size = 2000

num_chains = 1
chain = [(x * chrom_size, (x + 1) * chrom_size, 0) for x in range(num_chains)]

N_monomers = chain[-1][1]
pause_prob = 1 - 0.0025

translocator_initialization_steps = 10000
smcStepsPerBlock = 1
steps_per_block = 50 # about 2.5 ms per block
restartBondUpdaterEveryBlocks = 3000
save_every_x_blocks = 3000 # save every 7.5 seconds
total_saved_blocks = 1200
smcBondDist = 0.5
smcBondWiggleDist = 0.2
volume_density = 0.3
GPU_choice = 0

collision_rate = 0.03
truncated_potentials = 3

# 1D Boundary (CTCF) Configurations
ctcf_configs = [{
    'name': 'default',
    'sites_L': np.array([574, 694, 866, 1241, 1390, 1580, 1752, 1800]),
    'probs_L': np.array([0.6, 0.8 , 0.95, 0.1 , 0.6, 0.6, 0.8, 0.1 ]),
    'sites_R': np.array([200, 330, 724, 1425, 1433, 1604]),
    'probs_R': np.array([0.9, 0.3 , 0.95, 0.4 , 0.3, 0.4 ])
}]

# 3D Sticky Element Configurations
sticky_configs = [{
    'name': 'default',
    'sites': np.array([250, 372, 540, 745, 775, 833, 961, 1202, 1330, 1640, 1722]),
    'types': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
}]

# matchmaker locations
mm_configs = [{
    'name': 'default',
    'sites': np.array([456, 1054, 1507]),
}]


# Interaction Matrix (Type 0 = background, Type 1 = sticky)
EP_interaction_energy = 3.0
interaction_matrix = np.array([
    [0.0, 0.0],
    [0.0, EP_interaction_energy]
])


simulation_sweep_space = {
    'processivity': [300],
    'separations': [240],
    'ctcf_boost_factor': [4],
    'longlived_fraction': [0],
    'longlived_boost_factor': [20],
    'dsb_boost_factor': [0],
    'llp': [2],
    'h': [0.4],
    #'loading_b': [2, 4, 8, 16],  # targeted loading bias
    #'loading_width': [4, 8, 16, 32] # width of targeted loading in kb,
    'loading_width': [16], # width of targeted loading in kb,
    'cohesins_per_mm':[64] # if none, that means no bias
}


if __name__ == '__main__':

    # Pre-build the tiled arrays for 1D and 3D
    stall_L, stall_R, boundary_coords = build_tiled_boundaries(ctcf_configs, chrom_size, region_size)
    monomer_types = build_tiled_monomer_types(sticky_configs, chrom_size, region_size)

    positions_to_sample = np.arange(chrom_size)
    block_to_save_all = []

    # Iterate through the cleanly generated parameter grid
    for config in generate_param_grid(simulation_sweep_space):

        # Build the dynamic directory name based on current sweep parameters
        run_name = (f"sim_proc{config['processivity']}_sep{config['separations']}"
                    f"_ctcf{config['ctcf_boost_factor']}_n{config['cohesins_per_mm']}"
                    f"_h{config['h']}_llp{config['llp']}_w{config['loading_width']}")
        save_folder = os.path.join(save_base_folder, run_name)

        if not os.path.exists(save_folder):
            os.makedirs(save_folder)

        print(f"\n{'=' * 50}")
        print(f"Starting Run: {run_name}")
        print(f"Config: {config}")
        print(f"{'=' * 50}")

        # Construct loading biases (using the swept loading_b)
        # if config['loading_b'] != 1:
        #     loading_loc_biases = build_tiled_mm(mm_configs, chrom_size, region_size, config['loading_b'], config['loading_width'])
        # else:
        #     loading_loc_biases = None

        if config['cohesins_per_mm'] is not None:
            loading_loc_biases, loading_bias = build_tiled_mm_permm(mm_configs, chrom_size, region_size,
                                                      config['cohesins_per_mm'], config['loading_width'])
        else:
            loading_loc_biases = None
            loading_bias = 1
        # Pack SMC parameters for this run
        smc_params = {
            'processivity': config['processivity'],
            'separations': config['separations'],
            'ctcf_boost_factor': config['ctcf_boost_factor'],
            'longlived_fraction': config['longlived_fraction'],
            'longlived_boost_factor': config['longlived_boost_factor'],
            'dsb_boost_factor': config['dsb_boost_factor'],
            'loading_loc_biases': loading_loc_biases
        }

        # Save a text log of the exact parameters used for this specific run
        with open(os.path.join(save_folder, "run_parameters.txt"), "w") as f:
            f.write(str(config))
            f.write("\n Loading Bias:")
            f.write(str(loading_bias))

        # Launch the simulation
        #test_sweep(2, region_size, stall_L, stall_R, monomer_types, loading_loc_biases)
        run_simulation(
            smc_params=smc_params,
            stallLeftArray=stall_R,  # Map standard L/R to translocator logic
            stallRightArray=stall_L,
            boundary_coords=boundary_coords,
            loading_loc_biases=loading_loc_biases,
            PAUSEPROB=pause_prob,
            L=chrom_size,
            N_monomers=N_monomers,
            monomer_types=monomer_types,
            interaction_matrix=interaction_matrix,
            translocator_initialization_steps=translocator_initialization_steps,
            smcStepsPerBlock=smcStepsPerBlock,
            smcBondDist=smcBondDist,
            smcBondWiggleDist=smcBondWiggleDist,
            steps_per_block=steps_per_block,
            volume_density=volume_density,
            block=0,  # Initial step
            chain=chain,
            extra_bonds=[],
            save_folder=save_folder,
            save_every_x_blocks=save_every_x_blocks,
            total_saved_blocks=total_saved_blocks,
            restartBondUpdaterEveryBlocks=restartBondUpdaterEveryBlocks,
            GPU_choice=GPU_choice,
            overwrite=True,
            positions_to_sample=positions_to_sample,
            block_to_save_all=block_to_save_all,
            save_length=10,
            colrate=collision_rate,
            trunc=truncated_potentials,
        )