from covid19_sim.simulation import *

def sample_experiment():
    def sim1():
        return Simulation(
            label = 'No intervention',
            initial_state = [0.97, 0.03, 0],
            width = 30, height = 30,
            normalize_infected = True,
        )

    def sim2():
        return Simulation(
            label = '30% general quarantine',
            initial_state = [0.97, 0.03, 0],
            width = 30, height = 30,
            normalize_infected = True,
            general_quarantine = True,
            general_quarantine_percent = 0.3,
            quarantine_duration = 9999
        )
    return [sim1, sim2]