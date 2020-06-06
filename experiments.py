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

# EXPERIMENTS USED IN ARTICLES
"""
Article: Air travel can be used to fight covid-19

What are the effects of testing travelers
before flights and isolating them if positive?
"""
def air_travel_experiment_1a():
    """
    For version 0.2
    Air travel disabled. 
    """
    def sim1():
        return Simulation(
            label = 'Air travel disabled',
            initial_state = [0.97, 0.03, 0],
            width = 100, height = 100,
            normalize_infected = True,
            air_travel = [False, 0, 0],
        )
    return [sim1]

def air_travel_experiment_1b():
    """
    For version 0.2
    Air travel enabled. Test before travel, isolate if positive.
    """
    def sim1():
        return Simulation(
            label = 'Air travel enabled. Test before travel, isolate if positive.',
            initial_state = [0.97, 0.03, 0],
            width = 100, height = 100,
            normalize_infected = True,
            air_travel = [True, 0.01, 0.3],
            test_before_travel = True,
            isolate_travelers = True,
            diagnostic_test = Test(
                accuracy = 0.98, 
                false_positive_probability = 0.01,
                use_double_test = False,
                require_double_positive = False
            ),
            quarantine_duration_after_positive_test = 14,
        )
    return [sim1]

"""
Article: Testing can bring back freedom of travel

What are the effects of long distance travel and testing travelers
before travel and isolating them if positive?
"""
def travel_experiment_1a():
    """
    For version 0.1
    Long distance travel disabled. 
    """
    def sim():
        return Simulation(
            label = 'Air travel disabled',
            initial_state = [0.97, 0.03, 0],
            normalize_infected = True,
            width = 100, height = 100,
            long_distance_travel = [False, 0, 0],
        )
    
    return [sim]

def travel_experiment_1b():
    """
    For version 0.1
    Long distance travel enabled. Travelers tested before travel AND 
    isolated if positive. High accuracy test applied once.
    """
    def sim():
        return Simulation(
            label = 'Air travel enabled. Travelers tested before \
            travel and isolated if positive',
            initial_state = [0.97, 0.03, 0],
            normalize_infected = True,
            width = 100, height = 100,
            long_distance_travel = [True, 0.01, 0.3],
            test_before_travel = True,
            isolate_travelers = True,
            quarantine_duration_after_positive_test = 14,
            diagnostic_test = Test(
                accuracy = 0.98, 
                false_positive_probability = 0.01,
                use_double_test = False, 
                require_double_positive = False
            ),
        )
    return [sim]

def travel_experiment_1c():
    """
    For version 0.1
    Long distance travel enabled. Travelers tested before travel AND 
    isolated if positive. Lower accuracy test applied twice.
    """
    def sim():
        return Simulation(
            label = 'Air travel enabled. Travelers tested before \
            travel and isolated if positive',
            initial_state = [0.97, 0.03, 0],
            normalize_infected = True,
            width = 100, height = 100,
            long_distance_travel = [True, 0.01, 0.3],
            test_before_travel = True,
            isolate_travelers = True,
            quarantine_duration_after_positive_test = 14,
            diagnostic_test = Test(
                accuracy = 0.80, 
                false_positive_probability = 0.01,
                use_double_test = True, 
                require_double_positive = False
            ),
        )
    return [sim]