from covid19_sim.simulation import Simulation, Test

sim = Simulation(
    initial_state = [0.97, 0.03, 0],
    normalize_infected = True,
    width = 100, height = 100,
    air_travel = [False, 0, 0],
    general_quarantine = True,
    general_quarantine_percent = 0.3,
    quarantine_duration = 9999,
    intervention = [True, 14, 21],
)

"""
sim = Simulation(
    initial_state = [0.97, 0.03, 0],
    normalize_infected = True,
    width = 100, height = 100,
    air_travel = [True, 0.01, 0.3],
    test_before_travel = True,
    isolate_travelers = True,
    general_quarantine = True,
    general_quarantine_percent = 0.3,
    quarantine_duration = 9999,
    diagnostic_test = Test(
        accuracy = 0.98, 
        false_positive_probability = 0.01,
        use_double_test = False, 
        require_double_positive = False
    ),
    intervention = [True, 14, 21],
    quarantine_duration_after_positive_test = 14,
)
"""

sim.run(generate_graph = True)