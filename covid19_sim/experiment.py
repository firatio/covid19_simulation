import matplotlib.pyplot as plt
from matplotlib import colors
import math
import numpy as np
from textwrap import wrap
import statistics
import os
import time
import json
import datetime
from pathlib import Path

desc = []   # description for each experiment
ti = []     # average total infected for each experiment
tfq = []    # average total unnecessarily quarantined for each experiment
distributions = []  # several stats for each experiment
sim_properties = [] # properties of simulations in each experiment

def create_bar_chart(results):
    """Create bar chart for the experiments."""
    labels = results[0]
    width = math.ceil( 120 / len(labels) )
    if width < 20:
        width = 20
    labels = [ '\n'.join(wrap(l, width)) for l in labels ]

    ti_means = results[1]
    tfq_means = results[2]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, ti_means, width, label='Total infected')
    rects2 = ax.bar(x + width/2, tfq_means, width, label='Total unnecessarily quarantined')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('% of population')
    ax.set_title('Comparison of different strategies')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend()

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    
    fig.tight_layout()
    plt.ylim((0, 100))
    plt.show()

def run_sim(sim_function, n=100):
    """Run the given simulation n times."""
    last_time = time.time()
    # list of results for each simulation
    stats = []

    total_infected = total_unnecessarily_quarantined = 0
    total_unnecessary_quarantine_days = total_tests = epidemic_duration = 0
    total_infected_list = []
    for i in range(1, n+1):
        if i % 25 == 0:
            #os.system('cls')
            print(str(i) + '/' + str(n), flush='True')
            now = time.time()
            print(str(now - last_time) + ' s', flush='True')
            last_time = now
        sim = sim_function()
        sim.run()

        population_size = len(sim.population)
        total_infected += sim.total_infected
        total_infected_list.append(sim.total_infected/population_size)
        total_unnecessarily_quarantined += sim.total_unnecessarily_quarantined
        total_unnecessary_quarantine_days += sim.total_unnecessary_quarantine_days
        total_tests += sim.total_tests
        epidemic_duration += sim.epidemic_duration
        stats.append([
            sim.total_infected/population_size,
            sim.total_unnecessarily_quarantined/population_size,
            sim.total_unnecessary_quarantine_days,
            sim.total_tests,
            sim.epidemic_duration,
            sim.max_infected])
    
    avg_total_infected = round(100*total_infected / (n*population_size), 2)
    avg_total_unnecessarily_quarantined = round(100*total_unnecessarily_quarantined / (n*population_size), 2)
    avg_total_unnecessary_quarantine_days = math.ceil(total_unnecessary_quarantine_days / n)
    avg_total_tests = math.ceil(total_tests / n)
    avg_epidemic_duration = math.ceil(epidemic_duration / n)

    ti.append(avg_total_infected)
    tfq.append(avg_total_unnecessarily_quarantined)

    description = (
        sim.get_description()
        + ' Total test sessions: ' + str(avg_total_tests)
        + '. Total unnecessary quarantine days: ' + str(avg_total_unnecessary_quarantine_days)
        + '. Epidemic duration: ' + str(avg_epidemic_duration)
        + ' days. Total infected median: ' + str(round(statistics.median(total_infected_list)*100, 2)) + '%'
    )
    desc.append(description)

    # create sim description dictionary using the last simulation object
    sim_description_dict = {
        "label": sim.label,
        "width": sim.width,
        "height": sim.height,
        "initial_state": sim.initial_state,
        "normalize_infected": sim.normalize_infected,
        "infection_risk": sim.infection_risk,
        "min_disease_duration": sim.min_disease_duration,
        "max_disease_duration": sim.max_disease_duration,
        "asymptomatic_ratio": sim.asymptomatic_ratio,
        "diagnostic_test_info": sim.diagnostic_test.get_properties(),
        "long_distance_travel": sim.long_distance_travel,
        "air_travel": sim.air_travel,
        "intervention": sim.intervention,
        "general_quarantine": sim.general_quarantine,
        "general_quarantine_percent": sim.general_quarantine_percent,
        "quarantine_duration": sim.quarantine_duration,
        "quarantine_duration_after_negative_test": sim.quarantine_duration_after_negative_test,
        "quarantine_duration_after_positive_test": sim.quarantine_duration_after_positive_test,
        "testing_strategy": sim.testing_strategy,
        "screening_percentage": sim.screening_percentage,
        "max_test_percentage": sim.max_test_percentage,
        "test_repeat_period": sim.test_repeat_period,
        "test_before_travel": sim.test_before_travel,
        "isolate_travelers": sim.isolate_travelers,
        "probabilities_to_go_to_hospital": sim.probabilities_to_go_to_hospital,
        "probability_to_be_traced": sim.probability_to_be_traced,
        "probability_to_be_reached": sim.probability_to_be_reached,
        "similar_symptoms_probability": sim.similar_symptoms_probability,
        "probability_to_put_contacts_in_quarantine": sim.probability_to_put_contacts_in_quarantine,
    }
    sim_properties.append(sim_description_dict)
    
    distributions.append(stats)

def run_experiment(experiment, n=200):
    """Run the simulations n times and save results to a json file.
    
    experiment: a function that returns a list of simulations
    n: number of times each simulation will be run
    """
    sims = experiment()
    size = len(sims)
    i = 1
    for sim in sims:
        print('Simulation ' + str(i) + '/' + str(size), flush = 'True')
        run_sim(sim, n)
        i += 1
    
    # save experiment results as a json file
    results_folder = Path("results/")
    results = [desc, ti, tfq, distributions, sim_properties]
    dtnow = datetime.datetime.now()
    filename = experiment.__name__
    stamped_filename = results_folder / (dtnow.strftime('%Y%m%d_%H%M%S_') + filename + '.json')
    with open(stamped_filename, 'w') as outfile:
        json.dump(results, outfile)