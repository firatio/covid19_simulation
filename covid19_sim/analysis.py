from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import math
from textwrap import wrap
import json
import glob
import os

def monte_carlo_permutation_test(group_1, group_2, n_random):
    """
    An implementation of 2-tailed monte carlo permutation test.
    Compares the absolute value of the difference in means of 2 groups
    and returns a p-value.

    group_1: list for group 1
    group_2: list for group 2
    n_random: number of random samples generated from both groups
    """
    group_1_size, counter = len(group_1), 0
    observed_difference = np.abs(np.mean(group_1) - np.mean(group_2))
    pooled_group = np.concatenate([group_1, group_2])

    for i in range(n_random):
        np.random.shuffle(pooled_group)
        if np.abs(np.mean(pooled_group[:group_1_size]) - np.mean(pooled_group[group_1_size:])) >= observed_difference:
            counter += 1

    return counter / n_random

class Metric(Enum):
    TOTAL_INFECTED = 1
    MAX_INFECTED = 2
    TOTAL_UNNECESSARILY_QUARANTINED = 3
    DAYS_LOST_TO_UNNECESSARY_QUARANTINE = 4
    EPIDEMIC_DURATION = 5
    TOTAL_TEST_SESSIONS = 6

def create_custom_bar_chart(results, **kwargs):
    """Create bar chart for the experiments."""
    metrics = kwargs.get('metrics', [Metric.TOTAL_INFECTED, Metric.TOTAL_UNNECESSARILY_QUARANTINED])
    custom_labels = kwargs.get('custom_labels', [])

    #labels = results['labels']
    labels = ['Simulation ' + str(i) for i in range(1, len(results['labels']) + 1)]
    number_of_experiments = len(labels)
    if custom_labels and len(custom_labels) == number_of_experiments:
        labels = custom_labels
        
    width = math.ceil( 120 / len(labels) )
    if width < 20:
        width = 20
    labels = [ '\n'.join(wrap(l, width)) for l in labels ]

    data = []
    for metric in metrics:
        if metric == Metric.TOTAL_INFECTED:
            label = 'Total infected'
            values = results[Metric.TOTAL_INFECTED]
        if metric == Metric.MAX_INFECTED:
            label = 'Max infected'
            values = results[Metric.MAX_INFECTED]
        if metric == Metric.TOTAL_UNNECESSARILY_QUARANTINED:
            label = 'Total unnecessarily quarantined'
            values = results[Metric.TOTAL_UNNECESSARILY_QUARANTINED]
        if metric == Metric.EPIDEMIC_DURATION:
            label = 'Epidemic duration'
            values = results[Metric.EPIDEMIC_DURATION]
        data.append({
            'label': label,
            'values': values
        })
    
    x = np.arange(number_of_experiments)  # label locations
    width = 0.2 # the width of the bars

    fig, ax = plt.subplots()
    fig.set_size_inches(15, 5.5)
    metric_bars = []
    number_of_metrics = len(data)
    for i in range(number_of_metrics):
        metric_bars.append(ax.bar(x + (2*i - number_of_metrics + 1)*(width/2), data[i]['values'], width, label=data[i]['label']))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('% of population', fontsize=18)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_title('Comparison of different strategies', fontsize = 18)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=18)
    ax.legend(prop={'size': 18})

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        size = 18,
                        ha='center', va='bottom')
    for bars in metric_bars:
        autolabel(bars)
    
    fig.tight_layout()
    plt.xlim((-0.4, number_of_experiments - 0.3))
    plt.ylim((0, 109))
    plt.show()

def run_analysis(file_list):
    """Open files, load data and analyze them"""
    if not file_list:
        list_of_files = glob.glob('results/*.json')
        latest_file = max(list_of_files, key=os.path.getctime)
        if not latest_file:
            raise Exception('No files given to analyze!')
        file_list = [latest_file]
    
    results = []
    data = []
    for filename in file_list:
        file = open(filename, 'r')
        file_data = json.loads(file.read()) 
        file.close()
        data.append(file_data)
    
    # convert data in files into one single list
    n_elements = len(data[0])   # number of elements in each file_data
    for i in range(n_elements):
        merged = []
        for file_data in data:
            merged += file_data[i]
        results.append(merged)

    # ANALYSIS
    distributions = results[3]
    number_of_experiments = len(distributions)
    sim_properties = results[4]

    # check if population sizes of simulation data match
    width, height = sim_properties[0]['width'], sim_properties[0]['height']
    for i in range(number_of_experiments):
        if sim_properties[i]['width'] != width or sim_properties[i]['height'] != height:
            print('WARNING! POPULATION SIZES OF SIMULATIONS DO NOT MATCH!')

    # todo... hardcoded right now, change this
    total_infected_list = [[] for i in range(number_of_experiments)]
    epidemic_duration_list = [[] for i in range(number_of_experiments)]
    max_infected_list = [[] for i in range(number_of_experiments)]
    for j in range(number_of_experiments):
        for i in range(len(distributions[j])):
            total_infected_list[j].append(distributions[j][i][0])
            epidemic_duration_list[j].append(distributions[j][i][-2])
            max_infected_list[j].append(distributions[j][i][-1])

    avg_total_infected_list, avg_epidemic_duration_list, avg_max_infected_list = [], [], []
    for i in range(number_of_experiments):
        population_size = sim_properties[i]['width'] * sim_properties[i]['height']
        avg_total_infected_list.append(round(100*np.mean(total_infected_list[i]), 2))
        avg_epidemic_duration_list.append(round(np.mean(epidemic_duration_list[i]), 2))
        avg_max_infected_list.append(round(100*np.mean(max_infected_list[i])/population_size, 2))

    # compare results of experiments with each other
    for i in range(number_of_experiments):
        for j in range(i + 1, number_of_experiments):
            print()
            print('experiment ' + str(i + 1) + ' vs ' + str(j + 1))

            print('test for total infected: ' + str(avg_total_infected_list[i]) + ' vs ' + str(avg_total_infected_list[j]))
            p = monte_carlo_permutation_test(total_infected_list[i], total_infected_list[j], 5000)
            print('p = ' + str(p))

            print('test for epidemic duration: ' + str(avg_epidemic_duration_list[i]) + ' vs ' + str(avg_epidemic_duration_list[j]))
            p = monte_carlo_permutation_test(epidemic_duration_list[i], epidemic_duration_list[j], 5000)
            print('p = ' + str(p))

            print('test for max infected: ' + str(avg_max_infected_list[i]) + ' vs ' + str(avg_max_infected_list[j]))
            p = monte_carlo_permutation_test(max_infected_list[i], max_infected_list[j], 5000)
            print('p = ' + str(p))
    
    results_dict = {}
    results_dict['labels'] = results[0]
    results_dict[Metric.TOTAL_INFECTED] = results[1]
    results_dict[Metric.TOTAL_UNNECESSARILY_QUARANTINED] = results[2]
    results_dict[Metric.MAX_INFECTED] = avg_max_infected_list
    results_dict[Metric.EPIDEMIC_DURATION] = avg_epidemic_duration_list

    labels = []
    for sim in sim_properties:
        if not 'label' in sim:
            labels = []
            break
        labels.append(sim['label'])

    create_custom_bar_chart(results_dict, 
        metrics = [Metric.TOTAL_INFECTED, Metric.MAX_INFECTED], 
        custom_labels = labels
        )
    #create_custom_bar_chart(results_dict, metrics = [Metric.TOTAL_INFECTED, Metric.MAX_INFECTED, Metric.TOTAL_UNNECESSARILY_QUARANTINED, Metric.EPIDEMIC_DURATION])
    #create_custom_bar_chart(results_dict)