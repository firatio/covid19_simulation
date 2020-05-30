# COVID-19 SIMULATION

This is an agent-based covid-19 simulation that can be used to discover which strategies can contain the covid-19 outbreak and minimize damage to our freedoms.

Although the simulation is based on a toy model, the positive effects of certain strategies should be applicable to more realistic models.

Discussions of some experiments run with this simulation can be found at https://www.forastelo.org/covid-19. New experiments and articles will be added in the future.

# Notes
- The model is very simple, ignores many details and is still being improved.
- The code is work in progress, needs lots of refactoring and optimization. Some portions of the code are not thoroughly tested and should not be used. See 'How to use the simulation' section.

# Requirements
The code has been run successfully with the following configuration:
- python >= 3.6.5 (64-bit)
- matplotlib 3.2.1
- numpy 1.18.2

# How to use the simulation
## Running a single simulation
- Navigate to the root of the repository
- Change the parameters of the sample simulation in [run_single.py](run_single.py) file
- Run the script
```bash
python run_single.py
```

Simulation will show how the state of the society changes in time and will create a summary graph with important metrics.

## Running an experiment and analyzing the results
- Modify the sample experiment in the [experiments.py](experiments.py) file
    - The resulting graph is optimized for 2 simulations, so it is better to compare only 2 simulations in a single experiment
    - As the world size and number of simulations get larger, it takes much more time to run an experiment
- Navigate to the root of the repository and run the experiment script
```bash
python run_experiment.py
```

Once the simulations end, a json output file will be generated in the results folder.

- Execute 'run_analysis.py' script to see an analysis of the results
```bash
python run_analysis.py
```
The results will be displayed as a bar chart. This script analyzes the latest file by default. In order to analyze other files, you need to change the script and pass a list of file names as the argument.

# How to configure a simulation
- More information will be here soon.
- Read the comments in the __init__method of the Simulation class in [simulation.py](covid19_sim/simulation.py) file to learn more about the parameters.
- There will be more sample scenarios in the future.
- Do NOT use TestingStrategy.SCAN and TestingStrategy.SELF_CHECK as they are not yet well tested.

# Assumptions, rules and definitions
- Each day, each person interacts with up to 8 nearby contacts unless they are in quarantine. Each person's interactions happen one after another but each person is selected randomly from the population.
- If long distance travel is enabled, locations of some people are swapped with each other to simulate long distance traveling.
- There is a probability for an infected person to infect another person in each interaction.
- By default, 20% of all cases are assumed to be asymptomatic.
- Deaths and fatality of the disease are ignored. It is assumed that everything else being equal, minimizing the number of infected people and maximum infection rate will minimize deaths.
- Infected people recover from the disease between 14 and 21 days. Both symptomatic and asymptomatic cases become recovered/immune when they recover.
- When a symptomatically infected person recovers, the person is 'identified as immune' by the society. However, when an asymptomatic case recovers, they are NOT 'identified as immune' because there is no way to know whether or not they had the disease.
- Right now, when a person is tested positive but has no symptoms (asymptomatic case), they are NEVER identified as immune by the society even if they are removed from quarantine.
- Contact tracing
    - There is a probability for a symptomatically infected person to go to hospital. There is also a probability for uninfected people to go to hospital with similar symptoms.
    - There is a probability to trace the contacts of a person. If the probability is lower than 100%, some of the contacts will not be traced at all.
    - There is a probability to reach the traced contacts of a person on a given day. A traced contact may or may not be available for testing on a given day.
    - There is a probability for a traced contact to be put in quarantine. When a contact is traced, it may or may not be possible (or desirable) to put the contact in quarantine until testing takes place.
    - If an infected person does long distance travel and goes to hospital, the person's contacts before the travel are not traced.
    - When a person goes to hospital, they are put in a sick queue. Their contacts are put in a contacts queue.
    - People in the sick queue are tested first. If the queue is empty, people in the contacts queue are tested. If a contact is positive, their contacts are traced.
    - Everyday, there is a maximum testing capacity. If the number of people in the queues are larger than the testing capacity, some people may not be tested on that day.
    - If scanning option is on, all the people who are not in any test queues and who are not identified as immune are put in a screening queue. Initially, people are ordered starting with the top left corner of the map, going downwards row by row.
    - The screening queue is an ordered list. People who are not tested and not in quarantine have higher priority. Previously tested people are also put in the queue after a certain time that can be specified as a parameter. Those whose last test dates are further in the past are given priority.
    - People in the screening queue are tested only if there is excess testing capacity on a given day.
- People who are infected 1 day ago are counted as new infections. Obviously, this information is not visible to the society, but the assumption is that no matter how society counts new infections, that number will be proportional to the real number of new infections.
- Unnecessarily quarantined people
    - When an uninfected person is quarantined, this is counted as an 'unnecessary quarantine'. Although it can be argued that isolating uninfected people is not unnecessary because it prevents them from catching the disease and later infecting others, this term makes sense from an ideal viewpoint: if the society could perfectly detect the infected people, they would be put in quarantine without the uninfected people having to give up their freedoms. So, it is important to keep in mind that 'unnecessary quarantine' is defined from that viewpoint. This metric measures the ineffectiveness of the society to detect the infected people and chase the virus. Therefore, the smaller is this metric, the better is the strategy in finding and isolating infected people.
- Unnecessary quarantine days
    - This is the total number of days lost to unnecessary quarantine. The larger is this metric, the more uninfected people lose their freedoms because of quarantines.
- By experimenting with different parameters and rules, it is possible to optimize the outcome and reduce infection rates and days lost to unnecessary quarantine.

# License
GNU General Public License v3.0 or later

See [LICENSE](LICENSE) to see the full text.