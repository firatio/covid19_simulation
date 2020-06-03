from covid19_sim.person import *

from random import seed
from random import randint
import random
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import colors
import math
import numpy as np
from textwrap import wrap

import statistics

def happens(probability):
    """compute if event with given probability happens"""
    if random.random() < probability:
        return True
    return False

class TestingStrategy(Enum):
    SCAN = 's'
    CONTACT_TRACE = 'ct'
    CONTACT_TRACE_WITH_SCANNING = 'ctws'
    SELF_CHECK = 'sc'


class Test:
    """Represents tests such as diagnostic tests, serological tests"""
    def __init__(self, **kwargs):
        self.__accuracy = kwargs.get('accuracy', 1)
        self.__false_positive_probability = kwargs.get('false_positive_probability', 0)

        # If True, the test is applied twice
        self.__use_double_test = kwargs.get('use_double_test', False)

        # Applies only if double test is used. If True, both tests
        # must be positive to interpret the result as positive.
        # By default, result is positive if one of two tests is
        # positive, which produces a higher accuracy but also a higher
        # false positive probability.
        self.__require_double_positive = kwargs.get('require_double_positive',
                                                        False)

        def calculate_effective_odds():
            false_negative_probability = 1 - self.__accuracy

            if self.__use_double_test:
                if not self.__require_double_positive:
                    return [1 - false_negative_probability**2,
                            1 - (1 - self.__false_positive_probability)**2]
                else:
                    return [self.__accuracy**2,
                            self.__false_positive_probability**2]
            
            # return single test results by default
            return [self.__accuracy, self.__false_positive_probability]

        effective_odds = calculate_effective_odds()
        self.__effective_accuracy = effective_odds[0]
        self.__effective_false_positive_probability = effective_odds[1]

    def get_accuracy(self):
        return self.__effective_accuracy

    def get_false_positive_probability(self):
        return self.__effective_false_positive_probability

    def get_properties(self):
        return [self.__accuracy, self.__false_positive_probability, self.__use_double_test, self.__require_double_positive]

class Simulation:
    """Simulation of how the state of a population evolves in time."""
    def __init__(self, **kwargs):
        self.label = kwargs.get('label', 'Simulation')
        self.population = []     # list of people
        # initial distribution of Susceptible, Infected
        # and Recovered/Immune as percentages
        self.initial_state = kwargs.get('initial_state', [0.99, 0.01, 0]) # S, I, R
        # if True, days since infection for infected people is a normal
        # distribution, otherwise it is 0 for all of them
        self.normalize_infected = kwargs.get('normalize_infected', True)
        
        # world size
        self.width = kwargs.get('width', 10)
        self.height = kwargs.get('height', 10)
        self.world = [ [ 0 for j in range(self.height) ] for i in range(self.width) ]

        self.max_days = 1000    # max simulation length
        self.current_day = 0

        # FOR GRAPHS AND STATISTICS
        # Stores the most recent state of the population
        # S-I-R-Q-UQ-T [Susceptible, Infected, Recovered, Quarantined,
        # Unnecessarily quarantined, number of tests]
        self.last_snapshot = [0, 0, 0, 0, 0, 0]
        self.historical_snapshots = []
        self.population_history = []
        self.generate_graph = False

        # Behavioral parameters
        #self.contacts_per_day = 8    # maximum number of contacts per day
        
        # DISEASE PARAMETERS
        # infection risk per interaction
        self.infection_risk = kwargs.get('infection_risk', 0.029)
        self.min_disease_duration = 14  # minimum days before recovery
        self.max_disease_duration = 21  # maximum days before recovery
        # percentage of infected people who show no symptoms while infected
        self.asymptomatic_ratio = kwargs.get('asymptomatic_ratio', 0.2)

        # PARAMETERS TO CONTROL THE OUTBREAK
        # If True, there will be a general quarantine
        self.general_quarantine = kwargs.get('general_quarantine', False)
        # what percentage of the population will be in the general quarantine?
        self.general_quarantine_percent = kwargs.get('general_quarantine_percent', 0)
        # if it is safe to do so, people will be removed from quarantine
        # after so many days
        self.quarantine_duration = kwargs.get('quarantine_duration', 9999)
        # People will be removed from quarantine after so many days if they
        # are tested negative
        self.quarantine_duration_after_negative_test = kwargs.get('quarantine_duration_after_negative_test', 1)
        # People will be removed from quarantine after so many days if they
        # are tested positive. This applies to certain scenarios in which
        # a testing strategy is not used. In such cases, this parameter can be 
        # used to remove them from quarantine after a certain period.
        self.quarantine_duration_after_positive_test = kwargs.get('quarantine_duration_after_positive_test', 9999)

        # todo: make Intervention a class with meaningful properties
        # [apply intervention?, number of days new infected increasing, number of days new infected decreasing]
        self.intervention = kwargs.get('intervention', [False, None, None])

        # long distance travel
        # NOTE: this will be modified in later versions!
        # [is_on?, traveler_percentage, min_distance_percentage]
        self.long_distance_travel = kwargs.get('long_distance_travel', [False, 0, 0])

        # air travel. based on a new model, different from long_distance_travel
        # [is_on?, traveler_percentage, min_distance_percentage]
        self.air_travel = kwargs.get('air_travel', [False, 0, 0])

        # TESTING STRATEGY PARAMETERS
        self.testing_strategy = kwargs.get('testing_strategy', None)
        # (used only in TestingStrategy.SCAN) how much of the testable will be screened?
        self.screening_percentage = kwargs.get('screening_percentage', 0)
        # maximum tests that can be run in a day, expressed as a 
        # percentage of population
        self.max_test_percentage = kwargs.get('max_test_percentage', 0)
        # how many days after the test will it be repeated?
        self.test_repeat_period = kwargs.get('test_repeat_period', 9999)
        # if True long distance travelers are tested before travel
        self.test_before_travel = kwargs.get('test_before_travel', False)
        # if True long distance travelers are tested after travel
        self.test_after_travel = kwargs.get('test_after_travel', False)
        # if True long distance travelers are isolated if tested positive
        self.isolate_travelers = kwargs.get('isolate_travelers', True)

        # CONTACT TRACING PARAMETERS
        # Some infected people with symptoms will go to hospitals
        # or test centers. ASSUMPTION: People will go to hospital with
        # different probabilities on different days.
        # [day 1, day 2, etc..]
        self.probabilities_to_go_to_hospital = kwargs.get(
            'probabilities_to_go_to_hospital',
            [0.1, 0.2, 0.3, 0.6, 0.8, 0.6, 0.3, 0.1]
        )
        # some contacts may not be traced at all. If
        # a person has 8 contacts, assume that only a certain
        # percentage of them can be traced.
        # for example, 0.5 means that on average only 50% of the contacts can be traced
        self.probability_to_be_traced = kwargs.get('probability_to_be_traced', 0.5)
        # when a person is identified for testing, there is a
        # probability that the person cannot be reached on a 
        # given day
        self.probability_to_be_reached = kwargs.get('probability_to_be_reached', 0.37)
        # Probability for anyone to have symptoms similar to the disease
        # without actually being infected. This parameter can also be used
        # to specify 'willingness of healthcare system to test those with mild 
        # symptoms'. If the system is willing to do so, this probability can 
        # be set higher so that more people go to hospitals. This probability
        # is also related to the nature of the disease. If the disease has
        # many symptoms that are similar to other diseases existing in the 
        # society, this probability should be higher.
        self.similar_symptoms_probability = kwargs.get('similar_symptoms_probability', 0)
        # probability to put contacts in quarantine when they are traced
        self.probability_to_put_contacts_in_quarantine = kwargs.get('probability_to_put_contacts_in_quarantine', 0.5)
        # highest priority. for people who go to hospital with symptoms
        self.sick_test_queue = []
        # 2nd priority. for people who are traced as contacts of sick people
        self.contacts_test_queue = []
        # lowest priority. for people who are to be screened
        self.screening_test_queue = []

        # DIAGNOSTIC TEST PARAMETERS
        self.diagnostic_test = kwargs.get('diagnostic_test', 
            Test(
                accuracy = 0.98, 
                false_positive_probability = 0.01,
                use_double_test = False, 
                require_double_positive = False
            )
        )
        self.test_accuracy = self.diagnostic_test.get_accuracy()
        self.false_positive_probability = self.diagnostic_test.get_false_positive_probability()
        
        # SEROLOGICAL TEST PARAMETERS
        self.serological_test = kwargs.get('serological_test', 
            Test(
                accuracy = 0.98, 
                false_positive_probability = 0.01,
                use_double_test = False, 
                require_double_positive = False
            )
        )
        self.serological_accuracy = self.serological_test.get_accuracy()
        self.serological_false_positive_probability = self.serological_test.get_false_positive_probability()

        # SCAN STRATEGY PARAMETERS
        self.apply_diagnostic_test = kwargs.get('apply_diagnostic_test', True)
        self.apply_serological_test = kwargs.get('apply_serological_test', False)
        # Which test will be applied first?
        self.apply_diagnostic_test_first = kwargs.get('apply_diagnostic_test_first', True)
        # sampling method to choose test subjects: if not random, 
        # it is ordered by last test date ascending
        self.use_random_sampling_method = kwargs.get('use_random_sampling_method', False)

        # FLAGS and CONTROL PARAMETERS
        self.general_quarantine_applied = False
        self.days_decreasing_new_infections = 0
        self.days_increasing_new_infections = 0

        # STATISTICAL PARAMETERS
        self.total_infected = 0
        self.max_infected = 0
        self.total_unnecessarily_quarantined = 0
        self.total_unnecessary_quarantine_days = 0
        self.total_identified_immune = 0
        self.total_falsely_assumed_immune = 0
        self.daily_tests = 0
        self.total_tests = 0
        self.total_tested_people = 0
        self.epidemic_duration = 0
        self.new_infections = []
        self.new_infections_ma = [] # moving average
        self.general_quarantines_applied = 0
        self.general_quarantines_lifted = 0
        
        # seed random number generator
        seed(datetime.now())

    def initialize(self):
        """Build initial population"""
        size = self.width * self.height
        infected = self.total_infected = int(self.initial_state[1] * size)
        asymptomatic = int(infected * self.asymptomatic_ratio)
        symptomatic = infected - asymptomatic
        recovered = int(self.initial_state[2] * size)
        susceptible = size - (infected + recovered)

        # susceptible
        for i in range(susceptible):
            self.population.append(Person(status = Status.SUSCEPTIBLE))

        # infected = asymptomatic + symptomatic
        # create a distribution for days_since_infected (dsi)
        dsi = []
        if self.normalize_infected:
            mu, sigma = 6, 3
            sample = np.random.normal(mu, sigma, infected)
            # convert to integers
            for val in sample:
                r = int(round(val))
                if r < 1:
                    r = 1
                if r > 12:
                    r = 12
                dsi.append(r)
        # create infected people
        for i in range(asymptomatic):
            p = Person(status = Status.INFECTED,
                       infection_type = InfectionType.ASYMPTOMATIC)
            if self.normalize_infected:
                p.days_since_infected = dsi.pop()
            self.population.append(p)
        for i in range(symptomatic):
            p = Person(status = Status.INFECTED,
                       infection_type = InfectionType.SYMPTOMATIC)
            if self.normalize_infected:
                p.days_since_infected = dsi.pop()
            self.population.append(p)

        # recovered/immune
        for i in range( recovered ):
            p = Person(status = Status.RECOVERED_IMMUNE)
            # todo: Later, make it optional to mark some
            # people as identified immune in the beginning. For
            # now, recovered/immune people will not be known to 
            # the society.
            #p.identified_immune = True
            #self.total_identified_immune += 1
            self.population.append(p)

        # randomly distribute susceptible, infected and recovered people
        random.shuffle(self.population)

        # if there is a general quarantine, put some people in quarantine
        if self.general_quarantine:
            self.apply_general_quarantine(self.general_quarantine_percent)

        # convert population to 2d array
        for x in range(self.width):
            for y in range(self.height):
                index = self.height*x + y
                p = self.population[index]
                p.x = x
                p.y = y
                self.world[x][y] = p

    def apply_general_quarantine(self, percentage, **kwargs):
        """Randomly put a percentage of the population in quarantine."""
        disable_long_distance_travel = kwargs.get('disable_long_distance_travel', False)
        target = random.sample(
            self.population,
            int(percentage * len(self.population))
            )
        for p in target:
            p.in_quarantine = True
            if p.status != Status.INFECTED:
                p.was_unnecessarily_quarantined = True
                self.total_unnecessarily_quarantined += 1
        
        self.general_quarantine_applied = True
        self.general_quarantines_applied += 1

        if disable_long_distance_travel:
            self.long_distance_travel[0] = False

    def lift_general_quarantine(self, count, **kwargs):
        """Randomly remove a certain number of people from quarantine."""
        enable_long_distance_travel = kwargs.get('enable_long_distance_travel', False)
        quarantined = []
        
        for p in self.population:
            # Those who are tested positive and not yet identified
            # as recovered are not removed
            if p.in_quarantine and not (p.tested_positive and not p.identified_immune):
                quarantined.append(p)
        if count > len(quarantined):
            target = quarantined
        else:
            target = random.sample(quarantined, count)
        #print(len(target))
        for p in target:
            p.remove_from_quarantine()

        self.general_quarantine_applied = False
        self.general_quarantines_lifted += 1

        if enable_long_distance_travel:
            self.long_distance_travel[0] = True

    def get_description(self):
        """generate a short description of the simulation based on parameters"""
        txt = '[' + str(self.width) + ', ' + str(self.height) + '] '

        if self.air_travel[0]:
            txt += 'Air travel enabled. '
        elif not self.long_distance_travel[0]:
            txt += 'Air travel disabled. '
        elif self.long_distance_travel[0]:
            txt += 'Long distance travel enabled. '
        
        if self.test_before_travel or self.test_after_travel:
            if self.test_before_travel:
                txt += 'Travelers tested before travel. '
            if self.test_after_travel:
                txt += 'Travelers tested after travel. '
            if self.isolate_travelers:
                txt += 'Travelers isolated if tested positive. '

        if self.general_quarantine:
            txt += 'General quarantine of ' + str(int(self.general_quarantine_percent * 100)) + '%. '

        if self.quarantine_duration >= 9999:
            txt += 'Indefinite quarantine. '
        else:
            txt += 'Quarantine duration: ' + str(self.quarantine_duration) + ' days. '

        if self.testing_strategy == TestingStrategy.CONTACT_TRACE:
            txt += 'Testing strategy: Contact trace. '
        elif self.testing_strategy == TestingStrategy.CONTACT_TRACE_WITH_SCANNING:
            txt += 'Testing strategy: Contact trace with scanning. '
        elif self.testing_strategy == TestingStrategy.SELF_CHECK:
            txt += 'Testing strategy: Self check. '
        elif self.testing_strategy == TestingStrategy.SCAN:
            txt += 'Testing strategy: Scan orderly. '

        if self.testing_strategy == TestingStrategy.CONTACT_TRACE \
            or self.testing_strategy == TestingStrategy.CONTACT_TRACE_WITH_SCANNING:
            txt += 'Max testing capacity: ' + str(round(self.max_test_percentage * 100, 2)) + '%. '
            txt += 'Probability to be traced: ' + str(int(self.probability_to_be_traced * 100)) + '%. '
            txt += 'Probability to be reached: ' + str(int(self.probability_to_be_reached * 100)) + '%. '
            txt += 'Probability to have similar symptoms: ' + str(int(self.similar_symptoms_probability * 100)) + '%. '
            if self.testing_strategy == TestingStrategy.CONTACT_TRACE_WITH_SCANNING:
                txt += 'Test repeat period: ' + str(int(self.test_repeat_period)) + '. '
            txt += 'Probability to put contacts in quarantine: ' \
                + str(self.probability_to_put_contacts_in_quarantine) + ' '

        if self.intervention[0]:
            txt += 'Apply quarantine after ' + str(self.intervention[1]) \
                    + ' days of increase. Lift quarantine after ' \
                    + str(self.intervention[2]) + ' days of decrease. '

        if self.screening_percentage > 0:
            txt += str(int(self.screening_percentage * 100)) + '% screening/day. '
            if self.test_repeat_period != 999:
                txt += 'Test repeat period: ' + str(int(self.test_repeat_period)) + '. '
            txt += 'Test randomly. ' if self.use_random_sampling_method else ''#'Test by last test date. '
            if self.quarantine_duration < 9999:
                txt += 'Quarantine duration: ' + str(self.quarantine_duration) + ' days. '

            if self.apply_diagnostic_test:
                txt += ('Test accuracy: '
                        + str(int(self.test_accuracy*100))
                        + '%. False positive prob: '
                        + str(int(self.false_positive_probability*100))
                        + '%. ')
                if self.diagnostic_test.use_double_test:
                    txt += 'Double diagnostic test. '
                    if self.diagnostic_test.require_double_positive:
                        txt += 'Double positive = positive. '
                    else:
                        txt += 'Single positive = positive. '
                else:
                    txt += 'Single diagnostic test. '
            
            if self.apply_serological_test:
                txt += 'Serological test.'

            if self.apply_diagnostic_test and self.apply_serological_test:
                if self.apply_diagnostic_test_first:
                    txt += 'Diagnostic test first.'
                else:
                    txt += 'Serological test first.'

        return txt

    def dump(self, p):
        """dump person's info"""
        print('________________')
        print( str(p.x) + ',' + str(p.y) )
        print(p.status)
        print(p.infection_type)
        print('days since infection: ' + str(p.days_since_infected))
        print('tested: ' + str(p.tested))
        print(p.last_test_date)
        print('tested positive: ' + str(p.tested_positive))
        print('in quarantine: ' + str(p.in_quarantine))

    def take_snapshot(self):
        """Take a daily snapshot of the population."""

        def compute_moving_average(time_series, day, moving_average_length):
            """Compute moving average up to the day in the time series
            with the given length. Day starts at 0."""
            if day < moving_average_length - 1:
                raise Exception('day is too small to compute moving average')
            if day > len(time_series) - 1:
                raise Exception('day is out of series')
            period = time_series[(day - moving_average_length + 1):(day + 1)]
            return sum(period)/len(period)
        
        # data to keep track of
        susceptible = 0
        infected = 0
        recovered = 0
        quarantined = 0
        unnecessarily_quarantined = 0
        new_infections = 0

        for p in self.population:
            if p.status == Status.SUSCEPTIBLE:
                susceptible += 1
            elif p.status == Status.INFECTED:
                infected += 1
                # How does society identify and count new infections?
                # ASSUMPTION: no matter what the method is, the number
                # will be proportional to the real new infections.
                # However, if identification relies on self-reporting,
                # like going to hospital for testing etc, the trend
                # may fluctuate and change with the average age of infection
                # in the society because probability to go to hospital
                # is not the same for each infection day!
                if p.days_since_infected == 1:
                    new_infections += 1
            else:
                recovered += 1
            if p.in_quarantine:
                quarantined += 1
                if p.status != Status.INFECTED:
                    unnecessarily_quarantined += 1
                    # record only first time unnecessary quarantine
                    if not p.was_unnecessarily_quarantined:
                        p.was_unnecessarily_quarantined = True
                        self.total_unnecessarily_quarantined += 1

        # record new infections
        self.new_infections.append(new_infections)

        # calculate 7-day moving averages for new infections and determine trend
        moving_average_length = 7
        if self.current_day >= moving_average_length - 1:
            ma = compute_moving_average(self.new_infections, self.current_day, moving_average_length)
            self.new_infections_ma.append(ma)
            # determine trend
            if self.current_day > moving_average_length - 1:
                difference = self.new_infections_ma[-1] - self.new_infections_ma[-2]
                # interpret 'difference = 0' differently depending on the trend
                if difference < 0 or (difference == 0 and self.days_decreasing_new_infections > 0):
                    self.days_decreasing_new_infections += 1
                    self.days_increasing_new_infections = 0 # reset
                elif difference > 0 or (difference == 0 and self.days_increasing_new_infections > 0):
                    self.days_increasing_new_infections += 1
                    self.days_decreasing_new_infections = 0 # reset

        # update max_infected if necessary
        if infected > self.max_infected:
            self.max_infected = infected

        self.last_snapshot[0] = susceptible
        self.last_snapshot[1] = infected
        self.last_snapshot[2] = recovered
        self.last_snapshot[3] = quarantined
        self.last_snapshot[4] = unnecessarily_quarantined
        self.last_snapshot[5] = self.daily_tests

        # keep track of total days lost because of unnecessary quarantine
        self.total_unnecessary_quarantine_days += unnecessarily_quarantined

        snapshot = self.last_snapshot.copy()
        self.historical_snapshots.append(snapshot)
       
        if self.generate_graph:
            # create a visual representation of the 2d world
            visual = [ [ 0 for j in range(self.height) ] for i in range(self.width) ]
            for x in range(self.width):
                for y in range(self.height):
                    p = self.world[x][y]
                    if p.status == Status.SUSCEPTIBLE:
                        visual[x][y] = 0
                    elif p.status == Status.INFECTED:
                        visual[x][y] = 1
                    elif p.status == Status.RECOVERED_IMMUNE:
                        visual[x][y] = 2
            self.population_history.append(visual)

    def run_diagnostic_test(self, p, **kwargs):
        """run a diagnostic test"""
        isolate_if_positive = kwargs.get('isolate_if_positive', True)
        remove_from_quarantine_if_negative = kwargs.get('remove_from_quarantine_if_negative', True)
        
        self.daily_tests += 1
        self.total_tests += 1
        if not p.tested:
            self.total_tested_people += 1

        p.tested = True
        p.last_test_date = self.current_day

        positive = (
            # infected people will trigger a positive response
            # based on test accuracy
            (p.status == Status.INFECTED and happens(self.test_accuracy))
            # those who are not infected may trigger a false positive
            or (p.status != Status.INFECTED
                and happens(self.false_positive_probability)
                )
        )

        if positive:
            p.tested_positive = True
            if isolate_if_positive:
                p.in_quarantine = True
            return True
        else:
            if p.in_quarantine:
                # For now, if p has symptoms, p is not removed
                # from quarantine automatically even if the test
                # is negative.

                # If p has no symptoms, remove from quarantine 
                # if the options says so.
                if not (p.status == Status.INFECTED
                        and p.infection_type == InfectionType.SYMPTOMATIC):
                    if remove_from_quarantine_if_negative:
                        p.in_quarantine = False
                        p.days_in_quarantine = 0  # reset quarantine count
            return False

    def run_tests(self):
        """Every day run tests to make decisions"""
        day = self.current_day
        self.daily_tests = 0

        def is_recovered_immune(p):
            """run a serological test"""
            if p.status == Status.RECOVERED_IMMUNE and happens(self.serological_accuracy):
                if not p.identified_immune:
                    self.total_identified_immune += 1
                    p.identified_immune = True
                p.in_quarantine = False
                p.days_in_quarantine = 0
                return True
            elif happens(self.serological_false_positive_probability):
                # false positive, although p is susceptible or infected,
                # serological test says p is immune
                p.falsely_identified_immune = True
                p.in_quarantine = False
                p.days_in_quarantine = 0
                self.total_falsely_assumed_immune += 1
                return True
            # negative in other cases
            return False

        def test_with_contact_tracing():
            """Test people who have been in contact with the positively
            tested people using contact tracing. Optionally, scan
            people periodically for infection."""

            # General method:
            # - Use test_queues as the ordered list of people to be tested
            # - self.sick_test_queue has the highest priority. If it
            # gets empty, test self.contacts_test_queue.
            # - If both queues get empty, and strategy involves scanning,
            # randomly or orderly (closest to origin) test those
            # who have not been tested yet,
            # giving priority to those not in quarantine. If the person is
            # positive, add their contacts to contacts test queue and continue
            # with the contacts test queue. If everyone has
            # been tested, test those who have not been confirmed
            # as recovered/immune and need to be retested periodically, 
            # and test until capacity is fulfilled.

            def get_probability_to_go_to_hospital(days_since_infection):
                """Give the probability of symptomatically infected 
                person to go to hospital. Probability becomes 0 after
                a certain time."""
                max_days = len(self.probabilities_to_go_to_hospital)
                if days_since_infection <= max_days:
                    probability = self.probabilities_to_go_to_hospital[days_since_infection - 1]
                else:
                    probability = 0
                return probability

            def goes_to_hospital(p):
                """Decide if a person will go to hospital for testing"""
                symptomatically_infected = (
                    p.status == Status.INFECTED
                    and p.infection_type == InfectionType.SYMPTOMATIC
                )
                # Some infected with symptoms and some not infected
                # with similar symptoms will go to hospitals
                # or test centers.
                return (
                    (
                    symptomatically_infected and
                    happens(get_probability_to_go_to_hospital(p.days_since_infected))
                    )
                    or
                    (
                    not symptomatically_infected
                    and happens(self.similar_symptoms_probability)
                    )
                )

            def remove_from_test_queue(test_queue, p):
                """remove person from the given test queue"""
                for i in range(len(test_queue)):
                    w = test_queue[i]
                    if w.x == p.x and w.y == p.y:
                        del test_queue[i]
                        break

            def update_test_queues():
                """create an ordered list of those who should be tested"""

                testable = []
                # skip those identified as immune and tested
                # and those who are already in the sick test queue
                for p in self.population:
                    if not (p.identified_immune and p.tested) \
                        and p.test_queue != TestQueueType.SICK_QUEUE:
                        testable.append(p)

                i = 0
                testable_size = len(testable)
                while i < testable_size:
                    p = testable[i]
                    if goes_to_hospital(p):
                        # If p goes for testing, remove them from
                        # testable and add them to sick_queue.
                        #print(p.days_since_infected)
                        # if p is already in other test queues, remove
                        if p.test_queue == TestQueueType.CONTACTS_QUEUE:
                            remove_from_test_queue(self.contacts_test_queue, p)
                        elif p.test_queue == TestQueueType.SCREENING_QUEUE:
                            remove_from_test_queue(self.screening_test_queue, p)

                        # add to sick test queue and put p in quarantine
                        # because if testing capacity is low, p may not
                        # be tested immediately
                        self.sick_test_queue.append(p)
                        p.test_queue = TestQueueType.SICK_QUEUE
                        p.in_quarantine = True
                        del testable[i]
                        testable_size -= 1
                    else:
                        i += 1
                        continue

                # if the strategy also includes scanning, add some to
                # screening test queue
                if self.testing_strategy == TestingStrategy.CONTACT_TRACE_WITH_SCANNING:
                    untested = [x for x in testable if not x.tested \
                        and x.test_queue == None]
                    tested = [x for x in testable if x.tested \
                        and x.test_queue == None]

                    # priority to those who are not in quarantine and are
                    # not tested yet
                    untested.sort(key=lambda x: x.in_quarantine, reverse=False)
                    for p in untested:
                        self.screening_test_queue.append(p)
                        p.test_queue = TestQueueType.SCREENING_QUEUE
                    
                    # for tested, priority to those who were tested
                    # further back in time
                    tested.sort(key=lambda x: x.last_test_date, reverse=False)
                    for p in tested:
                        # if tested recently and was negative, continue with next one
                        if not p.tested_positive and day - p.last_test_date < self.test_repeat_period:
                            continue
                        self.screening_test_queue.append(p)
                        p.test_queue = TestQueueType.SCREENING_QUEUE

            def trace_contacts(p):
                """Trace p's contacts and add them to contacts test queue
                if necessary"""

                def consider_for_test(c):
                    """Decide whether the contact c should be tested"""
                    if c.test_queue == TestQueueType.SICK_QUEUE \
                        or c.test_queue == TestQueueType.CONTACTS_QUEUE \
                        or (c.tested and c.identified_immune):
                        return

                    # There is a probability for the contact
                    # to be traced. If traced, add contact to test queue.
                    if happens(self.probability_to_be_traced):
                        c.traced_on = self.current_day

                        # remove from screening test queue if necessary
                        if c.test_queue == TestQueueType.SCREENING_QUEUE:
                            remove_from_test_queue(self.screening_test_queue, c)

                        self.contacts_test_queue.append(c)
                        c.test_queue = TestQueueType.CONTACTS_QUEUE

                        # If the contact can be informed and put to 
                        # self-isolation, this strategy works much better
                        # but successful tracing does not imply that it is
                        # possible to talk with the contact. So, there
                        # is a probability for this to happen.
                        if happens(self.probability_to_put_contacts_in_quarantine):
                            c.in_quarantine = True

                # if p's contacts are already traced, skip
                if p.contacts_traced:
                    return

                # consider p's contacts for testing
                x = p.x
                y = p.y
                
                # trace left side contacts
                if x > 0:
                    consider_for_test(self.world[x-1][y])
                    if y > 0:
                        consider_for_test(self.world[x-1][y-1])
                    if y < self.height-1:
                        consider_for_test(self.world[x-1][y+1])
                # trace right side contacts
                if x < self.width-1:
                    consider_for_test(self.world[x+1][y])
                    if y > 0:
                        consider_for_test(self.world[x+1][y-1])
                    if y < self.height-1:
                        consider_for_test(self.world[x+1][y+1])
                # trace north contact
                if y > 0:
                    consider_for_test(self.world[x][y-1])
                # trace south contact
                if y < self.height-1:
                    consider_for_test(self.world[x][y+1])
                
                p.contacts_traced = True

            def run_daily_tests():
                """Run the tests and perform contact tracing"""
                number_of_tests = 0
                max_tests = math.floor(len(self.population) * self.max_test_percentage)
                # list for people who could not be reached for testing
                unreachable = []

                update_test_queues()

                while number_of_tests < max_tests:
                    # first test sick people
                    if self.sick_test_queue:
                        p = self.sick_test_queue.pop(0)
                    elif self.contacts_test_queue:
                        p = self.contacts_test_queue.pop(0)
                        # GUARD: if p was traced and never tested, it is 
                        # still possible that p cannot be reached for testing
                        if p.traced_on is not None:
                            if not p.tested:
                                if not happens(self.probability_to_be_reached):
                                    unreachable.append(p)
                                    continue
                    elif self.screening_test_queue:
                        p = self.screening_test_queue.pop(0)
                    else:
                        # if all queues are empty, stop
                        break
                   
                    # remove from whatever queue p was in
                    p.test_queue = None
                    # if p was tested positive, do not test again but trace contacts
                    if p.tested_positive:
                        trace_contacts(p)
                        continue
                    # in all other cases, test and if positive, trace contacts 
                    else:
                        number_of_tests += 1
                        if self.run_diagnostic_test(p):
                            trace_contacts(p)
                
                # add unreachable contacts to the beginning
                self.contacts_test_queue = unreachable + self.contacts_test_queue

            run_daily_tests()
            
        def test_with_scan_strategy():
            """
            UNTESTED!
            Test by scanning the entire population in an orderly or random fashion.
            """

            def generate_test_queue():
                """ generate a list of people to be tested: those ...
                who are not identified as recovered/immune and tested
                who are not tested or whose tests need to be repeated
                who are not tested positive with symptoms
                """
                testable = []
                for p in self.population:
                    test_criterion = not (p.identified_immune and p.tested) \
                        and (not p.tested or (p.tested and p.last_test_date <= day-self.test_repeat_period)) \
                        and not (p.tested_positive and p.status == Status.INFECTED and p.infection_type == InfectionType.SYMPTOMATIC)
                    if test_criterion:
                        testable.append(p)
                
                max_tests = math.floor(len(self.population) * self.screening_percentage)
                sample_size = max_tests if len(testable) > max_tests else len(testable)

                test_queue = []
                if self.use_random_sampling_method:
                    test_queue = random.sample(testable, sample_size)
                else:
                    # everyone must get tested first
                    untested = [x for x in testable if not x.tested]
                    #print( str(len(testable)) + ':' + str(len(untested)) )
                    tested = [x for x in testable if x.tested]

                    # priority to those who are not in quarantine
                    untested.sort(key=lambda x: x.in_quarantine, reverse=False)
                    for p in untested:
                        test_queue.append(p)
                    if len(test_queue) > sample_size:
                        test_queue = test_queue[:sample_size]
                    else:
                        difference = sample_size - len(test_queue)
                        tested.sort(key=lambda x: x.last_test_date, reverse=False)
                        for i in range(difference):
                            test_queue.append(tested[i])
                
                return test_queue

            #print('day: ' + str(day) + ' - tested today.....')
            number_of_tests = 0
            test_queue = generate_test_queue()
            #for p in sample:
            for p in test_queue:
                #self.dump(p)
                #print('in quarantine: '+str(p.in_quarantine))
                if not p.tested:
                    p.tested = True
                    self.total_tested_people += 1
                p.last_test_date = day
                number_of_tests += 1
                #self.total_tests += 1
                
                # which test will be applied first?
                if self.apply_diagnostic_test_first:
                    if self.apply_diagnostic_test:
                        if self.run_diagnostic_test(p):
                            continue
                    if self.apply_serological_test:
                        if is_recovered_immune(p):
                            continue
                else:
                    if self.apply_serological_test:
                        if is_recovered_immune(p):
                            continue
                    if self.apply_diagnostic_test:
                        if self.run_diagnostic_test(p):
                            continue

            #print('total tested people: ' + str(self.total_tested_people))
            #self.last_snapshot[5] = number_of_tests
            self.daily_tests += number_of_tests

        def test_with_self_check():
            """
            UNTESTED!
            No testing. People check themselves and self-isolate if
            they think they have the symptoms.
            """

            for p in self.population:
                # They will have the symptoms on different days and
                # self-isolate with different probabilities on
                # different days.
                probabilities_to_self_isolate = [0.3, 0.4, 0.5, 0.6, 0.8]
                if p.days_since_infected < 6:
                    probability_to_self_isolate = probabilities_to_self_isolate[p.days_since_infected-1]
                else:
                    probability_to_self_isolate = 1
                # self-isolation for 14 days
                if p.status == Status.INFECTED and InfectionType.SYMPTOMATIC:
                    if not p.in_quarantine:
                        if p.days_since_infected < 14 and happens(probability_to_self_isolate):
                            p.in_quarantine = True
                    # but 14 days later, p may lift self-isolation without being immune
                    elif p.days_since_infected >= 14:
                            p.in_quarantine = False
                            p.days_in_quarantine = 0
                elif p.status != Status.INFECTED and not p.identified_immune:
                    # false positive: some will self-isolate unnecessarily
                    false_positive_probability = 0.01
                    if happens(false_positive_probability):
                        p.in_quarantine = True

        # run tests using the specified strategy
        if (self.testing_strategy == TestingStrategy.CONTACT_TRACE
            or self.testing_strategy == TestingStrategy.CONTACT_TRACE_WITH_SCANNING):
            test_with_contact_tracing()
        elif self.testing_strategy == TestingStrategy.SELF_CHECK:
            test_with_self_check()
        else:
            test_with_scan_strategy()

    def swap_people(self, p1, p2):
            """Swap p1 and p2's locations."""
            p1_x = p1.x
            p1_y = p1.y
            # put p1 to p2's location
            p1.x = p2.x
            p1.y = p2.y
            self.world[p1.x][p1.y] = p1
            # put p2 to p1's location
            p2.x = p1_x
            p2.y = p1_y
            self.world[p2.x][p2.y] = p2

    def simulate_air_travel(self):
        """Simulates air travel in the society by selecting n people 
        not in quarantine, swapping their places with others who are not 
        in the list and are at least D distance away."""
        traveler_percentage = self.air_travel[1]
        min_distance_percentage = self.air_travel[2]
        to_be_moved = []
        # number of air travelers
        n = int(traveler_percentage*len(self.population))

        not_in_quarantine = [p for p in self.population if not p.in_quarantine]
        free = len(not_in_quarantine)
        if n*2 > free:
            n = math.floor(free/2)

        # randomly select some as to be moved
        random.shuffle(not_in_quarantine)
        for i in range(n):
            to_be_moved.append(not_in_quarantine.pop(0))

        # min distance is min distance percentage of world's diagonal length
        min_distance_square = min_distance_percentage*(self.width**2 + self.height**2)

        sim_population_size = self.width * self.height
        # for now, probability that a passenger is infected is equal to 
        # the percentage of infected in the general population
        # Later: 
        # -- this probability should be based on the local percentage 
        # of infected people.
        probability_passenger_is_infected = self.last_snapshot[1] / sim_population_size
        probability_passenger_is_recovered = self.last_snapshot[2] / sim_population_size

        # testing will reduce number of infected passengers boarding a plane
        if self.test_before_travel:
            probability_passenger_is_infected = probability_passenger_is_infected * (1 - self.diagnostic_test.get_accuracy())

        # each infected passenger will be assumed to infect 8 other susceptible
        # passengers during the flight, limited by susceptible passengers
        # Later:
        # -- infection during flight should be better modeled
        susceptible_passengers_ratio = 1 - (probability_passenger_is_infected + probability_passenger_is_recovered)
        newly_infected_ratio = probability_passenger_is_infected * 8
        if newly_infected_ratio > susceptible_passengers_ratio:
            newly_infected_ratio = susceptible_passengers_ratio
        probability_to_get_infected = newly_infected_ratio

        #print(self.current_day, flush=True)
        #print(total_flights_sim, flush=True)
        #print(probability_to_get_infected, flush=True)
        
        while to_be_moved:
            p1 = to_be_moved.pop(0)

            # if travelers are tested before travel and p1 is positive, 
            # continue with next traveler
            if self.test_before_travel and self.run_diagnostic_test(p1, isolate_if_positive = self.isolate_travelers):
                continue

            # find another traveler to swap with
            for i in range(len(not_in_quarantine)):
                p2 = not_in_quarantine[i]
                # if distance is greater than minimum distance,
                # swap their locations and remove p2 from list
                if (p1.x - p2.x)**2 + (p1.y - p2.y)**2 >= min_distance_square:
                    # if travelers are tested before travel and p1 is positive, 
                    # continue with next traveler
                    if self.test_before_travel and self.run_diagnostic_test(p2, isolate_if_positive = self.isolate_travelers):
                        continue

                    # If a traveler is infected, traveler remains 
                    # infected. If the traveler is susceptible,
                    # there is a risk for the traveler to get infected during
                    # air travel. If the risk happens, traveler will be 
                    # infected at the new location. 
                    travelers = [p1, p2]
                    for p in travelers:
                        if p.status == Status.SUSCEPTIBLE:
                            if happens(probability_to_get_infected):
                                #print('infected during flight', flush=True)
                                p.status = Status.INFECTED
                                if happens(self.asymptomatic_ratio):
                                    p.infection_type == InfectionType.ASYMPTOMATIC
                    
                    self.swap_people(p1, p2)
                    del not_in_quarantine[i]
                    
                    # if travelers are tested after travel, test them
                    if self.test_after_travel:
                        self.run_diagnostic_test(p1, isolate_if_positive = self.isolate_travelers)
                        self.run_diagnostic_test(p2, isolate_if_positive = self.isolate_travelers)
                    
                    # if a traveler is found, end the search
                    break

    def simulate_interactions(self):
        """ Compute how infection propagates in the population.
        Take each person and compute whether they infect others or get
        infected by them. Each person will interact with n neighboring
        people (close proximity). Optionally, some people will travel
        to other locations."""
        
        def interact(person, contact):
            """Simulate an interaction between 2 people and decide
            whether one infects the other."""

            # neither of them must be in quarantine
            if person.in_quarantine or contact.in_quarantine:
                return
            
            # for a new infection to happen, only one of them must be infected
            infection_can_happen = (
                (person.status == Status.INFECTED and contact.status == Status.SUSCEPTIBLE) or
                (person.status == Status.SUSCEPTIBLE and contact.status == Status.INFECTED)
            )
            if not infection_can_happen:
                return
            
            # one of them is at risk
            if person.status == Status.INFECTED:
                person_at_risk = contact
            else:
                person_at_risk = person

            # there is a probability for transmission
            if happens(self.infection_risk):
                person_at_risk.status = Status.INFECTED
                # will the contact be symptomatic or asymptomatic?
                if happens(self.asymptomatic_ratio):
                    person_at_risk.infection_type = InfectionType.ASYMPTOMATIC
                else:
                    person_at_risk.infection_type = InfectionType.SYMPTOMATIC
                self.total_infected += 1

        def swap_people(p1, p2):
            """Swap p1 and p2's locations."""
            p1_x = p1.x
            p1_y = p1.y
            # put p1 to p2's location
            p1.x = p2.x
            p1.y = p2.y
            self.world[p1.x][p1.y] = p1
            # put p2 to p1's location
            p2.x = p1_x
            p2.y = p1_y
            self.world[p2.x][p2.y] = p2

        def simulate_long_distance_travel():
            """Select n people not in quarantine,
            swap their places with others who are not in the list and
            are at least D distance away"""
            traveler_percentage = self.long_distance_travel[1]
            min_distance_percentage = self.long_distance_travel[2]
            to_be_moved = []
            # number of long distance travelers
            n = int(traveler_percentage*len(self.population))

            not_in_quarantine = [p for p in self.population if not p.in_quarantine]
            free = len(not_in_quarantine)
            if n*2 > free:
                n = math.floor(free/2)

            # randomly select some as to be moved
            random.shuffle(not_in_quarantine)
            for i in range(n):
                to_be_moved.append(not_in_quarantine.pop(0))

            # min distance is min distance percentage of world's diagonal length
            min_distance_square = min_distance_percentage*(self.width**2 + self.height**2)

            while to_be_moved:
                p1 = to_be_moved.pop(0)

                # if travelers are tested before travel and p1 is positive, 
                # continue with next traveler
                if self.test_before_travel and self.run_diagnostic_test(p1, isolate_if_positive = self.isolate_travelers):
                    continue

                # find another traveler to swap with
                for i in range(len(not_in_quarantine)):
                    p2 = not_in_quarantine[i]
                    # if distance is greater than minimum distance,
                    # swap their locations and remove p2 from list
                    if (p1.x - p2.x)**2 + (p1.y - p2.y)**2 >= min_distance_square:
                        # if travelers are tested before travel and p1 is positive, 
                        # continue with next traveler
                        if self.test_before_travel and self.run_diagnostic_test(p2, isolate_if_positive = self.isolate_travelers):
                            continue

                        swap_people(p1, p2)
                        del not_in_quarantine[i]
                        
                        # if travelers are tested after travel, test them
                        if self.test_after_travel:
                            self.run_diagnostic_test(p1, isolate_if_positive = self.isolate_travelers)
                            self.run_diagnostic_test(p2, isolate_if_positive = self.isolate_travelers)
                        # if a traveler is found, end the search
                        break
        
        def close_interactions():
            """Simulate interactions of each person in a randomized way."""
            shuffled_population = self.population.copy()
            random.shuffle( shuffled_population )
            for p in shuffled_population:
                x = p.x
                y = p.y

                # interact with left side
                if x > 0:
                    interact( p, self.world[x-1][y] )
                    if y > 0:
                        interact( p, self.world[x-1][y-1] )
                    if y < self.height-1:
                        interact( p, self.world[x-1][y+1] )
                # interact with right side
                if x < self.width-1:
                    interact( p, self.world[x+1][y] )
                    if y > 0:
                        interact( p, self.world[x+1][y-1] )
                    if y < self.height-1:
                        interact( p, self.world[x+1][y+1] )
                # interact with south if there's someone
                if y > 0:
                    interact( p, self.world[x][y-1] )
                # interact with north if there's someone
                if y < self.height-1:
                    interact( p, self.world[x][y+1] )

        # simulate travel and interactions
        # air travel has priority over long distance travel
        # long distance travel and air travel are mutually exclusive
        if self.air_travel[0]:
            self.simulate_air_travel()
        elif self.long_distance_travel[0]:
            simulate_long_distance_travel()
        close_interactions()
    
    def progress_disease(self):
        """Simulate the progression of disease and decide who will 
        recover on that day."""
        # probability to recover after minimum disease duration
        recovery_probability = 0.7

        for p in self.population:
            if p.status == Status.INFECTED:
                p.days_since_infected += 1
            
                # Infected may recover after min disease duration. Infected 
                # recovers after max disease duration.
                if (p.days_since_infected > self.min_disease_duration and happens(recovery_probability)) \
                    or p.days_since_infected > self.max_disease_duration:
                    p.status = Status.RECOVERED_IMMUNE
                    
                    # ASSUMPTION: Symptomatic people can be identified as 
                    # Recovered/Immune when they recover. Both the patient
                    # and others can know that patient has recovered, so
                    # the person can be removed from quarantine.
                    identify_as_recovered = (p.infection_type == InfectionType.SYMPTOMATIC)
                    if identify_as_recovered:
                        if not p.identified_immune:
                            p.identified_immune = True
                            self.total_identified_immune += 1
                        p.remove_from_quarantine()

                    # ASSUMPTION: Asymptomatic people may not be identified
                    # as recovered/immune and removed from quarantine even
                    # if they recover. This must happen somewhere else.

    def progress_quarantine(self):
        """Remove some from quarantine if safe to do so."""
        day = self.current_day
        for p in self.population:
            if p.in_quarantine:
                # one more day has passed in quarantine
                p.days_in_quarantine += 1

                # RULES:
                # - If p does not show any symptoms, was tested negative
                # and some quarantine duration has passed, remove from
                # quarantine.
                # - If p does not show any symptoms but was tested positive,
                # p may still be removed from quarantine after some duration.
                # - If p has no symptoms and was not tested, remove from
                # quarantine if a specified quarantine duration has passed.
                # Note: Symptomatic people are removed from
                # quarantine when they recover.
                if p.infection_type != InfectionType.SYMPTOMATIC:
                    tested_negative = (
                        p.tested and not p.tested_positive
                        and day - p.last_test_date > self.quarantine_duration_after_negative_test
                    )
                    finished_quarantine = (
                        p.days_in_quarantine > self.quarantine_duration
                        or (
                            p.tested and p.tested_positive
                            and day - p.last_test_date > self.quarantine_duration_after_positive_test
                        )
                    )
                    if tested_negative or finished_quarantine:
                        p.remove_from_quarantine()

    def run_day(self):
        """Simulate what happens during a day."""
        self.run_tests()
        self.simulate_interactions()
        self.progress_disease()
        self.progress_quarantine()

    def perform_intervention(self):
        """If conditions are satisfied, perform the intervention."""
        if self.intervention[0]:
            apply_quarantine_condition = self.intervention[1]
            lift_quarantine_condition = self.intervention[2]

            if apply_quarantine_condition is not None:
                if not self.general_quarantine_applied:
                    if self.days_increasing_new_infections >= apply_quarantine_condition:
                        self.apply_general_quarantine(self.general_quarantine_percent)

            if lift_quarantine_condition is not None:
                if self.general_quarantine_applied:
                    if self.days_decreasing_new_infections >= lift_quarantine_condition:
                        self.lift_general_quarantine(
                            int(self.general_quarantine_percent * len(self.population))
                            )

    def create_graphs(self):
        """Create graphs for the simulation"""
        fig, ax = plt.subplots()
        cmap = colors.ListedColormap(['blue', 'red', 'green'])
        bounds = [0, 0.5, 1.5, 2.5]
        norm = colors.BoundaryNorm(bounds, cmap.N)

        for i in range(len(self.population_history)):
            ax.cla()
            ax.imshow(self.population_history[i], cmap=cmap, norm=norm)
            ax.set_title("Day {}".format(i))
            plt.pause(0.1)

        susceptible = []
        infected = []
        recovered = []
        quarantined = []
        unnecessarily_quarantined = []
        number_of_tests = []

        txt = self.get_description()
        if self.intervention[0]:
            txt += (
                'Quarantines applied: ' + str(self.general_quarantines_applied) + '.'
                + ' Quarantines lifted: ' + str(self.general_quarantines_lifted) + '.'
            )
        txt = '\n'.join(wrap(txt, 150))

        for snapshot in self.historical_snapshots:
            susceptible.append(snapshot[0])
            infected.append(snapshot[1])
            recovered.append(snapshot[2])
            quarantined.append(snapshot[3])
            unnecessarily_quarantined.append(snapshot[4])
            number_of_tests.append(snapshot[5])
        
        fig, ax = plt.subplots()
        fig.text(.5, 0.90, txt, ha='center')
        ax.plot(susceptible, label = 'susceptible')
        ax.plot(infected, label = 'infected')
        ax.plot(recovered, label = 'recovered')
        ax.plot(quarantined, label = 'quarantined')
        ax.plot(unnecessarily_quarantined, label = 'unnecessarily quarantined')
        ax.plot(number_of_tests, label = 'number of tests')
        ax.set_xlabel('Day', fontsize = 14)

        plt.legend()
        plt.show()

    def run(self, generate_graph = False):
        """Run the simulation day by day."""
        self.generate_graph = generate_graph
        self.initialize()
        self.take_snapshot()

        for day in range(1, self.max_days):
            self.current_day = day

            self.run_day()
            self.take_snapshot()
            self.perform_intervention()

            # if no infected stop sim
            if self.last_snapshot[1] == 0:
                break

        self.epidemic_duration = self.current_day
        
        if generate_graph:
            self.create_graphs()