from enum import Enum

class Status(Enum):
    SUSCEPTIBLE = 1
    INFECTED = 2
    RECOVERED_IMMUNE = 3

class InfectionType(Enum):
    ASYMPTOMATIC = 1
    SYMPTOMATIC = 2

class TestQueueType(Enum):
    SICK_QUEUE = 1
    CONTACTS_QUEUE = 2
    SCREENING_QUEUE = 3

class Person:
    def __init__(self, **kwargs):
        self.x = None
        self.y = None

        self.status = kwargs.get('status', Status.SUSCEPTIBLE)
        self.infection_type = kwargs.get('infection_type', None)
        self.days_since_infected = 0
        self.in_quarantine = False

        # Days in 'uninterrupted' quarantine. If person leaves
        # quarantine, this must be set to 0.
        self.days_in_quarantine = 0

        self.tested = False
        self.last_test_date = -999

        # most recent test result
        self.tested_positive = False

        # correctly identified as immune: person is immune and it is known
        self.identified_immune = False

        # falsely identified as immune: person is NOT immune but
        # it is believed to be immune
        self.falsely_identified_immune = False

        # True means person was quarantined despite not being infected
        self.was_unnecessarily_quarantined = False

        # True means person's contacts had been traced
        self.contacts_traced = False

        # the day when the person has been traced as a contact
        # of an infected person
        self.traced_on = None

        # used for contact tracing. Which test queue
        # is the person in?
        self.test_queue = None
        
    def remove_from_quarantine(self):
        """Remove person from quarantine"""
        self.in_quarantine = False
        self.days_in_quarantine = 0