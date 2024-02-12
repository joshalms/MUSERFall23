import os

# Directory for workplaces
directory = '/data/links/network/workplaces_loop'

# Links start on day 3
start_day = 3
end_day = 30

# HIGH-LEVEL OVERVIEW
# Before anything: create viral load profile for every person in population for if they get infected
# Iteratively:
    # Create rooms (meetings) -- describe two people that interact with each other on a given day
    # Create viral load list
    # Call aerosol model on room with people meeting (includes parameters like duration)
        # Area of the room comes after each day since we need to know how many people are in the room
    # Aerosol model infects rooms

# Define functions for processing meetings
def process_links_file(file_name, day):
    # Open the file and read the data
    with open(file_name, 'r') as file:
        data = [line.split() for line in file.readlines()]
        # Parse the data
        parsed_data = parse_data(data)
        # Simulate daily meetings
        meetings = simulate_daily_meetings(parsed_data, day)
        # Visualize the meetings
        graph(meetings)

# Iterate through each day
for day in range(start_day, end_day + 1):
    file_name = os.path.join(directory, f'links_{day}.txt')
    # Process the links file for the current day
    process_links_file(file_name, day)

# Parse the data
def parse_data(data):
    parsed_data = []
    for line in data:
        parts = line.split()
        parsed_data.append((int(parts[0]), int(parts[1]), float(parts[2])))
    return parsed_data

# Maybe use a class to keep track of agent IDs and connections; if so, build further
# Bottom line: want contact matrix
class Individual:
    def __init__(self, id):
        self.id = id
        self.connections = []

# Daily Meeting Simulation
def simulate_daily_meetings(parsed_data, day):
    meetings = {}
    for link in parsed_data:
        individual1, individual2, _ = link
        if individual1 not in meetings:
            meetings[individual1] = set()
        if individual2 not in meetings:
            meetings[individual2] = set()
        meetings[individual1].add(individual2)
        meetings[individual2].add(individual1)
    return meetings

# Visualize
def graph(meetings):
    pass  # Code to graph meetings using TBD

###