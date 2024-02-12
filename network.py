import os
import networkx as nx
import matplotlib.pyplot as plt

# Directory for workplaces
directory = './data/links/network/workplaces_loop'

# Links start on day 3
start_day = 3
end_day = 30

# HIGH-LEVEL OVERVIEW
# Before anything: create viral load profile for every person in population for if they get infected
# Iteratively:
    # Create rooms -- describe two people that interact with each other on a given day
    # Create viral load list
    # Call aerosol model on room with people meeting (includes parameters like duration)
        # Area of the room comes after each day since we need to know how many people are in the room
    # Aerosol model infects rooms

# Parse the data
def parse_data(data):
    parsed_data = []
    for line in data[:1000]: # Currently only going through first 1000 lines
        parsed_data.append((int(line[0]), int(line[1]), float(line[2])))
    return parsed_data

# Simulate rooms
def sim_rooms(parsed_data, day):
    rooms = {}
    for link in parsed_data:
        individual1, individual2, _ = link
        if individual1 not in rooms:
            rooms[individual1] = set()
        if individual2 not in rooms:
            rooms[individual2] = set()
        rooms[individual1].add(individual2)
        rooms[individual2].add(individual1)
    return rooms

# Visualize
def graph(rooms):
    G = nx.Graph()
    for individual, connections in rooms.items():
        for connection in connections:
            G.add_edge(individual, connection)

    pos = nx.spring_layout(G)  # positions for all nodes

    # Nodes
    nx.draw_networkx_nodes(G, pos, node_size=50)
    # Edges
    nx.draw_networkx_edges(G, pos, width=5)
    # Labels
    nx.draw_networkx_labels(G, pos, font_size=5)

    plt.axis("off")
    plt.show()
    print('Finished graphing')

# Proccess links file
def process_links_file(file_name, day):
    # Open the file and read the data
    with open(file_name, 'r') as file:
        data = [line.split() for line in file.readlines()]
        # Parse the data
        parsed_data = parse_data(data)
        # Simulate rooms
        rooms = sim_rooms(parsed_data, day)
        # Visualize rooms
        print('Graphing rooms, please wait...')
        graph(rooms)

# TEMP: run only start day
file_name = os.path.join(directory, f'links_{start_day}.txt')
process_links_file(file_name, start_day)

# # Iterate through each day
# for day in range(start_day, end_day + 1):
#     file_name = os.path.join(directory, f'links_{day}.txt')
#     # Process the links file for the current day
#     process_links_file(file_name, day)

# Maybe use a class to keep track of agent IDs and connections; if so, build further
# Bottom line: want contact matrix
class Individual:
    def __init__(self, id):
        self.id = id
        self.connections = []

###