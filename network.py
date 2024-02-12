import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Venue:
    def __init__(self, venue_id, capacity, size, venue_type, transmission_characteristics):
        self.venue_id = venue_id
        self.capacity = capacity
        self.size = size
        self.venue_type = venue_type
        self.transmission_characteristics = transmission_characteristics

def calculate_weighted_transmission_probability(venue_type_i, venue_type_j, distance, duration, transmission_characteristics_i, transmission_characteristics_j):
    # Implement logic to calculate weighted transmission probability
    # Consider factors such as venue characteristics, distance, duration, etc.
    # Using roomProb.py

def simulate_transmission_within_venue(G, venue_type, sigma, rho, distance_threshold, duration_threshold):
    for edge in G.edges():
        i, j = edge
        venue_type_i = G.nodes[i]['venue'].venue_type
        venue_type_j = G.nodes[j]['venue'].venue_type
        transmission_characteristics_i = G.nodes[i]['venue'].transmission_characteristics
        transmission_characteristics_j = G.nodes[j]['venue'].transmission_characteristics
        distance = np.random.uniform(0, transmission_characteristics_i['max_distance'])  # Simulate distance
        duration = np.random.uniform(0, transmission_characteristics_i['max_duration'])  # Simulate duration

        # Check if the interaction occurs within the venue and meets the duration threshold
        if venue_type_i == venue_type and venue_type_j == venue_type and distance < distance_threshold and duration > duration_threshold:
            transmission_prob = calculate_weighted_transmission_probability(
                venue_type_i, venue_type_j, distance, duration, transmission_characteristics_i, transmission_characteristics_j
            )
            if np.random.rand() < transmission_prob:
                G[i][j]['transmission'] = True

def seir_simulation(G, sigma, rho, gamma, simulation_days, distance_threshold, duration_threshold):
    for day in range(simulation_days):
        # Simulate transmission within each venue
        for venue in venues:
            layer = venue.venue_type  # Layer corresponds to venue_type
            for edge in G.edges():
                i, j = edge
                if G.nodes[i]['venue'].venue_type == layer and G.nodes[j]['venue'].venue_type == layer:
                    simulate_transmission_within_venue(G, venue_type=layer, sigma=sigma, rho=rho, distance_threshold=distance_threshold, duration_threshold=duration_threshold)

        # Update individual states based on SEIR model
        for node in G.nodes():
            state = G.nodes[node]['state']
            if state == 'S':
                # Transition from Susceptible to Exposed
                if np.random.rand() < sigma:
                    G.nodes[node]['state'] = 'E'
            elif state == 'E':
                # Transition from Exposed to Infectious
                if np.random.rand() < rho:
                    G.nodes[node]['state'] = 'I'
            elif state == 'I':
                # Transition from Infectious to Recovered
                if np.random.rand() < gamma:
                    G.nodes[node]['state'] = 'R'

# Example simulation
G = nx.Graph()
venues = [
    Venue(venue_id=i, capacity=np.random.randint(50, 200), size=np.random.randint(100, 500), venue_type="community", transmission_characteristics={'max_distance': 10, 'max_duration': 30}),
    # Add more venues with their characteristics
]

n = 1000  # Number of individuals
for i in range(n):
    venue_choices = np.random.choice(venues, size=np.random.randint(1, 4))
    for venue in venue_choices:
        G.add_node(i, venue=venue, state='S')

initial_infected = np.random.choice(range(n), size=10, replace=False)  # Initial infected individuals
for infected_node in initial_infected:
    G.nodes[infected_node]['state'] = 'E'  # Assume initial exposed individuals

sigma = 0.1  # Exposed to Infectious rate
rho = 0.03  # Transmission probability
gamma = 0.1  # Infectious to Recovered rate

# Set up layers based on venue types
layers = list(set(venue.venue_type for venue in venues))
G_layers = {layer: G.copy() for layer in layers}

# Simulate for each layer
for layer, G_layer in G_layers.items():
    seir_simulation(G_layer, sigma, rho, gamma, simulation_days=10, distance_threshold=5, duration_threshold=20)

# Visualize the Network (optional)
for layer, G_layer in G_layers.items():
    nx.draw(G_layer, with_labels=True)
    plt.title(f"Layer: {layer}")
    plt.show()


##############


# import EoN
# import networkx as nx
# from collections import defaultdict
# import matplotlib.pyplot as plt
# import random

# N = 100000
# print('generating graph G with {} nodes'.format(N))
# G = nx.fast_gnp_random_graph(N, 5./(N-1))

# # Add random variation in the rate of leaving exposed class and partnership transmission rate
# # No variation in recovery rate

# node_attribute_dict = {node: 0.5+random.random() for node in G.nodes()}
# edge_attribute_dict = {edge: 0.5+random.random() for edge in G.edges()}

# nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
# nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')

# # These individual and partnership attributes will be used to scale transition rates
# # When we define \texttt{H} and \texttt{J}, we provide the name of these attributes

# # More advanced techniques to scale the transmission rates are shown in online documentation
# # Two directed graphs, the first for spontaneous transitions and the following for induced transitions

# # Spontaneous: A node of status A becomes B without any neighbor's influence
# # Not sure if we need something like this, since incubation/AKA E-->I is already defined in induced
# H = nx.DiGraph() 
# H.add_edge('E', 'I', rate = 0.0, weight_label='expose2infect_weight')
# H.add_node('S') # This line is actually unnecessary
# H.add_edge('I', 'R', rate = 0.0)

# # Induced: An edge between status A and status B nodes suddenly becomes an edge between a status A and a status C node
# # because the status B node changes status due to the existence of the edge
# # In principle, both could change status, but the code currently only allows the second node to change status
# J = nx.DiGraph()
# J.add_edge(('I', 'S'), ('I', 'E'), rate = 0.1, weight_label='transmission_weight') # Pull from roomProb
# IC = defaultdict(lambda: 'S') 
# for node in range(100): # PEOPLE INITIALLY IN INFECTED STATE
#     IC[node] = 'I'

# return_statuses = ('S', 'E', 'I', 'R')
                   
# print('Doing Gillespie simulation')
# t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = float('Inf'))
                                                                                          
# print('Done with simulation, now plotting')
# plt.plot(t, S, label = 'Susceptible')
# plt.plot(t, E, label = 'Exposed')
# plt.plot(t, I, label = 'Infected')
# plt.plot(t, R, label = 'Recovered')
# plt.xlabel('$t$')
# plt.ylabel('Simulated numbers')
# plt.legend()
# plt.show()

# # Visualize the graph
# plt.figure(figsize=(10, 10))
# pos = nx.spring_layout(G)  # You can use different layout algorithms
# nx.draw(G, pos, with_labels=False, node_size=10, node_color='blue', edge_color='gray', alpha=0.5)
# plt.title("Random Graph Visualization")
# plt.show()
