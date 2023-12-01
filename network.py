import EoN
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import random

N = 100000
print('generating graph G with {} nodes'.format(N))
G = nx.fast_gnp_random_graph(N, 5./(N-1))

# Add random variation in the rate of leaving exposed class and partnership transmission rate
# No variation in recovery rate

node_attribute_dict = {node: 0.5+random.random() for node in G.nodes()}
edge_attribute_dict = {edge: 0.5+random.random() for edge in G.edges()}

nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')

# These individual and partnership attributes will be used to scale transition rates
# When we define \texttt{H} and \texttt{J}, we provide the name of these attributes

# More advanced techniques to scale the transmission rates are shown in online documentation
# Two directed graphs, the first for spontaneous transitions and the following for induced transitions

# Spontaneous: A node of status A becomes B without any neighbor's influence
# Not sure if we need something like this, since incubation/AKA E-->I is already defined in induced
H = nx.DiGraph() 
H.add_edge('E', 'I', rate = 0.6, weight_label='expose2infect_weight')
H.add_node('S') # This line is actually unnecessary
H.add_edge('I', 'R', rate = 0.1)

# Induced: An edge between status A and status B nodes suddenly becomes an edge between a status A and a status C node
# because the status B node changes status due to the existence of the edge
# In principle, both could change status, but the code currently only allows the second node to change status
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate = 0.1, weight_label='transmission_weight') # Pull from roomProb
IC = defaultdict(lambda: 'S') 
for node in range(100): # PEOPLE INITIALLY IN INFECTED STATE
    IC[node] = 'I'

return_statuses = ('S', 'E', 'I', 'R')
                   
print('Doing Gillespie simulation')
t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = float('Inf'))
                                                                                          
print('Done with simulation, now plotting')
plt.plot(t, S, label = 'Susceptible')
plt.plot(t, E, label = 'Exposed')
plt.plot(t, I, label = 'Infected')
plt.plot(t, R, label = 'Recovered')
plt.xlabel('$t$')
plt.ylabel('Simulated numbers')
plt.legend()
plt.show()

# Visualize the graph
plt.figure(figsize=(10, 10))
pos = nx.spring_layout(G)  # You can use different layout algorithms
nx.draw(G, pos, with_labels=False, node_size=10, node_color='blue', edge_color='gray', alpha=0.5)
plt.title("Random Graph Visualization")
plt.show()
