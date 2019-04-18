import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import os

# ------------------------------------------------------------------------------- HYPERPARAMETERS SIS
# mu = 1  # spontaneous recovery probability (infected -> susceptible)
# posible_betas = np.linspace(start=0, stop=1, num=101)  # infection probability of a susceptible individual when
# contacted by an infected
# ro_init = 0.2  # initial fraction of infected nodes
# Tmax = 1000  # number of timesteps of the simulation
# Ttrans = Tmax-100  # transitory steps
# Nrep = 50  # Number of repetitions

mu = 1  # spontaneous recovery probability (infected -> susceptible)
posible_betas = [0.18, 0.19, 0.90]  # infection probability of a susceptible individual when contacted by an infected
ro_init = 0.2  # initial fraction of infected nodes
Tmax = 50  # number of timesteps of the simulation
Ttrans = Tmax-10  # transitory steps
Nrep = 1  # Number of repetitions

path = './Networks'
for r, d, f in os.walk(path):
    for file in f:
        # --------------------------------------------------------------------------------------------- GENERATE NETWORK
        # Scale free graph
        print('Computing file '+file)
        graph = file
        G = nx.Graph(nx.read_pajek('Networks/' + graph))

        # Random graph
        # n = 100
        # p = 0.1
        # G = nx.gnp_random_graph(n, p, directed=False)

        # ------------------------------------------------------------------------------------------------- DRAW NETWORK
        fig = plt.figure()
        nx.draw_networkx(G, with_labels=False, node_size=20)
        limits = plt.axis('off')
        plt.suptitle('Network: ' + graph + '\n\n ', fontsize=16)
        fig.savefig('./Images/' + graph + '.png')

        # ------------------------------------------------------------------------------------------------ SIS ALGORITHM
        Nnodes = len(G.nodes)

        ro_over_beta = []
        ro_over_time = np.zeros((6,Tmax))
        itt_beta = 0
        for beta in posible_betas:
            p_Nrep_time = np.zeros((Nrep,Tmax))
            for tries in range(Nrep):
                init_array = list(G.nodes)
                random.shuffle(init_array)
                infected = init_array[0:round(ro_init * Nnodes)]
                suscepti = [node for node in init_array if node not in infected]

                for tstep in range(Tmax):
                    p_Nrep_time[tries,tstep] = len(infected)/Nnodes
                    recovered = []
                    nowinfect = []

                    # 1. For each infected node at time step t, we recover it with probability mu:
                    for ninf in infected:
                        if random.uniform(0, 1) < mu:
                            recovered.append(ninf)

                    # 2. For each susceptible node at time step t, we traverse all of its neighbors
                    # For each infected neighbor (at time step t), the reference node becomes infected with probability beta
                    for susnode in suscepti:
                        neighbours = list(G.adj[susnode])
                        infected_nei = [infnei for infnei in neighbours if infnei in infected]
                        counter = 0
                        enter = True
                        while (enter is True) and (counter < len(infected_nei)):
                            if random.uniform(0, 1) < beta:
                                nowinfect.append(susnode)
                                enter = False
                            counter = counter + 1

                    # 3. Actualize infected() and suscepti()
                    for item in recovered:
                        infected.remove(item)
                    infected = infected + nowinfect

                    suscepti = [node for node in init_array if node not in infected]

            if beta in [0.18, 0.19, 0.20, 0.25, 0.40, 0.90]:
                ro_over_time[itt_beta] = p_Nrep_time[0,:]
                itt_beta = itt_beta + 1
            ro_over_beta.append(np.mean(p_Nrep_time[:,Ttrans:]))

        # -------------------------------------------------------------------------------------------------------- PLOTS
        # Plot one repetition of ro over time for multiple betas
        fig = plt.figure(figsize=(10, 8))
        for iteration in range(ro_over_time.shape[0]):
            plt.plot(ro_over_time[iteration,:])
        plt.ylabel('Fraction of infected nodes', fontsize=14)
        plt.xlabel('Time t', fontsize=14)
        plt.ylim((0,1))
        plt.grid(color='gray', linestyle='-', linewidth=1)
        plt.xlim((0,Tmax))
        plt.suptitle('Network: '+graph[:-4] +' (mu='+str(mu)+', p(0)='+str(ro_init)+')',fontsize=17)
        fig.savefig('./Images/'+'['+graph+']'+'m'+str(mu)+'r'+str(ro_init)+'ro_over_time'+'.png')

        # Plot ro over beta
        fig = plt.figure(figsize=(10, 8))
        plt.plot(posible_betas, ro_over_beta, '*')
        plt.ylabel('Fraction of infected nodes', fontsize=14)
        plt.xlabel('Beta', fontsize=14)
        plt.ylim((0, 1))
        plt.grid(color='gray', linestyle='-', linewidth=1)
        plt.xlim((0, Tmax))
        plt.suptitle('Network: ' + graph[:-4] +' (mu='+str(mu)+', p(0)='+str(ro_init)+')', fontsize=17)
        fig.savefig('./Images/'+'[' + graph + ']'+'m'+str(mu)+'r'+str(ro_init)+'ro_over_beta'+'.png')

