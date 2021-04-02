"""copied code since previous script has interpreter issues"""
'''
Created on Feb 27, 2018
Author: <redacted for anonymised review>
'''

import random as rnd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('seaborn-talk')
default_path = os.getcwd()


class Individual:
    '''Class that defines an individual as an object'''

    def __init__(self,
                 x,
                 y):
        # Initialization
        self.x = x  # location (x, y)
        self.y = y

    def move(self, connection, max_x, max_y):
        '''the individual moves'''
        # an individual moves
        i = rnd.choice(*np.where(connection[self.x * max_y + self.y]))
        self.x, self.y = i // max_y, i % max_y


class Metapopulation:
    '''Contains the whole population, regulates daily affairs'''

    def __init__(self,
                 max_x,
                 max_y,
                 lambd,
                 K,
                 d,
                 m):
        """Initialization"""
        self.max_x = max_x  # number of grid cells along the first dimension of the landscape
        self.max_y = max_y  # number of grid cells along the second dimension of the landscape
        self.lambd = lambd  # Local optimal growth rate/optimal fecundity
        self.K = K  # Local Carrying capacity
        self.d = d  # dispersal propensity
        self.m = m  # dispersal mortality
        self.population = []  #population object: container for all individuals
        self.localsize = np.zeros((self.max_x, self.max_y)) #keeps track of locl population sizes in all patches
        self.immi = np.zeros((self.max_x, self.max_y)) #keeps track of the number of immigrant into every patch
        self.initialize_pop()
        self.connection_matrix()

    def initialize_pop(self):
        '''Initialize starting population'''
        startpop = 0.5 * self.max_x * self.max_y * K  # initial metapopulation size

        for _ in range(int(startpop)):
            x, y = rnd.randint(0, (self.max_x - 1)), rnd.randint(0, (self.max_y - 1))
            self.population.append(Individual(x, y))
            self.localsize[x, y] += 1

    def connection_matrix(self):
        '''defines the connections in the network of patches'''
        self.connection = []
        for x in range(self.max_x):
            for y in range(self.max_y):
                conn = [self.max_y * g + h for g in (x - 1, x, x + 1) for h in (y - 1, y, y + 1) if
                        (g, h) != (x, y) and 0 <= g < self.max_x and 0 <= h < self.max_y]

                self.connection.append([(i in conn) for i in range(self.max_x * self.max_y)])
        self.connection = np.array(self.connection)

    def lifecycle(self):
        '''all actions for all individuals during one time step (1 generation)'''

        self.oldpop = self.population[:] #make a copy of the population
        self.oldlocalsize = self.localsize[:]
        self.immi = np.zeros((self.max_x, self.max_y))  # to record immigrants

        rnd.shuffle(self.oldpop)

        self.new = []
        self.newsize = np.zeros((self.max_x, self.max_y))  # to record local population sizes
        # new generation of population
        self.population = []  #empty this copy of the population to, then, add indivduals of the next generation
        self.localsize = np.zeros((self.max_x, self.max_y))

        for ind in self.oldpop:
            # reproduce
            exp_off = r * (1 + (self.oldlocalsize[ind.x, ind.y] * (r - 1)) / a) ** -b  # Hassell
            for _ in range(np.random.poisson(max(0, exp_off))):
                new = Individual(ind.x, ind.y)
                # offspring moves
                if self.d > rnd.random():
                    new.move(self.connection, self.max_x, self.max_y)
                    # dispersal mortality
                    if rnd.random() > self.m:
                        self.population.append(new)
                        self.localsize[new.x, new.y] += 1
                        self.immi[new.x, new.y] += 1
                else:
                    self.population.append(new)
                    self.localsize[new.x, new.y] += 1


class Datacollector:
    '''extracts data during the metapopulation runs, calculates derived metrics and plots'''
    def __init__(self, maxtime, max_x, max_y):
        self.mt = maxtime
        self.max_x = max_x
        self.max_y = max_y
        self.sizesintime = np.zeros((maxtime, max_x, max_y))

    def collect(self, time, localsize):
        #collect the array that stores local population size for every time step of a run
        self.sizesintime[time] = localsize

    def plotlocal(self):
        '''plots local population sizes in time for every patch'''
        fig, axes = plt.subplots(self.max_x, self.max_y)
        for i, ax in enumerate(axes.flatten()):
            ax.plot(self.sizesintime[:, i // self.max_y, i % self.max_y])
            ax.set_xlabel('time')
            ax.set_title(i)

        fig.suptitle('local population sizes')
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.show()

    def variabilities(self, extent):
        '''returns alpha, beta and gamma variabilities from the localpopulation sizes per time step for a given window (extent)
        the extent enables to exclude dynamics that occur before stability in the system'''
        reshaped_sizes = np.reshape(self.sizesintime, (self.mt, self.max_x * self.max_y))[int(self.mt - extent):]
        total = np.sum(reshaped_sizes, axis=1)
        total_mean = np.mean(total)
        total_var = np.var(total) / total_mean ** 2
        covars = np.cov(np.transpose(reshaped_sizes))
        alpha = (np.sum(np.diag(covars) ** 0.5) / total_mean) ** 2
        gamma = np.sum(covars) / (total_mean ** 2)
        beta = alpha - gamma
        beta_ = alpha / gamma
        s_corner = np.mean(reshaped_sizes[:,(0,2,6,8)], axis = 0)
        s_side = np.mean(reshaped_sizes[:,(1,3,5, 7)],  axis=0)
        s_center = np.mean(reshaped_sizes[:, 4] , axis=0)
        alpha_corner = np.sum(np.std(reshaped_sizes[:, (0, 2, 6, 8)], axis=0) / s_corner)
        alpha_side = np.sum(np.std(reshaped_sizes[:, (1, 3, 5, 7)], axis=0) / s_side)
        alpha_center = np.std(reshaped_sizes[:, 4]) / s_center

        return alpha, gamma, beta, beta_, alpha_corner, alpha_side, alpha_center, np.mean(s_corner), np.mean(s_side), np.mean(s_center)


def run(r, K, d, m):
    '''one simulated landscape/metapopulation that returns the datacollector object'''
    meta = Metapopulation(max_x, max_y, r, K, d, m)
    data = Datacollector(maxtime, max_x, max_y)
    # simulate MAXTIME timesteps (print generation time and metapopulation size for quickly checking during runs)
    for timer in range(maxtime):
        meta.lifecycle()
        data.collect(timer, meta.localsize)
        """print(meta.localsize)
        unique = np.unique([ind.y+max_y*ind.x for ind in meta.population], return_counts = 1)
        print(np.reshape(unique[1], (max_x, max_y)))"""
        # print('generation ',timer)
        # print("popsize: {}\n".format(len(meta.population)))
        # print(f'{meta.localsize[0, 0]}, {meta.localsize[0, 1]}')
    # print(str(time.clock()))
    # data.plotlocal()
    # data.plottotal()
    print(meta.immi)
    return data

def varplot(data, data_mean, x, y):
    '''plot a variation metric conditional o treatment in the virtual experiment'''
    fig, ax = plt.subplots()
    plot1 = ax.scatter(data[x], data[y])
    plot1_1 = ax.plot(data_mean[x], data_mean[y])
    plt.savefig(f'Hassell/{x}/{y}.png')

def partplot(data, data_mean, x, y):
    '''plot local metric conditional on the type of position in the landscape (center, side or corner)'''
    disp = {'alphavar': r'$\alpha$-variation', 'size': 'population size'}
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    ax1.set_title(f'{disp[y]} corner')
    ax2.set_title(f'{disp[y]} side')
    ax3.set_title(f'{disp[y]} center')

    ax1.scatter(data[x], data[f'{y} corner'])
    ax1.plot(data_mean[x], data_mean[f'{y} corner'])
    ax2.scatter(data[x], data[f'{y} side'])
    ax2.plot(data_mean[x], data_mean[f'{y} side'])
    ax3.scatter(data[x], data[f'{y} center'])
    ax3.plot(data_mean[x], data_mean[f'{y} center'])
    plt.savefig(f'Hassell/{x}/{y}part.png')

def replicate_runs(n, v):
    '''define replicated virtual experiment for one of the parameters (dispersal propensity or dispersal mortality)'''
    x_s = []
    y_s = []
    a_s = []
    g_s = []
    b_s = []
    b__s = []
    a_cos = []
    a_sis = []
    a_ces = []
    s_cos = []
    s_sis = []
    s_ces = []
    #virtual treatment levels
    for var in np.arange(0.05, 0.85, 0.05) if v == 'disp_mort' else np.arange(0.05, 1, 0.05):
        #replication
        for i in range(n):
            x_s += [var]
            #start one simulation
            data = run(r, K, var if v == 'disp_prop' else d, var if v == 'disp_mort' else m)
            alph, gam, bet, bet_, a_corner, a_side, a_center, s_corner, s_side, s_center = data.variabilities(maxtime//2)
            a_s += [alph]
            g_s += [gam]
            b_s += [bet]
            b__s += [bet_]
            a_cos += [a_corner]
            a_sis += [a_side]
            a_ces += [a_center]
            s_cos += [s_corner]
            s_sis += [s_side]
            s_ces += [s_center]

            print(f'{var}, replicate {i}')
    #construct an easy-to-analyze dataframe
    data = pd.DataFrame({v: x_s,
                         'alphavar': a_s,
                         'gammavar': g_s,
                         'betavar': b_s,
                         'betavar_': b__s,
                         'alphavar corner' : a_cos,
                         'alphavar side' : a_sis,
                         'alphavar center' : a_ces,
                         'size corner' : s_cos,
                         'size side' : s_sis,
                         'size center' : s_ces})
    data_mean = data.groupby(v).mean().reset_index()
    print(data_mean)
    #save an output file
    data.to_csv('var_sims')
    #plots
    varplot(data, data_mean, v, 'alphavar')
    varplot(data, data_mean, v, 'gammavar')
    varplot(data, data_mean, v, 'betavar_')
    partplot(data, data_mean, v, 'size')
    partplot(data, data_mean, v, 'alphavar')





##############
# Parameters #
##############
maxtime = 500
r, K = 2, 100
b = 2
a = K * (r - 1) / (r ** (1 / b) - 1)
d = 0.25
m = 0.25
max_x, max_y = 3, 3

#replicated model runs varied in a gradient of either dispersal  ('disp_prop') or dispersal mortality ('disp_mort')
replicate_runs(20, 'disp_mort')
