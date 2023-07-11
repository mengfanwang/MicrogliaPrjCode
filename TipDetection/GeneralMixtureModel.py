from scipy.io import loadmat
import numpy as np
from pomegranate import *
from math import exp

# for step in range(6, 22):
#     print(step)
#     distance_map = loadmat('D:/dropbox/distance_map/' + str(step) + '/1')
#     distance_map = distance_map['distance_map']
#     distance_map[distance_map < 0.01] = 0
#     distance_map[distance_map > step - 0.01] = step - 0.01
#     distance_map = distance_map / step
#     if step == 6:
#         distance_map_max = distance_map
#     else:
#         distance_map_max = np.maximum(distance_map, distance_map_max)
# x = distance_map_max[(distance_map_max < 1) & (distance_map_max > 0)]
# x = x.reshape(-1, 1)
x = loadmat('./distance_map_max.mat')
x = x['x']
x[x>0.89] = 0.89
x = x/np.max(x)
model = GeneralMixtureModel.from_samples([BetaDistribution, BetaDistribution], n_components=2, X=x)
print(model.distributions)
print([exp(model.weights[0]), exp(model.weights[1])])

# model_val1 = GammaDistribution(model.distributions[0].parameters[0], model.distributions[0].parameters[1])
# model_val2 = NormalDistribution(model.distributions[1].parameters[0], model.distributions[1].parameters[1])
# model_val = GeneralMixtureModel([model_val1, model_val2], weights=[exp(model.weights[0]), exp(model.weights[1])])
# x_val = model_val.sample(len(x))



import matplotlib.pyplot as plt
binEdge = np.linspace(0, 1, 100)
plt.hist(x, bins=binEdge, alpha=0.5, label='data')
# plt.hist(x_val, bins=binEdge, alpha=0.5, label='model')
plt.legend()
plt.show()