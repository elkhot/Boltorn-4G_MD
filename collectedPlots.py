import newB30M
import Bolope10M
import Bolope1M
import Bolope100K
import Bolope10K
import Bolope1K
import matplotlib.pyplot as plt

plt.style.use("ggplot")
plt.rcParams['axes.edgecolor'] = "#777777"
plt.rcParams['axes.facecolor'] = '#FFFFFF'


time_30M = []
for x in newB30M.time_30M:
    time_30M.append(x)

time_10M  = []
for x in Bolope10M.time_10M:
    time_10M.append(x)

time_1M = []
for x in Bolope1M.time_1M:
    time_1M.append(x)

time_100K = []
for x in Bolope100K.time_100K:
    time_100K.append(x)

time_10K = []
for x in Bolope10K.time_10K:
   time_10K.append(x)

time_1K = []
for x in Bolope1K.time_1K:
    time_1K.append(x)




averageDistances_30M = []
averageDistances = []
radiiOfGyration_30M =[]
radiiOfGyration = []
for x in newB30M.averageDistances_30M:
    averageDistances_30M.append(x)
    averageDistances.append(x)
for x in newB30M.radiiOfGyration_30M:
    radiiOfGyration_30M.append(x)
    radiiOfGyration.append(x)

averageDistances_10M  = []
radiiOfGyration_10M = []
for x in Bolope10M.averageDistances_10M:
    averageDistances_10M.append(x)
    averageDistances.append(x)
for x in Bolope10M.radiiOfGyration_10M:
    radiiOfGyration_10M.append(x)
    radiiOfGyration.append(x)

averageDistances_1M = []
radiiOfGyration_1M = []
for x in Bolope1M.averageDistances_1M:
    averageDistances_1M.append(x)
    averageDistances.append(x)
for x in Bolope1M.radiiOfGyration_1M:
    radiiOfGyration_1M.append(x)
    radiiOfGyration.append(x)

averageDistances_100K = []
radiiOfGyration_100K = []
for x in Bolope100K.averageDistances_100K:
    averageDistances_100K.append(x)
    averageDistances.append(x)
for x in Bolope100K.radiiOfGyration_100K:
    radiiOfGyration_100K.append(x)
    radiiOfGyration.append(x)
    
averageDistances_10K = []
radiiOfGyration_10K = []
for x in Bolope10K.averageDistances_10K:
    averageDistances_10K.append(x)
    averageDistances.append(x)
for x in Bolope10K.radiiOfGyration_10K:
    radiiOfGyration_10K.append(x)
    radiiOfGyration.append(x)
    
averageDistances_1K = []
radiiOfGyration_1K = []
for x in Bolope1K.averageDistances_1K:
    averageDistances_1K.append(x)
    averageDistances.append(x)
for x in Bolope1K.radiiOfGyration_1K:
    radiiOfGyration_1K.append(x)
    radiiOfGyration.append(x)
    

pairs_30M = [[t,r] for t, r in zip(time_30M, averageDistances_30M)]

pairsRg_30M = [[t,r] for t, r in zip(time_30M, radiiOfGyration_30M)]

pairs_10M = [[t,r] for t, r in zip(time_10M, averageDistances_10M)]

pairsRg_10M = [[t,r] for t, r in zip(time_10M, radiiOfGyration_10M)]

pairs_1M = [[t,r] for t, r in zip(time_1M, averageDistances_1M)]

pairsRg_1M = [[t,r] for t, r in zip(time_1M, radiiOfGyration_1M)]

pairs_100K = [[t,r] for t, r in zip(time_100K, averageDistances_100K)]

pairsRg_100K = [[t,r] for t, r in zip(time_100K, radiiOfGyration_100K)]

pairs_10K = [[t,r] for t, r in zip(time_10K, averageDistances_10K)]

pairsRg_10K = [[t,r] for t, r in zip(time_10K, radiiOfGyration_10K)]

pairs_1K = [[t,r] for t, r in zip(time_1K, averageDistances_1K)]

pairsRg_1K = [[t,r] for t, r in zip(time_1K, radiiOfGyration_1K)]

plt.figure(figsize=(27,27))
plt.subplot(321)
plt.scatter(*zip(*pairs_1K), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 1K steps', pad=9, weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)

plt.subplot(322)
plt.scatter(*zip(*pairs_10K), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 10K steps', pad=9, weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)

plt.subplot(323)
plt.scatter(*zip(*pairs_100K), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 100K steps', pad=9, weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)

plt.subplot(324)
plt.scatter(*zip(*pairs_1M), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 1M steps', pad=9, weight="bold", fontsize=15)
plt.xlabel('Time (s)',labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)

plt.subplot(325)
plt.scatter(*zip(*pairs_10M), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 10M steps',pad=6, weight="bold", fontsize=15)
plt.xlabel('Time (s)',labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)

plt.subplot(326)
plt.scatter(*zip(*pairs_30M), color='b', s=7, marker="o")
plt.title('Atom-center of mass distance at 30M steps',pad=6, weight="bold", fontsize=15)
plt.xlabel('Time (s)',labelpad=2)
plt.ylabel('distance(nm)', labelpad=2)


plt.tight_layout(pad=20)

plt.show()



plt.figure(figsize=(27,27))
plt.subplot(321)
plt.scatter(*zip(*pairsRg_1K), color='k', s=7, marker="o")
plt.title('Rg at 10K steps',pad=6, weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('Rg(nm)', labelpad=2)

plt.subplot(322)
plt.scatter(*zip(*pairsRg_10K), color='k', s=7, marker="o")
plt.title('Rg at 10K steps',pad=6 ,weight="bold", fontsize=15)
plt.xlabel('Time (s)',labelpad=2)
plt.ylabel('Rg(nm)', labelpad=2)

plt.subplot(323)
plt.scatter(*zip(*pairsRg_100K), color='k', s=7, marker="o")
plt.title('Rg at 100K steps',pad=6,weight="bold" , fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('Rg(nm)', labelpad=2)

plt.subplot(324)
plt.scatter(*zip(*pairsRg_1M), color='k', s=7, marker="o")
plt.title('Rg at 1M steps',pad=6,weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('Rg(nm)',labelpad=2)

plt.subplot(325)
plt.scatter(*zip(*pairsRg_10M), color='k', s=7, marker="o")
plt.title('Rg at 10M steps',pad=6,weight="bold", fontsize=15)
plt.xlabel('Time (s)',labelpad=2)
plt.ylabel('Rg(nm)',labelpad=2)

plt.subplot(326)
plt.scatter(*zip(*pairsRg_30M), color='k', s=7, marker="o")
plt.title('Rg at 30M steps',pad=6,weight="bold", fontsize=15)
plt.xlabel('Time (s)', labelpad=2)
plt.ylabel('Rg(nm)',labelpad=2)
plt.tight_layout(pad=20)

plt.show()



time = []
for x in newB30M.time_30M:
    time.append(x)

for x in Bolope10M.time_10M:
    time.append(x)


for x in Bolope1M.time_1M:
    time.append(x)

for x in Bolope100K.time_100K:
    time.append(x)

for x in Bolope10K.time_10K:
    time.append(x)


for x in Bolope1K.time_1K:
    time.append(x)


pairs = [[t,r] for t, r in zip(time, averageDistances)]

pairsRg = [[t,r] for t, r in zip(time, radiiOfGyration)]

plt.figure(figsize=(20,20))
plt.scatter(*zip(*pairs) ,color='k', s=7, marker="o")
plt.title('The average atom-center of mass distance vs. log time', pad='4', weight='bold', fontsize=30)
plt.xscale('log')
plt.xlabel('log time',fontsize=20)
plt.ylabel('distance (nm)', fontsize=20)
plt.tight_layout(pad=10)
plt.show()
 



plt.figure(figsize=(20,20))
plt.scatter(*zip(*pairsRg) ,color='k', s=7, marker="o")
plt.title('The radius of gyration vs log time',pad='4', weight='bold',fontsize=30)
plt.xscale('log')
plt.xlabel('log time',fontsize=15)
plt.ylabel('radius of gyration (nm)', fontsize=20)
plt.tight_layout(pad=10)
plt.show()
 