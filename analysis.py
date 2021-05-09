import matplotlib.pyplot as plt
import numpy as np
import math


def LeastSquareMethod(xlst, ylst):
    n = len(xlst)
    avgx = sum(xlst)/n
    avgy = sum(ylst)/n
    sumxy = sum([xlst[i]*ylst[i] for i in range(n)])
    sumx2 = sum([xi*xi for xi in xlst])
    sumy2 = sum([yi*yi for yi in ylst])
    
    k = (sumxy/n - avgx * avgy) / (sumx2/n - avgx * avgx)
    b = avgy - k * avgx
    r = (sumxy - n*avgx*avgy)/math.sqrt((sumx2-n*avgx*avgx)*(sumy2-n*avgy*avgy))
    return k,b,r


felst = []
nitrolst = []
fepeaklst = []
nitropeaklst = []

fp = open("Fe.asc","r")
for i in range(512):
    s = fp.readline().split()
    felst.append(int(s[1]))
fp.close()

fp = open("nitro.asc","r")
for i in range(512):
    s = fp.readline().split()
    nitrolst.append(int(s[1]))
fp.close()

fepeaklst = [37, 76, 115, 146, 185, 225, 287, 326, 366, 395, 435, 474]
nitropeaklst = [120, 151, 361, 391]
plt.plot([i for i in range(512)], felst)
plt.scatter(fepeaklst, [felst[xi] for xi in fepeaklst], marker="x", color="red")
plt.ylim(0,10000)
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("fespec.png")
plt.savefig("fespec.eps")
#plt.show()
plt.clf()

plt.plot([i for i in range(512)], nitrolst)
plt.scatter(nitropeaklst, [nitrolst[xi] for xi in nitropeaklst], marker="x", color="red")
plt.ylim(0,10000)
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("nitrospec.png")
plt.savefig("nitrospec.eps")
#plt.show()
plt.clf()

plt.plot([i for i in range(512)], felst)
plt.plot([i for i in range(512)], nitrolst)
plt.scatter(fepeaklst, [felst[xi] for xi in fepeaklst], marker="x", color="red")
plt.scatter(nitropeaklst, [nitrolst[xi] for xi in nitropeaklst], marker="x", color="red")
plt.ylim(0,10000)
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("jointspec.png")
plt.savefig("jointspec.eps")
#plt.show()
plt.clf()

k = [fepeaklst[5]-fepeaklst[0], fepeaklst[6]-fepeaklst[11]]
k = [10.656/ki for ki in k]
print("K VALUE (MM/S)")
print(k[0],k[1])
print("")

c = [fepeaklst[0]+fepeaklst[1]+fepeaklst[4]+fepeaklst[5],fepeaklst[6]+fepeaklst[7]+fepeaklst[10]+fepeaklst[11]]
c = [ci/4.0 for ci in c]
print("CENTROID (ADDRESS)")
print(c[0],c[1])
print("")

zsshift = [-0.185/ki for ki in k]
print("ZERO SPEED SHIFT (ADDRESS)")
print(zsshift[0],zsshift[1])
print("")

zs = [c[i]+zsshift[i] for i in range(2)]
print("ZERO SPEED (ADDRESS)")
print(zs[0],zs[1])
print("")

fezeeman = [fepeaklst[:6],fepeaklst[6:]]
indlst = [i+1 for i in range(6)]
lsmlst = [LeastSquareMethod(indlst,fezeemani) for fezeemani in fezeeman]
print("LEAST SQUARED ANALYSIS")
print("\tINTERVAL (ADDRESS)")
print(lsmlst[0][0],lsmlst[1][0])
print("\tR VALUE")
print(lsmlst[0][2],lsmlst[1][2])
print("\tINTERVAL (MM/S)")
print(lsmlst[0][0]*k[0],-lsmlst[1][0]*k[1])
print("")

avgde = (lsmlst[0][0]*k[0]-lsmlst[1][0]*k[1])/2.0*4.80766e-8
print("AVERAGE DELTA E (EV)")
print(avgde)

deg = [fepeaklst[3]+fepeaklst[4]-fepeaklst[1]-fepeaklst[2],fepeaklst[7]+fepeaklst[8]-fepeaklst[9]-fepeaklst[10]]
deg = [deg[i]/2.0*k[i] for i in range(2)]
avgdeg = np.average(deg)
print("GROUND STATE ENERGY DIFFERENCE IN MM/S")
print(deg,avgdeg)
dee = [-fepeaklst[3]+fepeaklst[5]-fepeaklst[0]+fepeaklst[2],fepeaklst[6]-fepeaklst[8]+fepeaklst[9]-fepeaklst[11]]
dee = [dee[i]/4.0*k[i] for i in range(2)]
avgdee = np.average(dee)
print("EXCITED STATE ENERGY DIFFERENCE IN MM/S")
print(dee,avgdee)
print("GROUND STATE ENERGY DIFFERENCE IN EV")
print(avgdeg*4.80766e-8)
print("EXCITED STATE ENERGY DIFFERENCE IN EV")
print(avgdee*4.80766e-8)

cn = [nitropeaklst[1]+nitropeaklst[0],nitropeaklst[3]+nitropeaklst[2]]
cn = [ci/2.0 for ci in cn]
dn = [cn[i] - c[i] for i in range(2)]
ddn = [dn[i]*k[i]*4.80766e-8 for i in range(2)]
print(ddn,np.average(ddn))

td = [nitropeaklst[1] - nitropeaklst[0], nitropeaklst[2] - nitropeaklst[3]]
dtd = [td[i]*k[i]*4.80766e-8 for i in range(2)]
print(dtd,np.average(dtd))

noise = np.average(nitrolst[200:300])
plt.plot([i for i in range(512)], [xi - noise for xi in nitrolst])
for i in range(4):
    plt.plot([j for j in range(512)],[(nitrolst[nitropeaklst[i]]-noise)/2 for j in range(512)],linestyle="--",label="%d"%i)
plt.legend()
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("nitro-spec.png")
plt.savefig("nitro-spec..eps")
#plt.show()
plt.clf()

xxlst = []
for i in range(4):
    nitrolst = [ni-noise for ni in nitrolst]
    target = nitrolst[nitropeaklst[i]]/2
    for j in range(511):
#        print(nitrolst[j+1],target,nitrolst[j])
        if (((nitrolst[j+1]<target) and (nitrolst[j]>target)) or ((nitrolst[j+1]>target) and (nitrolst[j]<target))):
            _k = (nitrolst[j+1]-nitrolst[j])
            xxlst.append((nitrolst[j+1]-nitrolst[j])/_k+j) 
print(xxlst)

hwlst = [xxlst[2*j+1]-xxlst[2*j] for j in range(4)]
print(hwlst,np.average(hwlst))

noise = np.average(felst[240:270])
plt.plot([i for i in range(512)], [xi - noise for xi in felst])
for i in range(12):
    plt.plot([j for j in range(512)],[(felst[fepeaklst[i]]-noise)/2 for j in range(512)],linestyle="--",label="%d"%i)
plt.legend()
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("fe-spec.png")
plt.savefig("fe-spec.eps")
#plt.show()
plt.clf()
xxlst = []
for i in range(12):
    felst = [ni-noise for ni in felst]
    target = felst[fepeaklst[i]]/2
    for j in range(511):
#        print(felst[j+1],target,felst[j])
        if (((felst[j+1]<target) and (felst[j]>target)) or ((felst[j+1]>target) and (felst[j]<target))):
            _k = (felst[j+1]-felst[j])
            xxlst.append((felst[j+1]-felst[j])/_k+j) 
print(xxlst)
hwlst = [xxlst[2*j+1]-xxlst[2*j] for j in range(11)]
hwlst.remove(1.0)
hwlst.remove(1.0)
hwlst.remove(1.0)
print(hwlst,np.average(hwlst))