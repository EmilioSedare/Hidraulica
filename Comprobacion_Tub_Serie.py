import numpy as np

D=np.array([200,150,250])/1000
L=np.array([184,393,200])
km=np.array([7.1,11.2,11.2])
ks=np.array([0.0015,0.0015,0.0015])/1000
Qn=np.array([94,87,50])/1000
#Qf=87/1000
n_efi = 75/100
H=0
Z2=31.7




Ht=H-Z2
# Hallando Viscosidad y densidad:
T=15
vis = (1.14 - 0.031 * (T - 15) + 0.00068 * ((T - 15) ** 2)) * 10 ** (-6)  # m2/s
rho = 1000.1 + 0.0209 * T - 0.006 * T ** (2) + 2 * 10 ** (-5) * T ** 3  # kg/m3
g = 9.81  # m/s2
for i in range(len(L)):
    Qt=[]
    for n in range(len(Qn)):
        Qtt=np.sum(Qn[(0+n):(len(Qn))])
        Qt.append(Qtt)
    Qt = np.array((Qt))

    Vi = []
    for n in range(len(Qn)):
        V = (4*Qt[n])/(np.pi*D[n]**2)
        Vi.append(V)
    Vi=np.array((Vi))

    frr = []
    for n in range(len(L)):
        # Iteradores
        Re = (Vi[n] * D[n]) / vis
        error = 1
        Tol = 0.000001
        x_0 = 0.75
        if Re <= 2200:
            fr = 64 / Re
        else:
            while error >= Tol:
                f = -2 * np.log10(ks[n] / (3.7 * D[n]) + 2.51 * x_0 / Re) - x_0
                df = -2 / (np.log(10)) * ((2.51 / Re) / (ks[n] / (3.7 * D[n]) + 2.51 * x_0 / Re)) - 1
                x_i = x_0 - (f / df)
                error = abs(x_i - x_0)
                fr = 1 / (x_i ** 2)
                x_0 = x_i
            frr.append(fr)
    frr = np.array((frr))
    hf = frr * (L / D) * (Vi ** 2) / (2 * g)
    hm=km*Vi**2/(2*g)

htf=np.sum(hf)+np.sum(hm)
print(htf)
HT=abs(Ht)+htf
print(HT)


Pot = rho * np.sum(Qn) * g * HT / (n_efi * 1000)


print("Q= ", Qt)
print("V= ",Vi)
print("hf= ",hf)
print("hm= ",hm)
print("ht= ",HT)
print("po= ",Pot)
