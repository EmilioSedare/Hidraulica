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





''' 
err=Ht
Hti,hf0=0,0
Toll=0.000001
while err >= Toll:
    #tuberia inicial
    hf0=(Ht-Hti)*(L[0]/D[0]**5)/(np.sum(L/D**5))+hf0
    Vi=-2*((2*g*D[0]*hf0)**0.5/(L[0])**0.5)*np.log10(ks[0] / (3.7 * D[0]) + (2.51 * vis*L[0]**0.5) /(D[0]*(2*g*D[0]*hf0)**0.5))
    Qi=np.pi/4*D[0]**2*Vi
    hm0=km[0]*Vi**2/(2*g)
    Di=D[1:,]
    Li=L[1:,]
    kmi=km[1:,]
    ksi=ks[1:,]
    for i in range(len(Li)):
        #i+1 TUB
        Qii=np.zeros(len(Qn))
        for n in range(len(Qn)):
            Qii[n] = Qi - np.sum(Qn[0:(n+1)])
        Vii = []
        for n in range(len(Qn)):
            V=(4/np.pi)*Qii[n]/(Di[n]**2)
            Vii.append(V)

        frr = []
        for n in range(len(Li)):
            # Iteradores
            Re = (Vii[n] * Di[n]) / vis
            error = 1
            Tol = 0.000001
            x_0 = 0.75
            if Re <= 2200:
                fr = 64 / Re
            else:
                while error >= Tol:
                    f = -2 * np.log10(ksi[n] / (3.7 * Di[n]) + 2.51 * x_0 / Re) - x_0
                    df = -2 / (np.log(10)) * ((2.51 / Re) / (ksi[n] / (3.7 * Di[n]) + 2.51 * x_0 / Re)) - 1
                    x_i = x_0 - (f / df)
                    error = abs(x_i - x_0)
                    fr = 1 / (x_i ** 2)
                    x_0 = x_i
            frr.append(fr)
        hfi=[]
        hmi=[]
        for n in range(len(Li)):
            hf = frr[n] * (Li[n] / Di[n]) * (Vii[n] ** 2) / (2 * g)
            hm = kmi[n] * (Vii[n] ** 2) / (2 * g)
            hfi.append(hf)
            hmi.append(hm)

    Hti=hf0+hm0+np.sum(hfi)+np.sum(hmi)
    #Hti=Ht+Hti
    #Pot = rho * Q * g * Hbomba / (n * 1000)
    dh=abs(Ht-Hti)
    err=dh

Q=np.append(Qi,Qii)
V=np.append(Vi,Vii)
hf=np.append(hf0,hfi)
hm=np.append(hm0,hmi)

print("Q= ", Q)
print("V= ",V)
print("hf= ",hf)
print("hm= ",hm)
'''