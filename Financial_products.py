import numpy as np
import matplotlib.pyplot as plt
from math import *
import random
from mpl_toolkits import mplot3d
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from tkinter import *

#-------------- Macros --------------

M1 = 99
N1 = 49
N = 99
M = 4999
T = 0.5
K =  10
K1 = 3
r = 0.1
sigma = 0.5
L = 20
delta_s = L/(N+1)
delta_t = T/(M+1)
S = np.linspace(0,L,N+2)
t = np.linspace(0,T,M+2)
S1 = np.linspace(0,L,N1+2)
t1 = np.linspace(0,T,M1+2)
gamma = 0.2316419
a1 = 0.319381530
a2 = -0.356563782
a3 = 1.781477937
a4 = -1.821255978
a5 = 1.330274429
Nme = 100
S0=10


#------------------------------------


#------------------------- Methode de differences finis -------------------------------------

def Max(a,b):
    if a > b:
        return a
    else :
        return b


def PayOff_Call_Euro(k,s):
    return Max(0, s-k)




def BS_Dirichlet_Call_Euro():
    V = np.zeros((M+2,N+2))
    for i in range(N+2):
        V[M+1][i] = Max(S[i]-K,0)
    for j in range(M+1):
        V[j][N+1] = L - K*exp(-r*(T-t[j]))
        V[j][0] = 0
    for n in range(M+1,0,-1):
        for i in range(1,N+1):
            V[n-1][i] = V[n][i] + delta_t*(r*S[i]*((V[n][i+1] - V[n][i-1])/(2*delta_s)) + 0.5*pow(sigma,2)*pow(S[i],2)*((V[n][i+1] + V[n][i-1] -2*V[n][i])/pow(delta_s,2)) - r*V[n][i])
    return V



def Graph_Dirichlet_DF():
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N+2)
    t2 = np.linspace(0,T,M+2)
    V1 = BS_Dirichlet_Call_Euro()
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(t2, S2, V1, rstride=1, cstride=1, cmap='jet')
    show()



def d1(S2,K1,sigma1,r1,t2,T1):
    return (log(S2/K1) + (r1 + (pow(sigma1,2)/2))*(T1 - t2))/(sigma*sqrt(T1 - t2))


def d2(S2,K1,sigma1,r1,t2,T1):
    return (log(S2/K1) + (r1 - (pow(sigma1,2)/2))*(T1 - t2))/(sigma*sqrt(T1 - t2))


def ka(x):
    return 1/(1+gamma*x)


def Nn(x):
    if x >= 0:
        N =0
        kaa = ka(x)
        N = pow(e,(-pow(x,2)/2))/sqrt(2*pi)
        N = N*(a1*kaa + a2*pow(kaa,2) + a3*pow(kaa,3) + a4*pow(kaa,4) + a5*pow(kaa,5))
        return 1 - N
    return 1 - Nn(-x)




def BS_Th_Call(S2,K1,sigma1,r1,t2,T1): 
    if t2 == T1:
        return Max(S2 - K1, 0)
    else :
        return S2*Nn(d1(S2,K1,sigma1,r1,t2,T1)) - K1*pow(e,-r1*(T1 - t2))*Nn(d2(S2,K1,sigma1,r1,t2,T1))



def BS_Th_Tab_Call(K1,sigma1,r1,T1):
    T = []
    T += [0]
    for i in range(1, N+2):
        T +=  [BS_Th_Call(S[i],K1,sigma1,r1,0,T1)]
    return T


def Graph_BS_Th_Call():
    S0_Tab = []
    Pay_off = []
    for i in range(N+2):
        Pay_off += [PayOff_Call_Euro(K, S[i])]
    Price = BS_Th_Tab_Call(K,sigma,r,T)
    plt.plot(S, Price)
    plt.plot(S, Pay_off)
    plt.show()

def BS_Neumann_Call_Euro():
    V = np.zeros((M+2,N+2))
    for i in range(N+2):
        V[M+1][i] = Max(S[i]-K,0)
    for n in range(M+1,0,-1):
        for i in range(1,N+1):
            V[n-1][i] = V[n][i] + delta_t*(r*S[i]*((V[n][i+1] - V[n][i-1])/(2*delta_s)) + 0.5*pow(sigma,2)*pow(S[i],2)*((V[n][i+1] + V[n][i-1] -2*V[n][i])/pow(delta_s,2)) - r*V[n][i])
        V[n-1][0] = V[n-1][1]
        V[n-1][N+1] = V[n-1][N] + delta_s
    
    return V
    


def Graph_Neumann_Call_DF():
    V1 = BS_Neumann_Call_Euro()
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N+2)
    t2 = np.linspace(0,T,M+2)
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(t2, S2, V1, rstride=1, cstride=1, cmap='jet')
    show()



def BS_Neumann_Put_Euro():
    V = np.zeros((M+2,N+2))
    for i in range(N+2):
        V[M+1][i] = Max(K - S[i],0)
    for n in range(M+1,0,-1):
        for i in range(1,N+1):
            V[n-1][i] = V[n][i] + delta_t*(r*S[i]*((V[n][i+1] - V[n][i-1])/(2*delta_s)) + 0.5*pow(sigma,2)*pow(S[i],2)*((V[n][i+1] + V[n][i-1] -2*V[n][i])/pow(delta_s,2)) - r*V[n][i])
        V[n-1][0] = V[n-1][1] - delta_s
        V[n-1][N+1] = V[n-1][N] 
    
    return V
    



def Graph_Neumann_Put_DF():
    V1 = BS_Neumann_Put_Euro()
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N+2)
    t2 = np.linspace(0,T,M+2)
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(S2, t2, V1, rstride=1, cstride=1, cmap='jet')
    show()


def PayOff_Butterfly(k,s):
    if s <= k:
        return 0
    elif (s > k and s <= 2*k):
        return s - k
    elif (s > 2*k and s <= 3*k):
        return 3*k - s
    else:
        return 0



def BS_Neumann_Butterfly_Euro():
    V = np.zeros((M+2,N+2))
    for i in range(N+2):
        V[M+1][i] = PayOff_Butterfly(K1,S[i])
    for n in range(M+1,0,-1):
        for i in range(1,N+1):
            V[n-1][i] = V[n][i] + delta_t*(r*S[i]*((V[n][i+1] - V[n][i-1])/(2*delta_s)) + 0.5*pow(sigma,2)*pow(S[i],2)*((V[n][i+1] + V[n][i-1] -2*V[n][i])/pow(delta_s,2)) - r*V[n][i])
        V[n-1][0] = V[n-1][1]
        V[n-1][N+1] = V[n-1][N] 
    return V


def Graph_Neumann_Butterfly_DF():
    V1 = BS_Neumann_Butterfly_Euro()
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N+2)
    t2 = np.linspace(0,T,M+2)
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(S2, t2, V1, rstride=1, cstride=1, cmap='jet')
    show()


def BS_Put_American():
    V = np.zeros((M+2,N+2))
    for i in range(N+2):
        V[M+1][i] = Max(K-S[i],0)
    for n in range(M+1,0,-1):
        for i in range(1,N+1):
            V[n-1][i] = Max(V[n][i+1]*0.5*delta_t*(pow(sigma,2)*(pow(S[i],2)/pow(delta_s,2)) + r*(S[i]/delta_s)) + V[n][i]*(1 - delta_t*(pow(sigma,2)*(pow(S[i],2)/pow(delta_s,2)) + r)) +V[n][i-1]*0.5*delta_t*(pow(sigma,2)*(pow(S[i],2)/pow(delta_s,2)) -r*(S[i]/delta_s)) , Max(K-S[i],0))
        V[n-1][0] = V[n-1][1]
        V[n-1][N+1] = V[n-1][N]     
    return V


def Graph_Put_American_DF():
    V1 = BS_Put_American()
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N+2)
    t2 = np.linspace(0,T,M+2)
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(S2, t2, V1, rstride=1, cstride=1, cmap='jet')
    show()


        
#--------------------------------------------------------------------------------------------------


#---------------------------- Methode de Monte Carlo --------------------------------------




def Prix_Call_S0_Fix(S0):
    somme = 0
    for i in range(Nme):
        S = S0*exp((r-(pow(sigma,2)/2))*T + sigma*sqrt(T)*np.random.randn(1)[0])
        gain = PayOff_Call_Euro(K,S)
        somme += gain
    return somme*exp(-r*T)/Nme



def Prix_Call_S0_Fix_Tab():
    Price = []
    for i in range(N1+2):
        Price += [Prix_Call_S0_Fix(S1[i])]
    return Price



def Graph_Price_Call_MC():
    S0_Tab = []
    Pay_off = []
    for i in range(N1+2):
        Pay_off += [PayOff_Call_Euro(K, S1[i])]
    Price = Prix_Call_S0_Fix_Tab()
    plt.plot(S1, Price)
    plt.plot(S1, Pay_off)
    plt.show()

#print(Graph_Price_Call_MC())

def Prix_Call_Stfixe_tfixe(t,st):
    gain = 0
    for i in range(Nme):
        sT = st*exp((r-(pow(sigma,2)/2))*(T - t) + sigma*sqrt(T-t)*np.random.randn(1)[0])
        gain += PayOff_Call_Euro(K,sT)
    price = exp(-r*(T-t))*gain/Nme
    return price

def Call_Euro_MC():
    V = np.zeros((M1+2,N1+2))
    for i in range(M1+2):
        for j in range(N1+2):
            V[i][j]= Prix_Call_Stfixe_tfixe(t1[i],S1[j])
    return V





def Graph_Call_MC():
    fig = figure()
    ax = Axes3D(fig)
    S2 = np.linspace(0,L,N1+2)
    t2 = np.linspace(0,T,M1+2)
    V1 = Call_Euro_MC()
    S2, t2 = np.meshgrid(S2, t2)
    ax.plot_surface(t2, S2, V1, rstride=1, cstride=1, cmap='jet')
    show()


#print(Graph_Call_MC())


#--------------------- Reduction de la variance --------------------------------


def Prix_Call_S0_Fix_Tab1():
    Price = []
    for i in range(N+2):
        Price += [Prix_Call_S0_Fix(S[i])]
    return Price



def Prix_Call_S0_Fix_Reduit(S0):
    somme = 0
    for i in range(Nme):
        wt = np.random.randn(1)[0]
        S = PayOff_Call_Euro(K,S0*exp((r-(pow(sigma,2)/2))*T + sigma*sqrt(T)*wt)) + PayOff_Call_Euro(K,S0*exp((r-(pow(sigma,2)/2))*T - sigma*sqrt(T)*wt))
        S = S/2
        somme += S
    return somme*exp(-r*T)/Nme

def Prix_Call_S0_Fix_Reduit_Tab():
    V = []
    for i in range(N+2):
        V += [Prix_Call_S0_Fix_Reduit(S[i])]
    return V

def Graph_Call_MC_1():
    Price1 = Prix_Call_S0_Fix_Tab1()
    Pay_off = []
    for i in range(N+2):
        Pay_off += [PayOff_Call_Euro(K, S[i])]
    plt.plot(S, Price1)
    plt.plot(S, Pay_off)
    plt.show()
    

def Graph_Prix_Compar_S0_Fix_Tab():
    Price1 = Prix_Call_S0_Fix_Tab1()
    Price2 = Prix_Call_S0_Fix_Reduit_Tab()
    Pay_off = []
    for i in range(N+2):
        Pay_off += [PayOff_Call_Euro(K, S[i])]
    plt.plot(S, Price1)
    plt.plot(S, Price2)
    plt.plot(S, Pay_off)
    plt.show()


#print(Graph_Prix_Compar_S0_Fix_Tab())


#def Calcul_Interv_Confiance():
    

#-------------------------------------------------------------------------------


#------------------------ Interface Graphique ---------------------------------

def open_window1():
    window1 = Tk()
    window1.title("Simulation par la methode de Monte Carlo")
    window1.geometry("800x800")
    window1.config(background='#41B77F')
    fn11 = Button(window1, text ="Call Européenne en 3D", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Call_MC)
    fn11.pack(pady=25)
    fn11 = Button(window1, text ="Call Européenne en 2D", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Call_MC_1)
    fn11.pack(pady=25)
    fn12 = Button(window1, text ="Call Européenne avec réduction de variance", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Prix_Compar_S0_Fix_Tab)
    fn12.pack(pady=25)


def open_window2():
    window2 = Tk()
    window2.title("Simulation par la methode de différences finies")
    window2.geometry("800x800")
    window2.config(background='#41B77F')
    fn21 = Button(window2, text ="call Européenne avec conditions de Dirichlet", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Dirichlet_DF)
    fn21.pack(pady=25)
    fn22 = Button(window2, text ="call Européenne avec conditions de Neumann", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Neumann_Call_DF)
    fn22.pack(pady=25)
    fn23 = Button(window2, text ="Put Européenne avec conditions de Neumann", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Neumann_Put_DF)
    fn23.pack(pady=25)
    fn24 = Button(window2, text ="Butterfly Européenne avec conditions de Neumann", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Neumann_Butterfly_DF)
    fn24.pack(pady=25)
    fn25 = Button(window2, text ="Put Américain avec conditions de Neumann", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_Put_American_DF)
    fn25.pack(pady=25)



    
window = Tk()

window.title("Financial Models")
window.geometry("1800x1000")
window.config(background='#41B77F')
label_title = Label(window, text="Bienvenue sur l'application de modélisation financière", font=("Courrier", 40), bg='#41B77F', fg='White')
label_title.pack()


fn1 = Button(window, text ="Simulation par la methode de Monte Carlo", font=("Courrier", 25), bg='white', fg='#41B77F',   command = open_window1)

fn2 = Button(window, text ="Simulation par la methode de différences finies", font=("Courrier", 25), bg='white', fg='#41B77F',   command = open_window2)

fn3 = Button(window, text ="Simulation par la methode analytique", font=("Courrier", 25), bg='white', fg='#41B77F',   command = Graph_BS_Th_Call)

fn1.pack(pady=25)
fn2.pack(pady=25)
fn3.pack(pady=25)

window.mainloop()



#------------------------------------------------------------------------------


def Calcul_Var():
    Var = 0
    S1 = 0
    for i in range(Nme):
        S = PayOff_Call_Euro(K,S0*exp((r-(pow(sigma,2)/2))*T + sigma*sqrt(T)*np.random.randn(1)[0]))
        S1 += S/Nme 
        Var += pow(S,2)/Nme
    Var -= pow(S1,2)
    return Var

print("La variance de S est : ", Calcul_Var())

Lp=0
Lf=0

for i in range(Nme):
    Lp += PayOff_Call_Euro(K,S0*exp((r-(pow(sigma,2)/2))*T + sigma*sqrt(T)*np.random.randn(1)[0]))
    Lf += PayOff_Call_Euro(K,S0*exp((r-(pow(sigma,2)/2))*T + sigma*sqrt(T)*np.random.randn(1)[0]))
Lp = Lp/Nme
Lf = Lf/Nme
Lp += (1.96*sqrt(Calcul_Var()))/Nme
Lf -= (1.96*sqrt(Calcul_Var()))/Nme

print("L'intervalle de confiance de S est : [ ", Lf, " , ", Lp ," ]")
    
