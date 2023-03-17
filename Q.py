from sympy import *
import numpy as np
import matplotlib.pyplot as plt
init_printing(use_unicode = true)

def sumUpowN(n,p,k,H):
    t = n*p
    s1 = sin(k*H/2)**(t)
    s2 = ((-1)**(t)/factorial(t))*diff(cot(k),k,t-1)
    return s1*s2

p = 3
k, H, a = symbols("k H a")
D = I * sin(k*H)/H
USquared = sin(k*H/2)/(k*H/2)
Dfunc = lambdify([k,H],D,"numpy")
sumUSquaredSymbol = sumUpowN(2,2,k,H)
sumUPow4 = sumUpowN(4,2,k,H)
S2 = (12/(k*a/2)**4)*(2-2*cos(k*a/2)-(k*a)/2*sin(k*a/2))
R = -I*S2**2 / k
# sumUSquaredFunction = lambdify([k,H],sumUSquaredSymbol,modules = ["numpy",{"cot": cot}])
# print(sumUSquaredFunction(1,1))
# print(sumUSquaredSymbol)
# print(simplify(sumUSquaredSymbol))
cotFunc = lambda x: 1/tan(x)
sumUSquaredSymbol = simplify(sumUSquaredSymbol)
sumUSquaredSymbol = sumUSquaredSymbol.subs(cot(k),1/tan(k))
sumUSquaredFunc = lambdify([k,H],sumUSquaredSymbol,modules = ["numpy",{"cot":cotFunc}])
# sumUSquaredFunc = lambdify()
sumUPow4Symbol = simplify(sumUPow4)
sumUPow4Func = lambdify([k,H],sumUPow4Symbol,"numpy")
USquaredFunc = lambdify([k,H],USquared,modules = ["numpy",{"cot":cotFunc}])

S2func = lambdify([k,a],S2)
Rfunc = lambdify([k,a],R)
print(sumUSquaredFunc(1,2))
print(sumUPow4Func(1,1))
print(Dfunc(1,1))
print(Rfunc(1,1))

def G(k,H,a):
    d = Dfunc(k,H)
    acc = 0
    for i in range(1,100):
        ktemp = 2*np.pi/i
        usquared = USquaredFunc(ktemp,H)
        r = Rfunc(ktemp,a)
        acc += r*usquared
    dSquared = np.abs(d) * np.abs(d)
    uSquared = sumUSquaredFunc(i,H)
    uSquaredSquared = uSquared * uSquared
    return (d*(acc)) / (dSquared * uSquaredSquared)

# print(G(1,1,1))
domain = []
yres = []
Zres = []
Qres = []
h = 5
for re in np.linspace(10,20,200):
    a = re*h/0.7
    P = 0
    # h = 0.7*a/re
    Z = 0
    for i in range(1,100):
        ktemp = 2*np.pi / i 
        RUacc = 0
        RsquaredAcc = 0
        for j in range(1,100):
            ktemp2 = 2*np.pi/j
            r = Rfunc(ktemp2,a)
            Usquared = USquaredFunc(ktemp2,h)
            RUacc += r*Usquared
            RsquaredAcc += np.abs(r)**2

        Ptemp = np.abs(np.abs(Dfunc(ktemp,h))**2)*G(ktemp,h,a)*(sumUSquaredFunc(ktemp,h)**2 - sumUPow4Func(ktemp,h))
        Z += np.abs(np.abs(Dfunc(ktemp,h))**2 * G(ktemp,h,a) * sumUPow4Func(ktemp,h) - 2*G(ktemp,h,a)*Dfunc(ktemp,h)*RUacc + RsquaredAcc)
        P += Ptemp
    Z /= 101**3
    P /= 101**3
    yres.append(np.sqrt(P)/1e88)
    Zres.append(np.sqrt(Z)/1e88)
    Qres.append(np.sqrt(P+Z)/1e88)
    domain.append(re/h)

plt.plot(domain,Qres,label = "$\sqrt{|Q|} p = 2$")
plt.plot(domain,yres,label = "$\sqrt{|P|} p = 2$")
plt.plot(domain,Zres,label = "$\sqrt{|Z|} p = 2$")



k, H, a = symbols("k H a")
D = I * sin(k*H)/H
USquared = sin(k*H/2)/(k*H/2)
Dfunc = lambdify([k,H],D,"numpy")
sumUSquaredSymbol = sumUpowN(2,1,k,H)
sumUPow4 = sumUpowN(4,1,k,H)
S2 = (12/(k*a/2)**4)*(2-2*cos(k*a/2)-(k*a)/2*sin(k*a/2))
R = -I*S2**2 / k
# sumUSquaredFunction = lambdify([k,H],sumUSquaredSymbol,modules = ["numpy",{"cot": cot}])
# print(sumUSquaredFunction(1,1))
# print(sumUSquaredSymbol)
# print(simplify(sumUSquaredSymbol))
cotFunc = lambda x: 1/tan(x)
sumUSquaredSymbol = simplify(sumUSquaredSymbol)
sumUSquaredSymbol = sumUSquaredSymbol.subs(cot(k),1/tan(k))
sumUSquaredFunc = lambdify([k,H],sumUSquaredSymbol,modules = ["numpy",{"cot":cotFunc}])
# sumUSquaredFunc = lambdify()
sumUPow4Symbol = simplify(sumUPow4)
sumUPow4Symbol = sumUPow4Symbol.subs(cot(k),1/tan(k))
sumUPow4Func = lambdify([k,H],sumUPow4Symbol,"numpy")
USquaredFunc = lambdify([k,H],USquared,modules = ["numpy",{"cot":cotFunc}])

S2func = lambdify([k,a],S2)
Rfunc = lambdify([k,a],R)
print(sumUSquaredFunc(1,2))
print(sumUPow4Func(1,1))
print(Dfunc(1,1))
print(Rfunc(1,1))

def G(k,H,a):
    d = Dfunc(k,H)
    acc = 0
    for i in range(1,100):
        ktemp = 2*np.pi/i
        usquared = USquaredFunc(ktemp,H)
        r = Rfunc(ktemp,a)
        acc += r*usquared
    dSquared = np.abs(d) * np.abs(d)
    uSquared = sumUSquaredFunc(i,H)
    uSquaredSquared = uSquared * uSquared
    return (d*(acc)) / (dSquared * uSquaredSquared)

# print(G(1,1,1))
domain = []
yres = []
Zres = []
Qres = []
h = 5
for re in np.linspace(10,20,200):
    a = re*h/0.7
    P = 0
    # h = 0.7*a/re
    Z = 0
    for i in range(1,100):
        ktemp = 2*np.pi / i 
        RUacc = 0
        RsquaredAcc = 0
        for j in range(1,100):
            ktemp2 = 2*np.pi/j
            r = Rfunc(ktemp2,a)
            Usquared = USquaredFunc(ktemp2,h)
            RUacc += r*Usquared
            RsquaredAcc += np.abs(r)**2

        Ptemp = np.abs(np.abs(Dfunc(ktemp,h))**2)*G(ktemp,h,a)*(sumUSquaredFunc(ktemp,h)**2 - sumUPow4Func(ktemp,h))
        Z += np.abs(np.abs(Dfunc(ktemp,h))**2 * G(ktemp,h,a) * sumUPow4Func(ktemp,h) - 2*G(ktemp,h,a)*Dfunc(ktemp,h)*RUacc + RsquaredAcc)
        P += Ptemp
    Z /= 101**3
    P /= 101**3
    yres.append(np.sqrt(P)/1e88)
    Zres.append(np.sqrt(Z)/1e88)
    Qres.append(np.sqrt(P+Z)/1e88)
    domain.append(re/h)

plt.plot(domain,Qres,label = "$\sqrt{|Q|} p = 1$")
plt.plot(domain,yres,label = "$\sqrt{|P|} p = 1$")
plt.plot(domain,Zres,label = "$\sqrt{|Z|} p = 1$")

plt.xlabel("$r_e/a$")
plt.ylabel("Percent")
plt.legend()
plt.savefig("P_Z_Errors.pdf",bbox_inches = "tight")
plt.show()


#TODO: Create function that evaluates G hat (k) Then compute and plot P, Z, Q = P+Z






    
