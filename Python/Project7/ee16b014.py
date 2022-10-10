from sympy import *
import pylab as p
from numpy import *
import scipy.signal as sp
s=symbols('s')

def lowpass(R1,R2,C1,C2,G,Vi,Ao):
       
        A= Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-Ao,Ao,1],[-(1/R1)-(1/R2)-(s*C1),1/R2,0,s*C1]])
        b= Matrix([0,0,0,-(Vi/(R1))])
        V=A.inv()*b
        return (A,b,V)

def coff(Vo):
	Vo=simplify(Vo)
	N,D=fraction(Vo)
	As=simplify(N)
	bs=simplify(D)
	As1=Poly(As,s)
	bs1=Poly(bs,s)
	Num=As1.all_coeffs()
	Den=bs1.all_coeffs()
	Num1=list(map(float,Num))
	Den1=list(map(float,Den)) 
	return Num1,Den1     
def u(t):
         return p.ones(len(t))     


A,b,V=lowpass(10**4,10**4,10**-11,10**-11,1.586,1,1000)
Vo=V[3]
w=p.logspace(0,8,801)
ss=1j*w
hf=lambdify(s,Vo,'numpy')
v=hf(ss)

p.loglog(w,abs(v))
p.xlabel('log(w)')
p.ylabel('log(|H(jw)|)')
p.title('magnitude of lowpass filter in loglog scale')
p.grid(True)
p.show()

   
Num1,Den1=coff(Vo)
H=sp.lti(Num1,Den1)
t=p.linspace(0,0.0001,3000)
                
t,y,svec=sp.lsim(H,u(t),t)
p.plot(t,y)
p.xlabel('t-->')
p.ylabel('Vo(t)')
p.title('step response for low pass filter')
p.show()

t1=p.linspace(0,0.00001,10000)
def x(t):
         return p.sin(2000*(p.pi)*t)+p.cos(2*1e10*(p.pi)*t)

t1,y2,svec=sp.lsim(H,x(t1),t1)
p.xlabel('t-->')
p.ylabel('Vo(t)')
p.title('Response of low pass filter for input sin(2000*(pi)*t)+cos(2*1e6*(pi)*t)')
p.plot(t1,y2)
p.show()


def highpass(R1,R3,C1,C2,G,Vi,Ao):
        A= Matrix([[0,0,1,-1/G],[-((s*C2*R3))/(1+(s*C2*R3)),1,0,0],[0,-Ao,Ao,1],[-(s*C1)-(s*C2)-(1/R1),s*C2,0,1/R1]])
        #A= Matrix([[0,0,1,-1/G],[-1/(1+1/(s*C2*R3),1,0,0],[0,-G,G,1],[-(s*C1)-(s*C2)-(1/R1),s*C2,0,1/R1]])
        b= Matrix([0,0,0,-(Vi*(s*C1))])
        V=A.inv()*b
        return (A,b,V)
       

A,b,V=highpass(10**4,10**4,10**-9,10**-9,1.586,1,1000)
Vo=V[3]
w=p.logspace(0,8,801)
ss=1j*w
hf=lambdify(s,Vo,'numpy')
v=hf(ss)
p.loglog(w,abs(v))
p.xlabel('log(w)')
p.ylabel('log(|H(jw)|)')
p.title('magnitude of highpass filter in loglog scale')
p.grid(True)
p.show()

Num1,Den1=coff(Vo)
H=sp.lti(Num1,Den1)
t=p.linspace(0,0.0001,300)

t,y,svec=sp.lsim(H,u(t),t)
p.plot(t,y)
p.xlabel('t-->')
p.ylabel('Vo(t)')
p.title('step response for highpass filter')
p.show()

def x(t):
         return (e**(-1000*t))*sin(100000*pi*t)
t1=p.linspace(0,0.0012,5000)
t1,y2,svec=sp.lsim(H,x(t1),t1)
p.xlabel('t-->')
p.ylabel('Vo(t)')
p.title('Response of highpass filter for input (e**(-10*t))*sin(200000*pi*t)')
p.plot(t1,y2)
p.show()


