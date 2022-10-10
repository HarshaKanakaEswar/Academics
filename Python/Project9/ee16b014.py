from pylab import*
import mpl_toolkits.mplot3d.axes3d as p3
t=linspace(-32*pi,32*pi,2049);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(2048)
wnd=fftshift(0.54+0.46*cos(2*pi*n/2047))
y=((cos(0.86*t))**3)# y=sin(1.25*t)

y[0]=0 
y1=fftshift(y) 
Y1=fftshift(fft(y))/2048.0

y=y*wnd
y=fftshift(y)
Y=fftshift(fft(y))/2048.0
w=linspace(-pi*fmax,pi*fmax,2049);w=w[:-1]
figure()

subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-8,8])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^3\left(0.86t\right)\times w(t)$""\n with hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)

figure()
subplot(2,1,1)
plot(w,abs(Y1),'b',lw=2)
xlim([-8,8])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^3\left(0.86t\right)\times w(t)$""\n without hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y1),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)







t=linspace(-pi,pi,129);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(128)
wnd=fftshift(0.54+0.46*cos(2*pi*n/127))
w0=1.34
d=pi/4
y=cos(w0*t+d)
y=y*wnd
y[0]=0
y=fftshift(y) 
Y=fftshift(fft(y))/128.0
w=linspace(-pi*fmax,pi*fmax,129);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b',w,abs(Y),'bo',lw=2)
xlim([-16,16])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos(w_0+\delta)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-16,16])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)




t=linspace(-pi,pi,129);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(128)
wnd=fftshift(0.54+0.46*cos(2*pi*n/128))
w0=1.34
d=pi/4
y=cos(w0*t+d)+0.1*randn(128)
y=y*wnd
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/128.0
w=linspace(-pi*fmax,pi*fmax,129);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b',w,abs(Y),'bo',lw=2)
xlim([-16,16])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos(w_0+\delta)\times w(t)$ with white gaussian noise")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-16,16])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)




t=linspace(-pi,pi,1025);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(1024)
wnd=fftshift(0.54+0.46*cos(2*pi*n/1024))

y=cos(16*(1.5+t/(2*pi))*t)# y=sin(1.25*t)
y=y*wnd
y[0]=0
y=fftshift(y) 
Y=fftshift(fft(y))/1024.0
w=linspace(-pi*fmax,pi*fmax,1025);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b',lw=2)
xlim([-100,100])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos\left(16\left(1.5+\frac{t}{2\pi}\right)t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-100,100])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)





t=linspace(-pi,pi,1025);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
A=zeros((64,16))
t2=arange(16)

y=cos(16*(1.5+t/(2*pi))*t)
for t1 in t2:
	A[:,t1]=fftshift(fft(fftshift(y[t1*64:(t1+1)*64]*wnd)))/64.0
print(A)
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
x=arange(0,16,1)
y=arange(64)
X,Y=meshgrid(x,y)

fig1=figure()
ax=p3.Axes3D(fig1)
title("time-frequency plot")

for i in range(16):
	ax.plot3D(X[:,i],w,abs(A[:,i]))
ylim(-50,50)
xlabel("time step")
ylabel("w")
show()
