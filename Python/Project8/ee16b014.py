from pylab import*
x=linspace(-4*pi,4*pi,513);
x=x[:-1]

y=sin(x)**3
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
ylim([0,0.75])
title(r"Spectrum of $\sin^3(t)$")
grid(True)
#ax.get_xticks(range(24))
xticks(np.arange(-10, 11, 1.0))
yticks(np.arange(0,0.837,0.1875))

subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
xticks(np.arange(-10, 11, 1.0))
show()




y=cos(x)**3
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^3(t)$")
grid(True)
#ax.get_xticks(range(24))
xticks(np.arange(-10, 11, 1.0))


subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
xticks(np.arange(-10, 11, 1.0))
show()




x=linspace(-8*pi,8*pi,1025);
x=x[:-1]

y=cos(20*x+5*cos(x))
Y=fftshift(fft(y))/1024.0
w=linspace(-64,64,1025)
w=w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-30,30])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos(20t+5cost)$")
grid(True)


subplot(2,1,2)
#plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-30,30])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
savefig("fig9-2.png")
show()   #/sqrt(2*pi)




x=linspace(-4*pi,4*pi,513);
x=x[:-1]


y=e**(-x*x/2.0)
Y=4*fftshift(fft(fftshift(y)))/512
w=linspace(-64,64,513)
w=w[:-1]

figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $e^\frac{-t^2}{2}}$")
grid(True)


subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-6)
#plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-8,8])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
savefig("fig9-2.png")
show()   #/sqrt(2*pi)
