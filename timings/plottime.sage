# plottime.sage

Lmagma = []
f = open("TimeMagmaHMF.log","r");
for line in f:
  v = line.strip().split(' ')
  N = ZZ(v[0])
  Psi = ZZ(v[1])
  tmagma = RR(v[2])
  if N > 200:
    Lmagma.append((N,tmagma))

Pmagma = points(Lmagma, scale="loglog", axes_labels=['level norm', 'time (seconds)'], legend_label="magma", color='red', marker='x', size=40)

Lfdom = []
Loneword = []
f = open("TimeCohom.log","r");
for line in f:
  v = line.strip().strip('[]').split(',')
  N = ZZ(v[0])
  Psi = ZZ(v[1])
  tfdom = (ZZ(v[2])*0.001)
  toneword = (ZZ(v[3])*0.001)
  if N > 200:
    Lfdom.append((N,tfdom))
    Loneword.append((N,toneword))

Pfdom = points(Lfdom, legend_label="fdom", color='blue')
Poneword = points(Loneword, legend_label="oneword", color='green', marker="+", size=40, frame=True)

(Pmagma + Pfdom + Poneword).save("time.png")

var('a,b')
model(x) = a+b*x
print("magma", find_fit([[log(p[0]),log(p[1])] for p in Lmagma],model))
print("fdom", find_fit([[log(p[0]),log(p[1])] for p in Lfdom],model))
print("oneword", find_fit([[log(p[0]),log(p[1])] for p in Loneword],model))

