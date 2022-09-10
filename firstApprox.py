import numpy as np
import math
import ephem
pi = math.pi
cos = np.cos; sin = np.sin; arccos = np.arccos; arcsin = np.arcsin
konst = 58.1324408700000; clight = 10065.320/konst; eps = 0.4089888210000
def sec(x):
    return 1/cos(x)

def getSunDist(time):
    s = ephem.Sun(time) 
    r = s.earth_distance; x = s.g_ra; y = s.g_dec
    sunDist  = [r * cos(y) * cos(x), r * cos(y) * sin(x), r * sin(y)]
    return sunDist
    
    
def getSunDist2(time):
    t = ephem.Date(time)
    s = ephem.Sun(time) 
    r = s.earth_distance;
    t += r/clight #yes, pyephem corrects itself for light travel time
    s = ephem.Sun(t); 
    r = s.earth_distance; x = s.g_ra; y = s.g_dec
    sunDist  = [r * cos(y) * cos(x), r * cos(y) * sin(x), r * sin(y)]
    return sunDist

def getRA(st):
    h, m, s = map(float, st.split(':'))
    return (h+m/60+s/3600)/12*pi

def getDEC(st):
    k = 1
    if st[0] == '-':
        k = -1
        st = st[1:]
    h, m, s = map(float, st.split(':'))
    return k * (h + m/60 + s/3600)/180 *pi
    
def updateSunDist(tt1, tt2, tt3):
    return (getSunDist(tt1), getSunDist(tt2), getSunDist(tt3))
        
def updateSunDist2(tt1, tt2, tt3):
    return (getSunDist2(tt1), getSunDist2(tt2), getSunDist2(tt3))


filename = input("Input file's name: ")
f = open(filename, 'r')
lat, lon = np.radians(np.array(list(map(float, f.readline().split()))))
coslat = cos(lat); sinlat = sin(lat)
tt11, tt12, obsstr11, obsstr12 = f.readline().split()
tt21, tt22, obsstr21, obsstr22 = f.readline().split()
tt31, tt32, obsstr31, obsstr32 = f.readline().split()
obs1, obs2, obs3 = (list(range(2)) for i in range(3))
sunDist1, sunDist2, sunDist3 = (list(range(3)) for i in range(3))
obs1[0] = getRA(obsstr11); obs1[1] =  getDEC(obsstr12)
obs2[0] = getRA(obsstr21); obs2[1] = getDEC(obsstr22)
obs3[0] = getRA(obsstr31); obs3[1] = getDEC(obsstr32)
times = [tt11+' ' +tt12, tt21 + ' ' +  tt22, tt31+ ' ' +tt32]
sunDist1, sunDist2, sunDist3 = updateSunDist(times[0], times[1], times[2])
tt1 = ephem.Date(times[0]); tt2 = ephem.Date(times[1]); tt3 = ephem.Date(times[2])
f.close()
t1 = (tt3 - tt2)/konst; t2 = (tt3 - tt1)/konst; t3 = (tt2 - tt1)/konst
b1 = t1/t2; b3 = t3/t2
a1 = b1; a3 = b3 #as a first approximation we equal the triangle ratios to sector ratios
Re = 6371/149600000; Rsinlat = Re * sinlat; Rcoslat = Re*coslat
def localST(t):
    obs = ephem.Observer()
    obs.lon = lon; obs.lat = lat; obs.date = t
    return obs.sidereal_time()
    
def topoCorr(t, x, y, z):
    lst = localST(t)/12 * pi
    x = x - Rcoslat*cos(lst)
    y = y - Rcoslat*sin(lst)
    z = z - Rsinlat
    return (x, y, z)
def calcAngle(sinx, cosx): #function that calculates angle taking into account both sine and cosine to get the correct quadrant
    if sinx > 0 and  cosx > 0:
        return arcsin(sinx)
    elif sinx > 0 and cosx < 0:
        return arccos(cosx)
    elif sinx < 0 and cosx < 0:
        return pi-arcsin(sinx)
    else: return arcsin(sinx)

def coeff(obs1): #calculates direction cosines from spherical coordinates' angles
    l = cos(obs1[0])*cos(obs1[1])
    m = sin(obs1[0])*cos(obs1[1])
    n = sin(obs1[1])
    return (l, m, n)

def koeff(a1, a3, obs1, obs2, obs3): 
    k1 = a1*obs1[0] - obs2[0] + a3*obs3[0]
    k2 = a1*obs1[1] - obs2[1] + a3*obs3[1]
    k3 = a1*obs1[2] - obs2[2] + a3*obs3[2]
    return[k1, k2, k3]

def arrs(a1, a3, c1, c2, c3, k): #structures variables into arrays for further calculation
    d0 = np.array([[c1[0]*a1, -c2[0], c3[0]*a3], 
                   [c1[1]*a1, -c2[1], c3[1]*a3],
                   [c1[2]*a1, -c2[2], c3[2]*a3]])
    d1 = np.array([[k[0], -c2[0], c3[0]*a3],            
                   [k[1], -c2[1], c3[1]*a3],
                   [k[2], -c2[2], c3[2]*a3]])
    d2 = np.array([[c1[0]*a1, k[0], c3[0]*a3],            
                   [c1[1]*a1, k[1], c3[1]*a3],
                   [c1[2]*a1, k[2], c3[2]*a3]])
    d3 = np.array([[c1[0]*a1, -c2[0], k[0]],            
                   [c1[1]*a1, -c2[1], k[1]],
                   [c1[2]*a1, -c2[2], k[2]]])
    return (d0, d1, d2, d3)

def deter(d): #returns determinant of a matrix
    return np.linalg.det(d)    
def sqsum(c): 
    return math.sqrt(c[0]**2 + c[1]**2 + c[2]**2)

def helDist(c, dist, sundist): #returns heliocentric equatorial coordinates from geocentric equatorial coordinates
    d1 = c*dist - sundist
    return d1
#these functions are defined for solving the equations for sector ratios
def fx(x):
    return 2*x
def gx(x):
    return x*(3*x - 2)
def fy(a, b, y):
    return a*sin(y)/((b - cos(y))**2)
def gy(a, b, y):
    siny = sin(y); cosy = cos(y); sin3y = siny**3; sin4y = sin3y*siny
    return a*(3*(y - siny*cosy)*cosy - 2*sin3y)/sin4y
def f(a, b, x, y):
    temp = x**2 - a/(b - cos(y))
    return temp
def g(a, b, x, y):
    siny = sin(y)
    temp = (x**2 - x)*x - a*(y - siny*cos(y))/(siny**3)
    return temp
def iterxy(a, b, x, y):
    fx_ = fx(x); fy_ = fy(a, b, y); gx_ = gx(x); gy_ = gy(a, b, y); f_ = f(a, b, x, y); g_ = g(a, b, x, y)
    denom = fx_*gy_-fy_*gx_
    h = (gy_*f_-fy_*g_)/denom
    k = (fx_*g_-gx_*f_)/denom
    return (x-h, y-k)


#this function returns heliocentric equatorial distances and their components + geocentric equatorial coordinates (cartesian coordinates)
def helioDist(a1, a3, sunDist1, sunDist2, sunDist3, t1, t2, t3):
    k = koeff(a1, a3, sunDist1, sunDist2, sunDist3)
    d0, d1, d2, d3 = arrs(a1, a3, c1, c2, c3, k)
    det0 = deter(d0)
    dist1 = deter(d1)/det0
    dist2 = deter(d2)/det0
    dist3 = deter(d3)/det0
    
    x1 = helDist(c1[0], dist1, sunDist1[0])
    x2 = helDist(c2[0], dist2, sunDist2[0])
    x3 = helDist(c3[0], dist3, sunDist3[0])
    y1 = helDist(c1[1], dist1, sunDist1[1])
    y2 = helDist(c2[1], dist2, sunDist2[1])
    y3 = helDist(c3[1], dist3, sunDist3[1])
    z1 = helDist(c1[2], dist1, sunDist1[2])
    z2 = helDist(c2[2], dist2, sunDist2[2])
    z3 = helDist(c3[2], dist3, sunDist3[2])
    r1 = sqsum([x1, y1, z1]); r2 = sqsum([x2, y2, z2]); r3 = sqsum([x3, y3, z3])
    return (x1,x2,x3,y1,y2,y3,z1,z2,z3,r1,r2,r3,dist1,dist2,dist3)

def triangleSectors(t1, t2, t3, r2):#returns triangle sector ratios
    b1 = t1/t2; b3 = t3/t2
    return (b1 + t1*t3*(1+b1)/(6*r2**3), b3 + t1*t3*(1+b3)/(6*r2**3))

def distSectRat(a1, a3, sunDist1, sunDist2, sunDist3, tt1, tt2, tt3, t1, t2, t3):
    
    x1, x2, x3, y1, y2, y3, z1, z2, z3, r1, r2, r3, dist1, dist2, dist3 = helioDist(a1, a3, sunDist1, sunDist2, sunDist3, t1, t2, t3)
    x1, y1, z1 = topoCorr(tt1, x1, y1, z1)
    x2, y2, z2 = topoCorr(tt2, x2, y2, z2)
    x3, y3, z3 = topoCorr(tt3, x3, y3, z3)
    #redoing r1, r2, r3
    r1 = sqsum([x1, y1, z1]); r2 = sqsum([x2, y2, z2]); r3 = sqsum([x3, y3, z3])
    
    f3 = 0.5 * arccos((x1*x2 + y1*y2+ z1*z2)/(r1*r2))
    f2 = 0.5 * arccos((x1*x3 + y1*y3 + z1*z3)/(r1*r3))
    f1 = 0.5 * arccos((x2*x3 + y2*y3 + z2*z3)/(r2*r3))
    cosf3 = cos(f3); cosf2 = cos(f2); cosf1 = cos(f1)

    M3 = t3/(2*math.sqrt((math.sqrt(r1*r2)*cosf3)**3))
    N3 = (r1 + r2)/(2*math.sqrt(r1*r2)*cosf3)
    M2 = t2/(2*math.sqrt((math.sqrt(r1*r3)*cosf2)**3))
    N2 = (r1 + r3)/(2*math.sqrt(r1*r3)*cosf2)
    M1 = t1/(2*math.sqrt((math.sqrt(r2*r3)*cosf1)**3))
    N1 = (r2 + r3)/(2*math.sqrt(r2*r3)*cosf1)

    R1 = 1; R2 = 1; R3 = 1;
    g1 = f1; g2 = f2; g3 = f3;
    m1sq = M1**2; m2sq = M2**2;m3sq = M3**2;
    for i in range(10):
        R1, g1 = iterxy(m1sq, N1, R1, g1)
        R2, g2 = iterxy(m2sq, N2, R2, g2)
        R3, g3 = iterxy(m3sq, N3, R3, g3)
    
    return x1,x2,x3,y1,y2,y3,z1,z2,z3,r1,r2,r3,dist1,dist2,dist3,R1,R2,R3

def cTime(tt1, tt2, tt3, dist1, dist2, dist3, R1, R2, R3): 
    #light travel time corrections, but preserving original times for later iterations
    newtt1 = tt1 - dist1/clight
    newtt2 = tt2 - dist2/clight
    newtt3 = tt3 - dist3/clight
    t1 = (newtt3 - newtt2)/konst
    t2 = (newtt3 - newtt1)/konst
    t3 = (newtt2 - newtt1)/konst
    b1 = t1/t2; b3 = t3/t2 #updating sector ratios
    s1, s2, s3 = updateSunDist2(newtt1, newtt2, newtt3)
    a1 = R2/R1 * b1 #calculating "real" triangle ratios
    a3 = R2/R3 * b3 #calculating "real" triangle ratios
    return a1, a3, t1, t2, t3, s1, s2, s3

def printElems(x1,x2,x3,y1,y2,y3,z1,z2,z3, r1, r2, r3, R1, t3):
    f3 = 0.5 * arccos((x1*x2 + y1*y2+ z1*z2)/(r1*r2))
    f2 = 0.5 * arccos((x1*x3 + y1*y3 + z1*z3)/(r1*r3))
    f1 = 0.5 * arccos((x2*x3 + y2*y3 + z2*z3)/(r2*r3))
    l = ((r1*r2*sin(2*f3)*R1)/t3)**2

    ecosv1 = l/r1 -1
    temp = (l /r1 - 1) * cos (2 * f2) - l/r3 + 1
    esinv1 = temp/sin(2*f2)

    e = math.sqrt(ecosv1**2+esinv1**2)
    sinv1 = esinv1/e
    cosv1 = ecosv1/e
    v1 = calcAngle(sinv1, cosv1)
    v2 = v1 + 2*f3
    v3 = v1 + 2*f2
    a = l /(1-e**2)

    px = (x1 * r3*sin(v3) - x3*r1*sin(v1))/(r1 * r3 *sin(2*f2))
    qx = (x3 * r1 *cos(v1) - x1 * r3*cos(v3))/(r1 * r3 *sin(2*f2))
    py = (y1 *r3 *sin(v3) - y3*r1*sin(v1))/(r1 * r3 *sin(2*f2))
    qy = (y3 * r1 * cos(v1) - y1 * r3 * cos(v3))/(r1 * r3 *sin(2*f2))
    pz = (z1 * r3 * sin(v3) - z3*r1*sin(v1))/(r1 * r3 *sin(2*f2))
    qz = (z3 * r1 * cos(v1) - z1 * r3 * cos(v3))/(r1 * r3 *sin(2*f2))
    sineps = sin(eps); coseps = cos(eps)
    sinwsini = pz*coseps - py*sineps
    coswsini = qz*coseps - qy*sineps
    w = np.arctan(sinwsini/coswsini)
    sinw = sin(w); cosw = cos(w); sini = sinwsini/sinw
    sinom = (py*cosw-qy*sinw)/coseps
    cosom = px*cosw-qx*sinw
    om = calcAngle(sinom, cosom)
    cosi=-1*(px*sinw + qx*cosw)/sinom
    i = calcAngle(sini, cosi)
    if sini < 0:
    	w = np.arctan(sinwsini/coswsini) + pi
    	sinw = sin(w); cosw = cos(w); 
    	sinom = (py*cosw-qy*sinw)/coseps
    	cosom = px*cosw-qx*sinw
    	om = calcAngle(sinom, cosom)
    	cosi=-1*(px*sinw + qx*cosw)/sinom
    	sini = -sini
    	i = calcAngle(sini, cosi)
    	    
    if w < 0: w += 2*pi
    P = math.sqrt(a**3) 
    i = np.degrees(i)
    w = np.degrees(w)
    om = np.degrees(om)
    cosE = (e+cos(v1))/(1+e*cos(v1))
    E = arccos(cosE)
    sinE = sin(E)
    T = tt1 - 365.25 *  P/(2*pi)*(E - e * sinE)
    T = round(T, 0) + 0.5
    dateString = str(ephem.Date(T)).split()[0]
    print('Orbital elements are:')
    print(f'a = {round(a, 3)}AU   also: P = {round(P, 2)} years          i = {round(i, 2)}\u00b0')
    print(f'e = {round(e, 4)}             om = {round(om, 1)}\u00b0')
    print(f'T = {dateString}        w = {round(w, 1)}\u00b0')
    

#getting diretion cosines in geocentric equatorial coordinate system for all three observations
c1 = coeff(obs1); c2 = coeff(obs2); c3 = coeff(obs3)
#iterating few times for heliocentric distances with approximated triangle ratios
a1 = b1; a3 = b3
for _ in range(15):
    x1, x2, x3, y1, y2, y3, z1, z2, z3, r1, r2, r3, dist1, dist2, dist3 = helioDist(a1, a3, sunDist1, sunDist2, sunDist3, t1, t2, t3)
    a1, a3 = triangleSectors(t1, t2, t3, r2)


#light travel time corrections, but preserving original times for later iterations
newtt1 = tt1 - dist1/clight
newtt2 = tt2 - dist2/clight
newtt3 = tt3 - dist3/clight
t1 = (newtt3 - newtt2)/konst
t2 = (newtt3 - newtt1)/konst
t3 = (newtt2 - newtt1)/konst
b1 = t1/t2; b3 = t3/t2 #updating sector ratios
a1, a3 = triangleSectors(t1, t2, t3, r2) #updating triangle ratios with updated times
sunDist1, sunDist2, sunDist3 = updateSunDist2(newtt1, newtt2, newtt3)

#iterating few times for heliocentric distances with sector ratios
for _ in range(20):
    x1,x2,x3,y1,y2,y3,z1,z2,z3,r1,r2,r3,dist1,dist2,dist3,R1,R2,R3= distSectRat(a1,a3,sunDist1,sunDist2,sunDist3,tt1,tt2,tt3,t1,t2,t3)
    a1, a3, t1, t2, t3, s1, s2, s3 = cTime(tt1, tt2, tt3, dist1, dist2, dist3, R1, R2, R3)
    sunDist1 = s1; sunDist2 = s2; sunDist3 = s3
	
printElems(x1, x2, x3, y1, y2, y3, z1, z2, z3, r1, r2, r3, R1, t3)
