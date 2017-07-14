import math

def sign(x):
    return float(abs(x))/x

# Calculates flow rate given flow coefficient and pressure differential
def Q(Cv, P_a, P_b, T):

    sign = 1
    
    if P_b > P_a:
        P_a, P_b = P_b, P_a
        sign = -1
    
    B = 0.3565 # mol / s
    
    if P_a > 2*P_b:
        return sign * B * 0.471 * Cv * P_a * math.sqrt(1.0/T)
    else:
        return sign * B * Cv * P_a * (1-2.0/3.0*(1-P_b/P_a)) * math.sqrt((1-P_b/P_a)/T)

# Calculates exit velocity based on tank pressure, volume, t-shirt dimensions and flow
def fire(P_0, V_0, dia, W, L, d, theta, f, Cv):

    P_atm = 14.70 # psi
    R = 1.2059 # (psi * L) / (mol * K)
    T = 293 # K
    dt = 0.00001 # s
    
    N_0 = P_0 * V_0 / (R * T) # mol

    A = dia * dia / 4 * math.pi # in^2

    d = max(1.0/32.0, d)

    m = W * 0.031081

    P_a = P_0
    N_a = N_0
    
    P_b = P_atm
    N_b = P_b * (A * d * 0.0163871) / (R * T)

    pos = 0
    vel = 0
    accel = 0

    t = 0

    while pos * 12.0 + d < L:

        dN_dt = Q(Cv, P_a, P_b, T)

        N_a -= dN_dt * dt
        N_b += dN_dt * dt

        P_a = N_a * R * T / V_0
        P_b = N_b * R * T / (A * (d + pos * 12.0) * 0.0163871)

        force = A * (P_b - P_atm) - math.sin(theta) * W
        force = force - sign(force) * min(f, abs(force))

        accel = force / m

        vel += accel * dt
        pos += vel * dt
        t += dt

    return {'v': vel, 't': t, 'P': P_a}

# Calculates trajectory of an object accounting for quadratic drag
def trajectory(v_0, h_0, theta, C_d, A, W):

    dt = 0.00001 # s

    m = W * 0.031081

    v_x = v_0 * math.cos(theta) # m/s
    v_y = v_0 * math.sin(theta) # m/s

    x = 0 # m
    y = h_0 # m

    t = 0

    while y >= 0:

        v_x += (-1.0/2.0 * 0.0023768924 * C_d * A * (v_x**2 + v_y**2) * math.cos(theta))/m * dt
        v_y += (-1.0/2.0 * 0.0023768924 * C_d * A * (v_x**2 + v_y**2) * math.sin(theta))/m * dt - 32.174 * dt 

        x += v_x * dt
        y += v_y * dt

        theta = math.atan(v_y / v_x)
        
        t+=dt

    return {'x': x, 't': t}

# Test Case:
# 11L tank at 120psi, 16 Cv valve
# 0.42lb t-shirt wrapped at 2.375" diameter, 24" long barrel
# 45 degree shot angle, 0.90 drag coefficient

v = fire(120, 11, 2.375, 0.42, 24, 1, math.pi/4, 0, 16)['v']
x = trajectory(v, 2.46, math.pi/4, 0.90, 0.03043388, 0.4188783)['x']

print("Exit velocity (ft/s): " + str(round(v,2)))
print("Distance (ft): " + str(round(x,2)))
print("Distance (yds): " + str(round(x/3.0,2)))
