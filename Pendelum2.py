from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from manim import *

#Dane środowiskowe
g = 9.81
L = 1.0
b = 1
m = 1.0

time = 3600*10

y0 = [np.pi/5, np.pi/5,0,0]
y0_2 = [y0[0]+0.001, y0[1]+0.001,0,0]
y0_3 = [y0[0]+0.0001, y0[1]+0.0001,0,0]
y0_4 = [y0[0]+0.00001, y0[1]+0.00001,0,0]
y0_5 = [y0[0]+0.000001, y0[1]+0.000001,0,0]

print(f"theta_0: {y0[0]*360/(2*np.pi)}, theta_0 {y0[1]*360/(2*np.pi)}")

def x1_fun(sol):
    return np.sin(sol.y[0])*L

def y1_fun(sol):
    return -np.cos(sol.y[0])*L

def x2_fun(sol):
    return L*(np.sin(sol.y[0])+np.sin(sol.y[1]))

def y2_fun(sol):
    return -L*(np.cos(sol.y[0])+np.cos(sol.y[1]))
# theta0 = np.pi/2
# omega0 = 0.0

# t = np.linspace(0,time,100)

# print(f"theta0: {theta0}, omega0: {omega0}")
# print(f"b: {b}, m: {m}, g: {g}, L: {L}")

# def pendulum(t, y):
#     theta, omega = y
#     dydt = [omega, -b/m * omega - g/L * np.sin(theta)]
#     return dydt

# test_y = [0.1, 0]  # Example initial condition
# print(pendulum(0, test_y))  # Should return a list with two elements

# y0 = [theta0, omega0]
# sol = solve_ivp(fun=pendulum, t_span=(0,50), y0=y0, method='RK45')

# plt.plot(sol.t, sol.y[0], label="theta(t)")
# #plt.plot(sol.t, sol.y[1], label="omega(t)")
# plt.legend()
# plt.show()

def pendulemDual(t, y):
    theta1, theta2, p1, p2 = y

    d=(m*L**2)*(16-9*np.cos(theta1-theta2)**2)

    theta1_dot = (6*(2*p1 - 3*np.cos(theta1-theta2)*p2))/d
    theta2_dot = (6*(8*p2 - 3*np.cos(theta1-theta2)*p1))/d

    p1_dot = (-1/2)*(m*L**2)*(theta1_dot*theta2_dot*np.sin(theta1-theta2)+3*g*np.sin(theta1)/L)
    p2_dot = (-1/2)*(m*L**2)*(-theta1_dot*theta2_dot*np.sin(theta1-theta2)+g*np.sin(theta2)/L)
    
    return [theta1_dot, theta2_dot, p1_dot, p2_dot]

t_span = (0, time)  # od 0 do 10 sekund
t_eval = np.linspace(0, time, time*100)  # wartości czasu dla rozwiązania

# Rozwiązywanie układu równań
sol = solve_ivp(pendulemDual, t_span, y0, method='RK45', t_eval=t_eval)
sol2 = solve_ivp(pendulemDual,t_span, y0_2, method='RK45',t_eval=t_eval)
sol3 = solve_ivp(pendulemDual,t_span, y0_3, method='RK45',t_eval=t_eval)
sol4 = solve_ivp(pendulemDual,t_span, y0_4, method='RK45',t_eval=t_eval)
sol5 = solve_ivp(pendulemDual,t_span, y0_5, method='RK45',t_eval=t_eval)

def text(y0):
    return r"($\theta_1)_{t=0}$: "+f"{y0[0]:.6f}"

def text2(y0):
    return r"($\theta_2)_{t=0}$: "+f"{y0[1]:.6f}"

fig, axes = plt.subplots(2, 1, figsize=(10, 4))

def energia(sol):
    theta1, theta2, p1, p2 = sol.y
    v1_sq = (L * np.gradient(theta1, sol.t)) ** 2
    v2_sq = (L * np.gradient(theta2, sol.t)) ** 2
    T = 0.5 * m * (v1_sq + v2_sq)
    U = m * g * (L * (1 - np.cos(theta1)) + L * (1 - np.cos(theta2)))
    return U

E1 = energia(sol)
E2 = energia(sol2)
E3 = energia(sol3)
E4 = energia(sol4)
E5 = energia(sol5)

# axes[2].plot(sol.t, E1, label='Energia całkowita y0', color='b')
# axes[2].plot(sol2.t, E2, label='Energia całkowita y0_2', color='r')
# axes[2].plot(sol3.t, E3, label='Energia całkowita y0_3', color='pink')
# axes[2].plot(sol4.t, E4, label='Energia całkowita y0_4', color='purple')
# axes[2].plot(sol5.t, E5, label='Energia całkowita y0_5', color='green')
# axes[2].legend()
# axes[2].set_xlabel('Czas [s]')
# axes[2].set_ylabel('Energia [J]')

axes[0].plot(sol.t, sol.y[0], label=text(y0),color="b")
axes[0].set_title(r"Kąt $\theta_1$")
axes[0].set_xlabel("t[s]")  # Podpis osi X
axes[0].set_ylabel(r"$\theta_1$[rad]")  # Podpis osi Y


axes[1].plot(sol.t, sol.y[1], label=text2(y0),color="b")

axes[1].set_title(r"Kąt $\theta_2$")
axes[1].set_xlabel("t[s]")  # Podpis osi X
axes[1].set_ylabel(r"$\theta_2$[rad]")  # Podpis osi Y

axes[0].plot(sol2.t, sol2.y[0], label=text(y0_2), color="r")
axes[1].plot(sol2.t, sol2.y[1], label=text2(y0_2), color="r")

axes[0].plot(sol3.t, sol3.y[0], label=text(y0_3),color="pink")
axes[1].plot(sol3.t, sol3.y[1], label=text2(y0_3),color="pink")

axes[0].plot(sol4.t, sol4.y[0], label=text(y0_4),color="purple")
axes[1].plot(sol4.t, sol4.y[1], label=text2(y0_4),color="purple")

axes[0].plot(sol5.t, sol5.y[0], label=text(y0_5),color="green")
axes[1].plot(sol5.t, sol5.y[1], label=text2(y0_5),color="green")

axes[0].legend()
axes[1].legend()
#plt.text(0,0,f"y0(Blue)={y0}, y0(RED)={y0_2}")
#plt.tight_layout()
plt.show()

#1
x1 = x1_fun(sol)
y1 = y1_fun(sol)

x2 = x2_fun(sol)
y2 = y2_fun(sol)

#2
x1_2 = x1_fun(sol2)
y1_2 = y1_fun(sol2)

x2_2 = x2_fun(sol2)
y2_2 = y2_fun(sol2)

#3
x1_3 = x1_fun(sol3)
y1_3 = y1_fun(sol3)

x2_3 = x2_fun(sol3)
y2_3 = y2_fun(sol3)

#4
x1_4 = x1_fun(sol4)
y1_4 = y1_fun(sol4)

x2_4 = x2_fun(sol4)
y2_4 = y2_fun(sol4)

#5
x1_5 = x1_fun(sol5)
y1_5 = y1_fun(sol5)

x2_5 = x2_fun(sol5)
y2_5 = y2_fun(sol5)

length = min(len(x1),len(x1_2),len(x1_3),len(x1_4),len(x1_5))-1
print(f"lenx1: {len(x1)}, lenx2: {len(x1_2)}, length: {length}")

class DoublePendulum(Scene):
    def construct(self):
        times = ValueTracker(0)
        pendulum1 = always_redraw(lambda: Line(ORIGIN, [x1[int(times.get_value())], y1[int(times.get_value())], 0], color=BLUE))
        pendulum2 = always_redraw(lambda: Line(pendulum1.get_end(), [x2[int(times.get_value())], y2[int(times.get_value())], 0], color=BLUE))

        bob1 = always_redraw(lambda: Dot(pendulum1.get_end(), color=BLUE))
        bob2 = always_redraw(lambda: Dot(pendulum2.get_end(), color=BLUE))
        pendulum_group = VGroup(pendulum1, pendulum2, bob1, bob2)

#2
        pendulum1_2 = always_redraw(lambda: Line(ORIGIN, [x1_2[int(times.get_value())], y1_2[int(times.get_value())], 0], color=RED))
        pendulum2_2 = always_redraw(lambda: Line(pendulum1_2.get_end(), [x2_2[int(times.get_value())], y2_2[int(times.get_value())], 0], color=RED))

        bob1_2 = always_redraw(lambda: Dot(pendulum1_2.get_end(), color=RED))
        bob2_2 = always_redraw(lambda: Dot(pendulum2_2.get_end(), color=RED))
        pendulum_group_2 = VGroup(pendulum1_2, pendulum2_2, bob1_2, bob2_2)

#3
        pendulum1_3 = always_redraw(lambda: Line(ORIGIN, [x1_3[int(times.get_value())], y1_3[int(times.get_value())], 0], color=PINK))
        pendulum2_3 = always_redraw(lambda: Line(pendulum1_3.get_end(), [x2_3[int(times.get_value())], y2_3[int(times.get_value())], 0], color=PINK))

        bob1_3 = always_redraw(lambda: Dot(pendulum1_3.get_end(), color=PINK))
        bob2_3 = always_redraw(lambda: Dot(pendulum2_3.get_end(), color=PINK))
        pendulum_group_3 = VGroup(pendulum1_3, pendulum2_3, bob1_3, bob2_3)

#4
        pendulum1_4 = always_redraw(lambda: Line(ORIGIN, [x1_4[int(times.get_value())], y1_4[int(times.get_value())], 0], color=PURPLE))
        pendulum2_4 = always_redraw(lambda: Line(pendulum1_4.get_end(), [x2_4[int(times.get_value())], y2_4[int(times.get_value())], 0], color=PURPLE_C))

        bob1_4 = always_redraw(lambda: Dot(pendulum1_4.get_end(), color=PURPLE))
        bob2_4 = always_redraw(lambda: Dot(pendulum2_4.get_end(), color=PURPLE_C))
        pendulum_group_4 = VGroup(pendulum1_4, pendulum2_4, bob1_4, bob2_4)

#5
        pendulum1_5 = always_redraw(lambda: Line(ORIGIN, [x1_5[int(times.get_value())], y1_5[int(times.get_value())], 0], color=GREEN))
        pendulum2_5 = always_redraw(lambda: Line(pendulum1_5.get_end(), [x2_5[int(times.get_value())], y2_5[int(times.get_value())], 0], color=GREEN_C))

        bob1_5 = always_redraw(lambda: Dot(pendulum1_5.get_end(), color=GREEN))
        bob2_5 = always_redraw(lambda: Dot(pendulum2_5.get_end(), color=GREEN_C))
        pendulum_group_5 = VGroup(pendulum1_5, pendulum2_5, bob1_5, bob2_5)

        self.add(pendulum_group_5)
        self.add(pendulum_group_4)
        self.add(pendulum_group_3)
        self.add(pendulum_group_2)
        self.add(pendulum_group)
        self.play(times.animate.set_value(length), rate_func=linear, run_time=time)