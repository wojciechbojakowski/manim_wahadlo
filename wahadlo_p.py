from manim import *
import numpy as np

# Parametry wahadła
g = 9.81  #[m/s^2]
l = 2.0   #[m]
gamma = 0.05  # współczynnik tarcia
dt = 0.02  #[s]
t_max = 10  #[s]

# Warunki początkowe
theta0 = np.radians(30) 
omega0 = 0.0 
times = np.arange(0, t_max, dt)

# Listy na wyniki
theta = [theta0]
omega = [omega0]
energy = [(0.5 * (l**2) * omega0**2 + g * l * (1 - np.cos(theta0)))]

# Metoda Eulera do rozwiązania równań ruchu
for _ in times[1:]:
    omega_new = omega[-1] + (-gamma * omega[-1] - (g / l) * np.sin(theta[-1])) * dt
    theta_new = theta[-1] + omega_new * dt
    omega.append(omega_new)
    theta.append(theta_new)
    kinetic = 0.5 * (l**2) * omega_new**2
    potential = g * l * (1 - np.cos(theta_new))
    energy.append(kinetic + potential)

# Klasa animacji wahadła
class PendulumAnimation(Scene):
    def construct(self):
        # Tworzenie osi i linii wahadła
        pivot = np.array([0, 1, 0])
        rod = Line(pivot, pivot + l * np.array([np.sin(theta[0]), -np.cos(theta[0]), 0]), color=WHITE)
        end = Dot(rod.get_end(), color=RED, radius=0.1)
        
        # Wykresy kąta i energii
        theta_graph = Axes(
            x_range=[0, t_max, 2],
            y_range=[-np.radians(40), np.radians(40), np.radians(10)],
            axis_config={"color": WHITE}
        ).scale(0.5).to_corner(UL)
        theta_plot = theta_graph.plot_line_graph(times, theta, add_vertex_dots=False, line_color=BLUE)

        energy_graph = Axes(
            x_range=[0, t_max, 2],
            y_range=[0, max(energy), max(energy)/4],
            axis_config={"color": WHITE}
        ).scale(0.5).to_corner(UR)
        energy_plot = energy_graph.plot_line_graph(times, energy, add_vertex_dots=False, line_color=YELLOW)
        
        # Grupa do aktualizacji
        pendulum = VGroup(rod, end)
        self.add(pendulum, theta_graph, theta_plot, energy_graph, energy_plot)

        
        # def update_pendulum(mob, dt):
        #     index = int(TAU/4 * dt)
        #     print(f"dt: {dt}, index: {index}")
        #     new_x = l * np.sin(theta[index])
        #     new_y = -l * np.cos(theta[index])
        #     new_pos = pivot + np.array([new_x, new_y, 0])
        #     rod.put_start_and_end_on(pivot, new_pos)
        #     bob.move_to(new_pos)

        # pendulum.add_updater(update_pendulum)
        self.add(pendulum)
        self.wait(t_max)
