from manim import *
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import tempfile
import os

# Parametry fizyczne
g = 9.81
L = 1.0
m = 1.0
t_max = 20
t_eval = np.linspace(0, t_max, 1000)

# Równania ruchu podwójnego wahadła
def deriv(t, y):
    θ1, θ2, ω1, ω2 = y
    Δ = θ2 - θ1

    denom = (m * (2 - np.cos(2*Δ)))
    dθ1 = ω1
    dθ2 = ω2

    dω1 = (m * g * np.sin(θ2) * np.cos(Δ) -
           m * np.sin(Δ) * (L * ω2**2 + L * ω1**2 * np.cos(Δ)) -
           (2 * m) * g * np.sin(θ1)) / (L * denom)

    dω2 = ((2 * np.sin(Δ)) *
           (L * ω1**2 * m + g * m * np.cos(θ1) +
            L * ω2**2 * m * np.cos(Δ))) / (L * denom)

    return [dθ1, dθ2, dω1, dω2]

class PodwojneWahadloChaosDemo(Scene):
    def construct(self):
        base_theta1 = np.pi / 4
        deltas = np.linspace(-0.001, 0.001, 5)
        y0s = [[base_theta1 + d, np.pi / 4, 0, 0] for d in deltas]
        sols = [solve_ivp(deriv, [0, t_max], y0, t_eval=t_eval) for y0 in y0s]
        colors = [RED, ORANGE, GREEN, BLUE, PURPLE]

        global times
        times = ValueTracker(0)

        # Do wykresu różnicy
        theta1_1 = sols[0].y[0]
        theta1_2 = sols[1].y[0]
        delta_theta1 = np.abs(theta1_1 - theta1_2)

        # Tworzenie elementów animacji
        for i, (sol, color) in enumerate(zip(sols, colors)):
            theta1, theta2 = sol.y[0], sol.y[1]

            # Użyj np.array dla pozycji
            x1 = L * np.sin(theta1)
            y1 = -L * np.cos(theta1)
            x2 = x1 + L * np.sin(theta2)
            y2 = y1 - L * np.cos(theta2)

            arm1 = always_redraw(lambda i=i, color=color: Line(
                start=ORIGIN,
                end= np.array([L * np.sin(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]),
                               -L * np.cos(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]), 0]),
                color=color, stroke_width=1.5
            ))

            bob1 = always_redraw(lambda i=i, color=color: Dot(
                point= np.array([L * np.sin(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]),
                                -L * np.cos(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]), 0]),
                radius=0.03, color=color
            ))

            arm2 = always_redraw(lambda i=i, color=color: Line(
                start= np.array([L * np.sin(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]),
                                 -L * np.cos(sols[i].y[0][min(int(times.get_value()), len(sols[i].y[0])-1)]), 0]),
                end= np.array([L * np.sin(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]),
                                 -L * np.cos(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]), 0]),
                color=color, stroke_width=1.5
            ))

            bob2 = always_redraw(lambda i=i, color=color: Dot(
                point= np.array([L * np.sin(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]),
                                 -L * np.cos(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]), 0]),
                radius=0.03, color=color
            ))

            # Trasa drugiej kulki
            trace = TracedPath(
                lambda i=i, color=color: np.array([L * np.sin(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]),
                                                   -L * np.cos(sols[i].y[1][min(int(times.get_value()), len(sols[i].y[1])-1)]), 0]),
                stroke_color=color,
                stroke_width=2,
                dissipating_time=0
            )

            self.add(trace, arm1, arm2, bob1, bob2)

        # Dodaj legendę
        legenda = VGroup(*[
            VGroup(
                Dot(color=c, radius=0.05),
                Text(f"Δθ₁ = {round(d, 5)}", font_size=24)
            ).arrange(RIGHT, buff=0.2)
            for c, d in zip(colors, deltas)
        ]).arrange(DOWN).to_corner(UR).shift(LEFT * 0.5)

        self.add(legenda)

        # Animacja czasu
        self.play(times.animate.set_value(len(t_eval) - 1), run_time=t_max, rate_func=linear)
        self.wait()

        # WYKRES delta_theta1
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(t_eval, delta_theta1, color="black")
        ax.set_title("Różnica między θ₁(0) i θ₁′(t)")
        ax.set_xlabel("Czas [s]")
        ax.set_ylabel("|Δθ₁(t)|")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmpfile:
            plt.savefig(tmpfile.name)
            wykres_path = tmpfile.name

        wykres = ImageMobject(wykres_path).scale(0.8).to_edge(DOWN)
        self.play(FadeIn(wykres))
        self.wait(3)

        os.remove(wykres_path)
