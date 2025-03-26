from manim import *

class PendulumWithGraph(Scene):
    def construct(self):
        times = ValueTracker(0)
        theta_max = PI/6
        l = 3
        w = np.sqrt(15/3)
        T = 2*PI/w

        refpoint = 3*UP

        theta = DecimalNumber().set_color(BLACK).move_to(15*RIGHT)
        theta.add_updater(lambda m: m.set_value(theta_max * np.sin(w*times.get_value())))
        self.add(theta)

        def get_line1(x,y):
            line = Line(start=ORIGIN+refpoint, end=x*RIGHT+y*UP+refpoint, color=BLUE)
            return line

        line = always_redraw(lambda: get_line1(l*np.sin(theta.get_value()), -l*np.cos(theta.get_value())))
        self.add(line)

        verticalline = DashedLine(start=line.get_start(), end=line.get_start()+3*DOWN+0.001*LEFT) # Add a small angle
        self.add(verticalline)

        def angle_arc(theta):
            if theta > 0:
                angle = Angle(line, verticalline, quadrant=(1, 1), other_angle=True, color=YELLOW, fill_opacity=0)
            else:
                angle = Angle(line, verticalline, quadrant=(1, 1), other_angle=False, color=YELLOW, fill_opacity=0)
            return angle

        angle = always_redraw(lambda: angle_arc(theta.get_value()))
        self.add(angle)

        arctext = MathTex(r"\theta").scale(0.5).add_updater(lambda m: m.next_to(angle, DOWN))
        self.add(arctext)

        def get_ball(x, y):
            dot = Dot(fill_color=BLUE, fill_opacity=1).move_to(x*RIGHT+y*UP+refpoint).scale(l)
            return dot

        ball = always_redraw(lambda: get_ball(l*np.sin(theta.get_value()), -l*np.cos(theta.get_value())))

        traced_path = TracedPath(ball.get_center, stroke_color=RED, stroke_width=2)
        self.add(ball, traced_path)

        dot_red = Dot(color=RED)
        dot_red.move_to(line.get_start()+3*DOWN)

        def update_dot_red(dot):
            if theta.get_value() > 0:
                dot.move_to(line.get_start()+3*DOWN)
            else:
                dot.move_to(line.get_start()+3*UP)

        dot_red.add_updater(update_dot_red)
        self.add(dot_red)

        # Create a graph of theta vs time
        axes = Axes(
            x_range=[0, 10, 1],
            y_range=[-1, 1, 0.2],
            axis_config={"color": WHITE},
            x_axis_config={"color": WHITE},
            y_axis_config={"color": WHITE}
        ).scale(0.6).to_edge(DOWN)

        graph = always_redraw(lambda: axes.plot(lambda t: theta_max * np.sin(w*t), color=GREEN))

        # Add a red dot on the graph
        dot_on_graph = always_redraw(lambda: Dot().move_to(axes.c2p(times.get_value() % T, theta_max * np.sin(w*times.get_value()))))

        self.add(axes, graph, dot_on_graph)
        
        self.play(times.animate.set_value(3*T), rate_func=linear, run_time=3*T)

# Render the scene
# if __name__ == "__main__":
#     script = f"{Path(__file__).resolve()}"
#     os.system(f"manim {script} PendulumWithGraph -pqh")