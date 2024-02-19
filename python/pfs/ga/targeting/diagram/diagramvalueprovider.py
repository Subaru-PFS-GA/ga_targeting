from ..photometry import Color, Magnitude
from ..diagram import ColorAxis, MagnitudeAxis

class DiagramValueProvider():
    """
    Mixin to access data values for plotting.
    """

    def get_magnitude(self, magnitude: Magnitude, DM=None, observed=False, mask=None):
        raise NotImplementedError()

    def get_color(self, color: Color, observed=False, mask=None):
        raise NotImplementedError()

    def get_diagram_values(self, axes, observed=False, mask=None):
        # Collect data
        x = []
        for ax in axes:
            if isinstance(ax, ColorAxis):
                x.append(self.get_color(ax.color, observed=observed, mask=mask))
            elif isinstance(ax, MagnitudeAxis):
                x.append(self.get_magnitude(ax.magnitude, observed=observed, mask=mask))
            elif isinstance(ax, Color):
                x.append(self.get_color(ax, observed=observed, mask=mask))
            elif isinstance(ax, Magnitude):
                x.append(self.get_magnitude(ax, observed=observed, mask=mask))
            else:
                raise NotImplementedError()
        return x