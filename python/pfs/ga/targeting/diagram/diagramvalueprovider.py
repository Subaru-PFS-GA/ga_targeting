from ..photometry import Color, Magnitude
from ..diagram import ColorAxis, MagnitudeAxis

class DiagramValueProvider():
    """
    Mixin to access data values for plotting.
    """

    def has_magnitude(self, magnitude: Magnitude, observed=False, dered=False):
        raise NotImplementedError()

    def get_magnitude(self, magnitude: Magnitude, DM=None, observed=False, mask=None):
        raise NotImplementedError()
    
    def has_color(self, color: Color, observed=False, dered=False):
        return self.has_magnitude(color.magnitudes[0], observed=observed, dered=dered) and \
               self.has_magnitude(color.magnitudes[1], observed=observed, dered=dered)

    def get_color(self, color: Color, observed=False, mask=None):
        raise NotImplementedError()
    
    def has_diagram_values(self, axes, observed=False):
        # Check if all magnitudes are available (for plotting, for example)
        for ax in axes:
            if isinstance(ax, ColorAxis):
                if not self.has_color(ax.color, observed=observed):
                    return False
            elif isinstance(ax, MagnitudeAxis):
                if not self.get_magnitude(ax.magnitude, observed=observed):
                    return False
            elif isinstance(ax, Color):
                if not self.has_color(ax, observed=observed):
                    return False
            elif isinstance(ax, Magnitude):
                if not self.get_magnitude(ax, observed=observed):
                    return False
            else:
                raise NotImplementedError()
        return True

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