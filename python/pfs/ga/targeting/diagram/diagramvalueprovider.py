from ..photometry import Color, Magnitude
from ..diagram import Diagram, ColorAxis, MagnitudeAxis, RaDecAxis

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
    
    def has_coords(self, ra=None, dec=None):
        raise NotImplementedError()
    
    def get_coords(self, ra=None, dec=None, mask=None, ctype=None):
        raise NotImplementedError()
    
    def has_diagram_values(self, axes, observed=False):
        if isinstance(axes, Diagram):
            axes = axes.axes

        if len(axes) == 2 and isinstance(axes[0], RaDecAxis) and isinstance(axes[1], RaDecAxis):
            return self.has_coords(ra=axes[0].coord, dec=axes[1].coord)
        else:
            # Check if all magnitudes are available (for plotting, for example)
            found = True
            for ax in axes:
                if isinstance(ax, (Color, ColorAxis)):
                    found &= self.has_color(ax.color, observed=observed)
                elif isinstance(ax, (Magnitude, MagnitudeAxis)):
                    found &= self.has_magnitude(ax.magnitude, observed=observed)
                else:
                    raise NotImplementedError()
            return found

    def get_diagram_values(self, axes, observed=False, mask=None):
        if isinstance(axes, Diagram):
            axes = axes.axes

        if len(axes) == 2 and isinstance(axes[0], RaDecAxis) and isinstance(axes[1], RaDecAxis):
            return self.get_coords(ra=axes[0].coord, dec=axes[1].coord, mask=mask)
        else:
            # Collect magnitudes and colors
            x = []
            for ax in axes:
                if isinstance(ax, ColorAxis):
                    x.append(self.get_color(ax.color, observed=observed, mask=mask))
                elif isinstance(ax, Color):
                    x.append(self.get_color(ax, observed=observed, mask=mask))
                elif isinstance(ax, MagnitudeAxis):
                    x.append(self.get_magnitude(ax.magnitude, observed=observed, mask=mask))
                elif isinstance(ax, Magnitude):
                    x.append(self.get_magnitude(ax, observed=observed, mask=mask))
                else:
                    raise NotImplementedError()
            
        return x