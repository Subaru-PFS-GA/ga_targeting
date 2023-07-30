from .projection import Projection

class TelescopeProjection(Projection):
    def __init__(self, pointing=None, orig=None):
        super(TelescopeProjection, self).__init__(pointing=pointing, orig=orig)

        if not isinstance(orig, TelescopeProjection):
            pass
        else:
            pass
