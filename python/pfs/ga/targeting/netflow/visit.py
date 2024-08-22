class Visit():
    def __init__(self, visit_idx, pointing_idx, pointing, visit_cost=None):
        """
        Initializes a visit object.
        
        Variables:
        ----------
        visit_idx : int
            Unique zero-based index of the visit within the list of visits.
        pointing_idx : int
            Unique zero-based index of the pointing within the list of pointings.
            A pointing can have multiple visits but a visit can belong to only one pointing.
        pointing : Pointing
            Pointing object that the visit belongs to.
        visit_cost : float
            Cost of the visit. Default value is None.
        """
        self.__visit_idx = visit_idx
        self.__pointing_idx = pointing_idx
        self.__pointing = pointing
        self.__visit_cost = visit_cost

    def __get_visit_idx(self):
        return self.__visit_idx
    
    visit_idx = property(__get_visit_idx)

    def __get_pointing_idx(self):
        return self.__pointing_idx
    
    pointing_idx = property(__get_pointing_idx)

    def __get_pointing(self):
        return self.__pointing
    
    pointing = property(__get_pointing)

    def __get_visit_cost(self):
        return self.__visit_cost
    
    def __set_visit_cost(self, value):
        self.__visit_cost = value

    visit_cost = property(__get_visit_cost, __set_visit_cost)