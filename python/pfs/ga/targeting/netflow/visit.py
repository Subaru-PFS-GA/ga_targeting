class Visit():
    def __init__(self, visit_idx, pointing_idx, pointing, visit_cost=None):
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