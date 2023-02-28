class Kriging:
    def __init__(self, model):
        self.__model = model
        self.__bandwidth_step = model.bandwidth_step
        self.__bandwidth = model.bandwidth
        self.__a = model.a
        self.__C0 = model.C0

    @property
    def model(self):
        return self.__model

    @property
    def bandwidth_step(self):
        return self.__bandwidth_step

    @property
    def bandwidth(self):
        return self.__bandwidth

    @property
    def a(self):
        return self.__a

    @property
    def C0(self):
        return self.__C0
