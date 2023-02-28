class PlotBase:
    def __init__(self, model, random_field, figsize=(10, 8)):
        self.__model = model
        self.__random_field = random_field
        self.__figsize = figsize
        self.__model_name = model.model_name
        self.__bandwidth_step = model.bandwidth_step
        self.__bandwidth = model.bandwidth
        self.__a = model.a
        self.__C0 = model.C0
        self.__size = len(random_field)
        self.__realization_number = len(random_field[:, 0])

    @property
    def model(self):
        return self.__model

    @property
    def random_field(self):
        return self.__random_field

    @property
    def figsize(self):
        return self.__figsize

    @property
    def model_name(self):
        return self.__model_name

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

    @property
    def size(self):
        return self.__size

    @property
    def realization_number(self):
        return self.__realization_number
