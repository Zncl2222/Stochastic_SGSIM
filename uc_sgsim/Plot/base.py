class PlotBase:
    def __init__(self, model, random_field, figsize=(10, 8)):
        self.model = model
        self.random_field = random_field
        self.figsize = figsize
        self.model_name = model.model_name
        self.bandwidth_step = model.bandwidth_step
        self.bandwidth = model.bandwidth
        self.a = model.a
        self.C0 = model.C0
        self.size = len(random_field)
        self.realization_number = len(random_field[:, 0])
