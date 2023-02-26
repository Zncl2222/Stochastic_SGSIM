class Kriging:
    def __init__(self, model):
        self.model = model
        self.bandwidth_step = model.bandwidth_step
        self.bandwidth = model.bandwidth
        self.a = model.a
        self.C0 = model.C0
