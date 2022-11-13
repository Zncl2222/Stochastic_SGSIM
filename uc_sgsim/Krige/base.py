class Kriging:
    def __init__(self, model):
        self.model = model
        self.hs = model.hs
        self.bw = model.bw
        self.a = model.a
        self.C0 = model.C0
