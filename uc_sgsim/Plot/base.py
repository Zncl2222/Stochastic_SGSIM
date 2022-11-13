class Plot_Base:
    def __init__(self, model, RandomField, figsize=(10, 8)):
        self.model = model
        self.RandomField = RandomField
        self.figsize = figsize
        self.model_name = model.model_name
        self.hs = model.hs
        self.bw = model.bw
        self.a = model.a
        self.C0 = model.C0
        self.size = len(RandomField)
        self.nR = len(RandomField[0])
