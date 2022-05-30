import numpy as np


        
class GUI_base:
    
    def __init__(self):

        self.a=1


    def thread(func,*args):
        
        t=threading.Thread(target=func,args=args)
        t.setDaemon(True)
        t.start()
