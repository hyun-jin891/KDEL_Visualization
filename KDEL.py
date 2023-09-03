from random import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import sys


ER_MEMBRANE = 19
GOLGI_MEMBRANE = 80
CHAPERONE_SPEED = 0.7
UNFOLED_PROTEIN_SPEED = 0.3
OTHER_PROTEIN_SPEED = 0.7
MEMBRANE_SPEED = 0.5
UNFOLDED_PROTEIN_GENERATION = 0.3
KDEL_RECEPTOR_GENERATION = 0.04
OTHER_RECEPTOR_GENERATION = 0.1
KDEL_RECEPTOR_AFFINITY = 0.3


class Resident_Soluble_Protein:
    def __init__(self, x, y):
        self._x = x
        self._y = y
    
    def setterXY(self, x, y):
        self._x = x
        self._y = y
    
    def getterX(self):
        return self._x
    
    def getterY(self):
        return self._y
    
    def move(self):
        while True:
            change_X = sample([-1, 1], 1)[0]
            change_Y = sample([-1, 1], 1)[0]
            if 2 <= self._y + change_Y <= 97:
                if (1 <= self._x + change_X < ER_MEMBRANE) or (GOLGI_MEMBRANE < self._x + change_X <= 98):
                    self._x += change_X
                    self._y += change_Y
                    return

class ER_Chaperone(Resident_Soluble_Protein):
    pass
    

class Unfolded_Protein(Resident_Soluble_Protein):
    pass


class Other_Protein(Resident_Soluble_Protein):
    pass
    


class Receptor:
    def __init__(self, x, y):
        self._x = x
        self._y = y
    
    def setterXY(self, x, y):
        self._x = x
        self._y = y
    
    def getterX(self):
        return self._x
    
    def getterY(self):
        return self._y
        
    def move(self):
        while True:
            change = sample([-1, 1], 1)[0]
            if 2 <= self._y + change <= 97:
                self._y += change
                return

class KDEL_Receptor(Receptor):
    pass
            
class Other_Receptor(Receptor):
    pass            
        

class Cell_World:
    def __init__(self, n1, n2, n3, n4, n5):
        self.__map = [[0 for i in range(100)] for j in range(100)]
        self.er_chaperone = []
        self.unfolded_protein = []
        self.other_protein = []
        self.kdel_receptor = []
        self.other_receptor = []
        
        for i in range(n1):
            weight = [0.8, 0.2]
            while True:
                random_x = randint(1, ER_MEMBRANE - 1)
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 1
                    self.er_chaperone.append(ER_Chaperone(random_x, random_y))
                    break
        
        for i in range(n2):
            weight = [0.8, 0.2]
            while True:
                random_x = choices([randint(1, ER_MEMBRANE - 1), randint(GOLGI_MEMBRANE + 1, 98)], weight)[0]
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 2
                    self.unfolded_protein.append(Unfolded_Protein(random_x, random_y))
                    break
        
        for i in range(n3):
            while True:
                random_x = sample([randint(1, ER_MEMBRANE - 1), randint(GOLGI_MEMBRANE + 1, 98)], 1)[0]
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 3
                    self.other_protein.append(Other_Protein(random_x, random_y))
                    break
        
        for i in range(n4):
            while True:
                random_x = GOLGI_MEMBRANE
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 4
                    self.kdel_receptor.append(KDEL_Receptor(random_x, random_y))
                    break
        
        for i in range(n5):
            while True:
                random_x = ER_MEMBRANE
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 5
                    self.other_receptor.append(Other_Receptor(random_x, random_y))
                    break
                    
    def getterMap(self):
        return self.__map
    
    def numUnfolded_Protein(self):
        num1 = 0
        num2 = 0
        
        for i in range(len(self.unfolded_protein)):
            if self.unfolded_protein[i] == None:
                continue
            else:
                if 1 <= self.unfolded_protein[i].getterX() < ER_MEMBRANE:
                    num1 += 1
                elif GOLGI_MEMBRANE < self.unfolded_protein[i].getterX() <= 98:
                    num2 += 1
        
        return num1, num2
    

    def runWorld(self):
        self.unfolded_protein = list(filter(None, self.unfolded_protein))
        p = uniform(0.0, 1.0)
        
        if p <= UNFOLDED_PROTEIN_GENERATION:
            while True:
                random_x = randint(0, ER_MEMBRANE - 1)
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 2
                    self.unfolded_protein.append(Unfolded_Protein(random_x, random_y))
                    break
        
        if p <= KDEL_RECEPTOR_GENERATION:
            while True:
                random_x = ER_MEMBRANE
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 4
                    self.kdel_receptor.append(KDEL_Receptor(random_x, random_y))
                    break
        
        
        if p <= OTHER_RECEPTOR_GENERATION:
            while True:
                random_x = ER_MEMBRANE
                random_y = randint(2, 97)
                
                if self.__map[random_y][random_x] == 0:
                    self.__map[random_y][random_x] = 5
                    self.other_receptor.append(Other_Receptor(random_x, random_y))
                    break
        
        
        for i in range(len(self.er_chaperone)):
            current_Chaperone = self.er_chaperone[i]
            current_x = current_Chaperone.getterX()
            current_y = current_Chaperone.getterY()
            
            if self.__map[current_y][current_x] == 6:
                flag = False

                if self.__map[current_y][current_x + 1] != 0:
                    if self.__map[current_y][current_x + 1] == 7:
                        current_Chaperone.setterXY(current_x + 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_Chaperone.setterXY(current_x + 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_Chaperone.getterX() == GOLGI_MEMBRANE:
                    while True:
                        current_Chaperone.move()
                        if self.__map[current_Chaperone.getterY()][current_Chaperone.getterX() + 4] == 0:
                           self.__map[current_y][current_x] = 0
                           self.__map[current_Chaperone.getterY()][current_Chaperone.getterX() + 4] = 1
                           break
                        else:
                            current_Chaperone.setterXY(GOLGI_MEMBRANE, current_y)
                    continue 
                
                self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 6
                self.__map[current_y][current_x] = 0
                continue
            
            elif self.__map[current_y][current_x] == 7:
                flag = False
                if self.__map[current_y][current_x - 1] != 0:
                    if self.__map[current_y][current_x - 1] == 6:
                        current_Chaperone.setterXY(current_x - 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_Chaperone.setterXY(current_x - 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_Chaperone.getterX() == ER_MEMBRANE:
                    while True:
                        current_Chaperone.move()
                        if self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] == 0:
                            self.__map[current_y][current_x] = 0
                            self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 1
                            break
                        else:
                            current_Chaperone.setterXY(ER_MEMBRANE, current_y)
                    continue
                
                self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 7
                self.__map[current_y][current_x] = 0
                continue
            
            factor = uniform(0.0, 1.0)
            
            if factor <= CHAPERONE_SPEED:
                current_Chaperone.move()
                if self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] != 0:
                    current_Chaperone.setterXY(current_x, current_y)
                else:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 1
            
            direction = [(-1, 0), (0, -1), (1, 0), (0, 1)]
            for j in range(4):
                search_locationX = current_Chaperone.getterX() + direction[j][0]
                search_locationY = current_Chaperone.getterY() + direction[j][1]
                if search_locationX == ER_MEMBRANE or search_locationX == GOLGI_MEMBRANE:
                    continue
                elif search_locationY < 0 or search_locationY > 99:
                    continue
                
                if self.__map[search_locationY][search_locationX] == 2:
                    self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 0
                    current_Chaperone.setterXY(search_locationX, search_locationY)
                    self.__map[current_Chaperone.getterY()][current_Chaperone.getterX()] = 1
                    break
        
        for i in range(len(self.unfolded_protein)):
            current_unfolded_protein = self.unfolded_protein[i]
            current_x = current_unfolded_protein.getterX()
            current_y = current_unfolded_protein.getterY()
            
            if self.__map[current_y][current_x] == 1:
                self.unfolded_protein[i] = None
                del current_unfolded_protein
                continue
            
            if self.__map[current_y][current_x] == 6:
                flag = False
                if self.__map[current_y][current_x + 1] != 0:
                    if self.__map[current_y][current_x + 1] == 7:
                        current_unfolded_protein.setterXY(current_x + 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_unfolded_protein.setterXY(current_x + 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_unfolded_protein.getterX() == GOLGI_MEMBRANE:
                    while True:
                        current_unfolded_protein.move()
                        if self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] == 0:
                           self.__map[current_y][current_x] = 0
                           self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] = 2
                           break
                        else:
                            current_unfolded_protein.setterXY(GOLGI_MEMBRANE, current_y)
                    continue 
                
                self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] = 6
                self.__map[current_y][current_x] = 0
                continue
            
            elif self.__map[current_y][current_x] == 7:
                flag = False
                if self.__map[current_y][current_x - 1] != 0:
                    if self.__map[current_y][current_x - 1] == 6:
                        current_unfolded_protein.setterXY(current_x - 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_unfolded_protein.setterXY(current_x - 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_unfolded_protein.getterX() == ER_MEMBRANE:
                    while True:
                        current_unfolded_protein.move()
                        if self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] == 0:
                            self.__map[current_y][current_x] = 0
                            self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] = 2
                            break
                        else:
                            current_unfolded_protein.setterXY(ER_MEMBRANE, current_y)
                    continue
                
                self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] = 7
                self.__map[current_y][current_x] = 0
                continue
            
            factor = uniform(0.0, 1.0)
            
            if factor <= UNFOLED_PROTEIN_SPEED:
                current_unfolded_protein.move()
                if self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] != 0:
                    current_unfolded_protein.setterXY(current_x, current_y)
                else:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_unfolded_protein.getterY()][current_unfolded_protein.getterX()] = 2
        
        
        for i in range(len(self.other_protein)):
            current_other_protein = self.other_protein[i]
            current_x = current_other_protein.getterX()
            current_y = current_other_protein.getterY()
            
            if self.__map[current_y][current_x] == 6:
                flag = False
                if self.__map[current_y][current_x + 1] != 0:
                    if self.__map[current_y][current_x + 1] == 7:
                        current_other_protein.setterXY(current_x + 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_other_protein.setterXY(current_x + 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_other_protein.getterX() == GOLGI_MEMBRANE:
                    while True:
                        current_other_protein.move()
                        if self.__map[current_other_protein.getterY()][current_other_protein.getterX()] == 0:
                           self.__map[current_y][current_x] = 0
                           self.__map[current_other_protein.getterY()][current_other_protein.getterX()] = 3
                           break
                        else:
                            current_other_protein.setterXY(GOLGI_MEMBRANE, current_y)
                    continue 
                
                self.__map[current_other_protein.getterY()][current_other_protein.getterX()] = 6
                self.__map[current_y][current_x] = 0
                continue
            
            elif self.__map[current_y][current_x] == 7:
                flag = False
                if self.__map[current_y][current_x - 1] != 0:
                    if self.__map[current_y][current_x - 1] == 6:
                        current_other_protein.setterXY(current_x - 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_other_protein.setterXY(current_x - 1, current_y)
                #self.__map[current_y][current_x] = 0
                if current_other_protein.getterX() == ER_MEMBRANE:
                    while True:
                        current_other_protein.move()
                        if self.__map[current_other_protein.getterY()][current_other_protein.getterX() + 4] == 0:
                            self.__map[current_y][current_x] = 0
                            self.__map[current_other_protein.getterY()][current_other_protein.getterX() + 4] = 3
                            break
                        else:
                            current_other_protein.setterXY(ER_MEMBRANE, current_y)
                    continue
                
                self.__map[current_other_protein.getterY()][current_other_protein.getterX()] = 7
                self.__map[current_y][current_x] = 0
                continue
            
            factor = uniform(0.0, 1.0)
            
            if factor <= OTHER_PROTEIN_SPEED:
                current_other_protein.move()
                if self.__map[current_other_protein.getterY()][current_other_protein.getterX()] != 0:
                    current_other_protein.setterXY(current_x, current_y)
                else:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_other_protein.getterY()][current_other_protein.getterX()] = 3
        
        for i in range(len(self.kdel_receptor)):
            current_kdel_receptor = self.kdel_receptor[i]
            current_x = current_kdel_receptor.getterX()
            current_y = current_kdel_receptor.getterY()
            
            if self.__map[current_y][current_x] == 6:
                flag = False

                if self.__map[current_y][current_x + 1] != 0:
                    if self.__map[current_y][current_x + 1] == 7:
                        current_kdel_receptor.setterXY(current_x + 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_kdel_receptor.setterXY(current_x + 1, current_y)
 
                #self.__map[current_y][current_x] = 0
                if current_kdel_receptor.getterX() == GOLGI_MEMBRANE:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_y][current_kdel_receptor.getterX()] = 4
                    continue
                
                self.__map[current_kdel_receptor.getterY()][current_kdel_receptor.getterX()] = 6
                self.__map[current_y][current_x] = 0
                continue
            
            elif self.__map[current_y][current_x] == 7:
                flag = False

                if self.__map[current_y][current_x - 1] != 0:
                    if self.__map[current_y][current_x - 1] == 6:
                        current_kdel_receptor.setterXY(current_x-3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_kdel_receptor.setterXY(current_x - 1, current_y)
                #self.__map[current_y][current_x] = 0
                if current_kdel_receptor.getterX() == ER_MEMBRANE:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_y][current_kdel_receptor.getterX()] = 4
                    continue
                
                self.__map[current_kdel_receptor.getterY()][current_kdel_receptor.getterX()] = 7
                self.__map[current_y][current_x] = 0
                continue
                        
            factor = uniform(0.0, 1.0)
            
            if factor <= MEMBRANE_SPEED:
                current_kdel_receptor.move()
                if self.__map[current_kdel_receptor.getterY()][current_kdel_receptor.getterX()] != 0:
                    current_kdel_receptor.setterXY(current_x, current_y)
                else:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_kdel_receptor.getterY()][current_kdel_receptor.getterX()] = 4
                       
            if current_kdel_receptor.getterX() == GOLGI_MEMBRANE:
                p = uniform(0.0, 1.0)
                if self.__map[current_kdel_receptor.getterY()][GOLGI_MEMBRANE + 1] == 1 and p <= KDEL_RECEPTOR_AFFINITY:
                    self.__map[current_kdel_receptor.getterY()][GOLGI_MEMBRANE] = 7
                    self.__map[current_kdel_receptor.getterY()][GOLGI_MEMBRANE + 1] = 7
                    
                    direction = [(0, -1), (0, 1), (1, -1), (1, 1)]
                    
                    for j in range(len(direction)):
                        search_locationX = current_kdel_receptor.getterX() + direction[j][0]
                        search_locationY = current_kdel_receptor.getterY() + direction[j][1]
                        if self.__map[search_locationY][search_locationX] != 0:
                            self.__map[search_locationY][search_locationX] = 7       
            
        
        for i in range(len(self.other_receptor)):
            current_other_receptor = self.other_receptor[i]
            current_x = current_other_receptor.getterX()
            current_y = current_other_receptor.getterY()
            
            if self.__map[current_y][current_x] == 6:
                flag = False

                if self.__map[current_y][current_x + 1] != 0:
                    if self.__map[current_y][current_x + 1] == 7:
                        current_other_receptor.setterXY(current_x + 3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_other_receptor.setterXY(current_x + 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_other_receptor.getterX() == GOLGI_MEMBRANE:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_y][current_other_receptor.getterX()] = 5
                    continue
                
                self.__map[current_other_receptor.getterY()][current_other_receptor.getterX()] = 6
                self.__map[current_y][current_x] = 0
                continue
            
            elif self.__map[current_y][current_x] == 7:
                flag = False

                if self.__map[current_y][current_x - 1] != 0:
                    if self.__map[current_y][current_x - 1] == 6:
                        current_other_receptor.setterXY(current_x-3, current_y + 2)
                        flag = True
                    else:
                        continue
                if flag == False:
                    current_other_receptor.setterXY(current_x - 1, current_y)

                #self.__map[current_y][current_x] = 0
                if current_other_receptor.getterX() == ER_MEMBRANE:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_y][current_other_receptor.getterX()] = 5
                    continue
                
                self.__map[current_other_receptor.getterY()][current_other_receptor.getterX()] = 7
                self.__map[current_y][current_x] = 0
                continue
                        
            factor = uniform(0.0, 1.0)
            
            if factor <= MEMBRANE_SPEED:
                current_other_receptor.move()
                if self.__map[current_other_receptor.getterY()][current_other_receptor.getterX()] != 0:
                    current_other_receptor.setterXY(current_x, current_y)
                else:
                    self.__map[current_y][current_x] = 0
                    self.__map[current_other_receptor.getterY()][current_other_receptor.getterX()] = 5
                       
            if current_other_receptor.getterX() == ER_MEMBRANE:
                if self.__map[current_other_receptor.getterY()][ER_MEMBRANE - 1] == 3:
                    self.__map[current_other_receptor.getterY()][ER_MEMBRANE] = 6
                    self.__map[current_other_receptor.getterY()][ER_MEMBRANE - 1] = 6
                    
                    direction = [(0, -1), (0, 1), (-1, -1), (-1, 1)]
                    
                    for j in range(len(direction)):
                        search_locationX = current_other_receptor.getterX() + direction[j][0]
                        search_locationY = current_other_receptor.getterY() + direction[j][1]
                        if self.__map[search_locationY][search_locationX] != 0:
                            self.__map[search_locationY][search_locationX] = 6    
               
def onClick(event):
    global pause
    pause ^= True
    anim.event_source.stop()

def update(i):
    print(miniCellWorld.numUnfolded_Protein())
    miniCellWorld.runWorld()
    matrix = np.matrix(miniCellWorld.getterMap())
    ax.cla()
    im = ax.matshow(matrix, cmap = cmap, norm = norm)
    plt.vlines(ER_MEMBRANE, 0, 99, color='black', linestyle='solid', linewidth=0.5)
    plt.vlines(GOLGI_MEMBRANE, 0, 99, color='black', linestyle='solid', linewidth=0.5)
    fig.tight_layout()
    
    return [im]
               
                
n1, n2, n3, n4, n5 = map(int, input().split())               


pause = False

miniCellWorld = Cell_World(n1, n2, n3, n4, n5)

matrix = np.matrix(miniCellWorld.getterMap())

fig, ax = plt.subplots()

y = np.arange(0, 99)

plt.vlines(ER_MEMBRANE, 0, 99, color='black', linestyle='solid', linewidth=0.5)
plt.vlines(GOLGI_MEMBRANE, 0, 99, color='black', linestyle='solid', linewidth=0.5)

fig.set_figheight(5)
fig.set_figwidth(5)

fig.canvas.mpl_connect('button_press_event', onClick)
fig.tight_layout()

cmap = colors.ListedColormap(['white', 'darkgreen', 'royalblue', 'black', 'aquamarine', 'darkgray', 'red', 'blue'])
norm = colors.BoundaryNorm(np.arange(len(matrix) + 1) - 0.5, len(matrix))
im = ax.matshow(matrix, cmap = cmap, norm = norm)

anim = animation.FuncAnimation(fig, update)

plt.show()
               
            
    
                
        






# %%
