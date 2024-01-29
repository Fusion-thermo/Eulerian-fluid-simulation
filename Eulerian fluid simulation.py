#tutoriel : https://www.youtube.com/watch?v=iKAVRgIrUOU
from Classes import *
from tkinter import *
#import tkinter.filedialog
from statistics import median, mean

#A faire

#vérifier s'il faut lancer plusieurs d'affilé une simulation pour que cela avance
#

height=500
width = 1.5 * height

# h = grid_spacing = taille d'un carré
scene=Scene()
#modele="hydrostatique"
modele="velocity"

def setObstacleCircle(x,y,rayon):
      #x, y et rayon sont donnés en fraction de numX, numY et numY respectivement pour que la taille du maillage n'influence pas
      #la taille et l'emplacement de l'obstacle
      x=int(x*scene.fluid.numX)
      y=int(y*scene.fluid.numY)
      rayon=int(rayon*scene.fluid.numY)

      for i in range(x-rayon,x+rayon+1):
            for j in range(y-rayon,y+rayon+1):
                  if (i-x)**2 + (j-y)**2 <= rayon**2:
                        scene.fluid.grid[i,j].s = 0


def setup(modele):
      scene.dt = 1 / 60
      scene.numIters=40
      scene.overRelaxation = 1.9

      numX = 100
      h = floor(width/numX)
      numY = floor(height / h)

      density = 1000 # kg/m3

      scene.fluid=Fluid(numX,numY,h,density,scene.overRelaxation)
      
      if modele=="hydrostatique":
      #     scene.showPressure = True
      #     scene.showSmoke = False
      #     scene.showStreamlines = False
      #     scene.showVelocities = False
            pass
      elif modele=="velocity":
            #définit la vitesse d'entrée du fluide
            inletVel=2
            for j in range(1,scene.fluid.numY-1):
            #for j in range(int(scene.fluid.numY/2 - scene.fluid.numY/8),int(scene.fluid.numY/2 + scene.fluid.numY/8)):
                  scene.fluid.grid[1,j].u = inletVel
            
            #met m à 0 sur une zone centre sur le milieu
            #je suppose que m n'existe que pour la première colonne
            pipeH = 0.1 * scene.fluid.numY
            minJ = floor(0.5*scene.fluid.numY - 0.5*pipeH)
            maxJ = floor(0.5*scene.fluid.numY + 0.5*pipeH)
            for j in range(minJ,maxJ):
                  scene.fluid.grid[1,j].m = 0
            
            setObstacleCircle(0.2, 0.5, 0.1)

            scene.gravity = 0.0
            # scene.showPressure = False
            # scene.showSmoke = True
            # scene.showStreamlines = False
            # scene.showVelocities = False




def simulate():
      scene.fluid.simulate(scene.dt,scene.gravity,scene.numIters)

def update(modele):
      simulate()
      #draw(modele)

def rgb_to_hex(rgb):
	return '#'+'%02x%02x%02x' % rgb

def draw(modele):
      if modele=="hydrostatique":
            #affichage de la pression, couleurs en fonction des valeurs min et max
            for i in range(1,scene.fluid.numX-1):
                  for j in range(1,scene.fluid.numY-1):
                        #ratio proche de 1 ==> valeur de la pression proche du maximum
                        ratio=(scene.fluid.grid[i,j].p - scene.fluid.min_p) / (scene.fluid.max_p - scene.fluid.min_p)
                        #print(ratio)
                        #color=rgb_to_hex((int(255*ratio),0,int(255*(1-ratio)))) #bleu et rouge
                        color=rgb_to_hex((int(255*ratio),int(255*ratio),int(255*ratio))) #nuances de gris
                        Canevas.create_rectangle(scene.fluid.h*(i-1), scene.fluid.h*(j-1), scene.fluid.h*i, scene.fluid.h*j,fill=color,outline='')
      elif modele=="velocity":
            #affichage de la vitesse, couleurs en fonction des valeurs min et max
            ratios=[]
            #scene.fluid.max_v=1*10**-8
            for i in range(1,scene.fluid.numX-1):
                  for j in range(1,scene.fluid.numY-1):
                        if scene.fluid.grid[i,j].s != 0:
                              #ratio proche de 1 ==> valeur proche du maximum
                              ratio=(sqrt(scene.fluid.grid[i,j].u**2 + scene.fluid.grid[i,j].v**2) - scene.fluid.min_v) / (scene.fluid.max_v - scene.fluid.min_v)
                              if ratio<0:
                                    ratio=0
                              elif ratio>1:
                                    ratio=1
                              ratios.append(ratio)
                              #print(ratio)
                              color=rgb_to_hex((int(255*ratio),0,int(255*(1-ratio)))) #bleu et rouge
                              #color=rgb_to_hex((int(255*ratio),int(255*ratio),int(255*ratio))) #nuances de gris
                              Canevas.create_rectangle(scene.fluid.h*(i-1), scene.fluid.h*(j-1), scene.fluid.h*i, scene.fluid.h*j,fill=color,outline='')
            print(mean(ratios),median(ratios),scene.fluid.min_v,scene.fluid.max_v)
     



fenetre=Tk()
#fenetre.attributes('-fullscreen', True)
fenetre.bind('<Escape>',lambda e: fenetre.destroy())

Canevas = Canvas(fenetre, width=width, height=height)
Canevas.pack()

Bouton1 = Button(fenetre, text = 'Quitter', command = fenetre.destroy)
Bouton1.pack()

setup(modele)
for i in range(10):
      update(modele)
draw(modele)

fenetre.mainloop()







