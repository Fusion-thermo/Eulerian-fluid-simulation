import numpy as np
from math import floor, sqrt
from copy import deepcopy

class Cell:
    def __init__(self,s=1):
        self.u=0 #horizontal velocity component
        self.v=0 #vertical velocity component
        self.s=s #1 pour les cellules à simuler, 0 pour les cellules de bord
        self.p=0 #pressure
        self.m = 1 #densité d'une cellule utilisée par advectSmoke

class Fluid:
    def __init__(self,numX,numY,grid_spacing,density,overRelaxation):
        self.overRelaxation = overRelaxation
        self.density=density
        self.numX=numX+2
        self.numY=numY+2
        self.h=grid_spacing
        #initialise la grille, d'abord les bords puis le centre
        self.grid=np.full((self.numX, self.numY), Cell())

        for i in range(self.numX):
            self.grid[i,0]=Cell(s=0)
            self.grid[i,self.numY-1]=Cell(s=0)
        for j in range(self.numY):
            self.grid[0,j]=Cell(s=0)
            self.grid[self.numX-1,j]=Cell(s=0)
        for i in range(1,self.numX-1):
            for j in range(1,self.numY-1):
                self.grid[i,j]=Cell()

    
    #Advection step
    
    def integrate(self,dt, gravity):
        for i in range(1,self.numX):
            for j in range(1,self.numY-1):
                #augmente la vitesse en fonction de la gravité si la case est une case de fluide est qu'elle n'a pas de bord en dessous
                if self.grid[i,j].s != 0 and self.grid[i,j+1] != 0:
                    self.grid[i,j].v += gravity * dt

    
    def solveIncompressibility(self,numIters,dt):
        #projection method
        #on peut sûrement améliorer le système en mettant en condition une valeur de divergence
        #au lieu d'un nombre d'itérations fixe
        cp = self.density * self.h / dt
        for iter in range(numIters):
            divs=0
            for i in range(1,self.numX-1):
                for j in range(1,self.numY-1):
                    if self.grid[i,j].s != 0:
                        #sum all the s values
                        sx0 = self.grid[i-1,j].s
                        sx1 = self.grid[i+1,j].s
                        sy0 = self.grid[i,j+1].s
                        sy1 = self.grid[i,j-1].s
                        s = sx0 + sx1 + sy0 + sy1
                        if s == 0:
                            continue
                        #compute the divergence
                        div = self.grid[i+1,j].u - self.grid[i,j].u + self.grid[i,j-1].v - self.grid[i,j].v
                        divs+=div
                        #multiply by the overRelaxation factor
                        p = -div / s * self.overRelaxation
                        self.grid[i,j].p += cp * p
                        #corrects the velocity component
                        self.grid[i,j].u -= sx0 * p
                        self.grid[i+1,j].u += sx1 * p
                        self.grid[i,j].v -= sy0 * p
                        self.grid[i,j-1].v += sy1 * p
            #print(divs/((self.numX-2)*(self.numY-2)))

    #rend chaque cellule de bordure égale à la cellule de fluide voisine                 
    def extrapolate(self):
        for i in range(0,self.numX):
            self.grid[i,0].u = self.grid[i,1].u
            self.grid[i,self.numY-1].u = self.grid[i,self.numY-2].u
        for j in range(0,self.numY):
            self.grid[0,j].v = self.grid[1,j].v
            self.grid[self.numX-1,j].v = self.grid[self.numX-2,j].v


    def sampleField(self,x, y, field):
        h=self.h
        h1= 1/h
        h2 = 0.5 * h

        x = max(min(x, self.numX * h), h)
        y = max(min(y, self.numY * h), h)

        dx = 0
        dy = 0

        if field=="U_FIELD":
            dy = h2
        elif field=="V_FIELD":
            dx = h2
        elif field=="S_FIELD":
            dx = h2
            dy = h2
        
        x0 = min(floor((x-dx)*h1), self.numX-1)
        tx = ((x-dx) - x0*h) * h1
        x1 = min(x0 + 1, self.numX-1)
        
        y0 = min(floor((y-dy)*h1), self.numY-1)
        ty = ((y-dy) - y0*h) * h1
        y1 = min(y0 + 1, self.numY-1)

        sx = 1.0 - tx
        sy = 1.0 - ty

        if field=="U_FIELD":
            var = sx*sy * self.grid[x0,y0].u + tx*sy * self.grid[x1,y0].u + tx*ty * self.grid[x1,y1].u + sx*ty * self.grid[x0,y1].u
        elif field=="V_FIELD":
            var = sx*sy * self.grid[x0,y0].v + tx*sy * self.grid[x1,y0].v + tx*ty * self.grid[x1,y1].v + sx*ty * self.grid[x0,y1].v
        elif field=="S_FIELD":
            #smoke field
            var = sx*sy * self.grid[x0,y0].m + tx*sy * self.grid[x1,y0].m + tx*ty * self.grid[x1,y1].m + sx*ty * self.grid[x0,y1].m
        
        return var
        
    def avgU(self,i,j):
        u = (self.grid[i,j+1].u + self.grid[i,j].u + self.grid[i+1,j+1].u + self.grid[i+1,j].u) / 4
        return u
    
    def avgV(self,i,j):
        v = (self.grid[i-1,j].v + self.grid[i,j].v + self.grid[i-1,j-1].v + self.grid[i,j-1].v) / 4
        return v

#advection
    def advectVel(self,dt):
        self.newGrid=deepcopy(self.grid)
        h = self.h
        h2 = 0.5 * h
        for i in range(1,self.numX):
            for j in range(0,self.numY-1):
                #u component
                if self.grid[i,j].s != 0 and self.grid[i-1,j].s != 0 and j<self.numY-1:
                    x=i*h
                    y=j*h +h2
                    u=self.grid[i,j].u
                    v=self.avgV(i,j)
                    x=x - dt*u
                    y=y - dt*v
                    u=self.sampleField(x,y,"U_FIELD")
                    self.newGrid[i,j].u = u
                #v component
                if self.grid[i,j].s != 0 and self.grid[i,j+1].s != 0 and i<self.numX-1:
                    x=i*h + h2
                    y=j*h
                    u=self.avgU(i,j)
                    v=self.grid[i,j].v
                    x=x - dt*u
                    y=y - dt*v
                    v=self.sampleField(x,y,"V_FIELD")
                    self.newGrid[i,j].v = v
        self.grid = deepcopy(self.newGrid)

    # pas nécessaire pour l'instant
    # def advectSmoke(self,dt):
    #     #retrace le parcours d'une pseudo particule en arrière
    #     self.newM=deepcopy(self.grid)

    #     h = self.h
    #     h2 = 0.5 * h
        
    #     for i in range(1,self.numX-1):
    #         for j in range(1,self.numY-1):
    #             if self.grid[i,j] != 0:
    #                 u = (self.grid[i,j].u + self.grid[i+1,j].u) * 0.5
    #                 v = (self.grid[i,j].v + self.grid[i,j-1].v) * 0.5
    #                 x = i*h + h2 - dt*u
    #                 y = j*h + h2 - dt*v

    #                 self.newGrid[i,j].m = self.sampleField(x,y,"S_FIELD")

    #     self.grid=deepcopy(self.newM)


    def simulate(self,dt, gravity, numIters):
        if gravity!=0:
            self.integrate(dt,gravity)

        for i in range(1,self.numX-1):
            for j in range(1,self.numY-1):
                self.grid[i,j].p = 0
        self.solveIncompressibility(numIters,dt)

        self.extrapolate()
        self.advectVel(dt)
        #self.advectSmoke(dt)

        pressures=[self.grid[i,j].p for i in range(1,self.numX-1) for j in range(1,self.numY-1)]
        self.min_p=min(pressures)
        self.max_p=max(pressures)
        velocity_magnitudes=[sqrt(self.grid[i,j].u**2 + self.grid[i,j].v**2) for i in range(1,self.numX-1) for j in range(1,self.numY-1)]
        self.min_v=min(velocity_magnitudes)
        self.max_v=max(velocity_magnitudes)

            


class Scene:
    def __init__(self):
        self.gravity = -9.81
        self.dt = 1 / 120
        self.numIters = 100
        self.frameNr = 0
        self.overRelaxation = 1.9
        self.obstacleX = 0
        self.obstacleY = 0
        self.obstacleRadius = 0.15
        self.paused = False
        self.sceneNr = 0
        self.showObstacle = False
        self.showStreamlines = False
        self.showVelocities = False
        self.showPressure = False
        self.showSmoke = True
        self.fluid = None