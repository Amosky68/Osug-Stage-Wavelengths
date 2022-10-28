from astropy.io import fits
import matplotlib.pyplot as plt
import os
import numpy
from math import pi,sqrt,exp



def getAverageImages(directoryPath):
    """ renvoie l'image moyenne d'une série situé dans un dossier """
    arraySize = (1096,2808)  # row , column

    averageArray = numpy.ndarray(arraySize,dtype=int)   # initialise le tableau
    allImagesPaths = os.listdir(directoryPath)  # liste les chemin d'accès de chaque timages du dossier


    for imageIndex,imagePath in enumerate(allImagesPaths) : # parcours toutes les images du dossier
        imagePath = directoryPath + '/' + imagePath
        file = fits.open(imagePath)     # ouvre les données

        imageData = numpy.array(file[0].data)   # récupère les données de celle-ci
        averageArray = numpy.add( averageArray , imageData) # fait l'addition de toutes les valeurs des donées et du tableau 
        print(f"{round(imageIndex*100/len(allImagesPaths),2)} %" , end="\r")    # affiche le poyrcentage de progression
        file.close()

    averageArray = numpy.divide(averageArray,len(allImagesPaths))   # divise toute les valeurs du tableau avec la quantité d'image pour en faire une moyenne

    return averageArray  # retourne la valeur finale



def showOneImage(imagePath) :
    """ Montre Une image basé sur son chemin d'accès """
    file = fits.open(imagePath) 
    imageData = file[0].data    # récupère les données de l'image en question

    plt.imshow(imageData , cmap='gray')     # affiche cette image
    file.close()
    plt.show()



def showImageWithArray(Array):
    """ Déssine l'image crée par une base de donnée """
    plt.imshow(Array , cmap='gray') 
    plt.show()



def showEveryForms():   
    """ 
    crée 3 images : 
    - Moyenne des images avec de la lumière (633 nm)
    - Moyenne des images sans lumière
    - soustration des deux images | une image avec le bruit de fond retiré 
    """

    lightDirectoryPath = '../Data/data_633nm'   
    darkDirectoryPath = '../Data/dark_500cit'


    averageDark = getAverageImages(darkDirectoryPath)
    averageLight = getAverageImages(lightDirectoryPath)
    withoutNoise = numpy.subtract(averageLight,averageDark)
    showImageWithArray(averageLight)
    showImageWithArray(averageDark)
    showImageWithArray(withoutNoise)



def getMediumValueAroundPixel(Radius,PixelCoordinates,Dataset): # revoie la valeur médiane des pixels autours d'un point et d'un rayons donné
    Radius = max(Radius,1)
    Pixelx,Pixely = PixelCoordinates

    valueList = []
    for y in range(Radius*2-1) :    # parcours tout les pixel dans le rayon correspondant
        for x in range(Radius*2-1):
            valueList.append(Dataset[y+Pixely][x+Pixelx])   # ajoute la valeur de ce pixel 

    valueList.sort()
    midle = int(len(valueList)/2) # trie la liste et revoie la valeur la plus au centre de celle ci
    return valueList[midle]



def getwaveLength(File,PixelCoordinatesList) :

    
    waveLength = []   # liste des longueurs d'ondes de chaques tuiles/pixel (319)
    file = fits.open(File)
    fileData = file[0].data

    for i,value in enumerate(PixelCoordinatesList) :    # crée la liste des médianes de longueurs d'onde de chaques tuiles (médiane pour être plus précis)
        waveLength.append(getMediumValueAroundPixel(8,value,fileData))

    return waveLength
        


def saveDataset(FilePath,Dataset):  # create a fits file with the dataset dans the file location

    data = fits.PrimaryHDU(Dataset) 
    data.writeto(FilePath)



def getTilesCoordinatesInformation(shape):

    tilesCoordinates = []   # crée la liste des coordonées de chaques tuiles (29*11)
    for i in range(shape[0] * shape[1]):    # car il y a 319 images
        y = i % shape[1]
        x = i // shape[1]
        tilesCoordinates.append((x,y))

    Offset = (13,13)
    tileSize = 96
    midlePixel = []
    for i,value in enumerate(tilesCoordinates) :    # crée la liste des coordonées des pixel du centre de chaque tuiles 
        x = int((value[0] + .5)*tileSize) + Offset[0]
        y = int((value[1] + .5)*tileSize) + Offset[1]
        midlePixel.append((x,y))

    return tilesCoordinates , midlePixel



def calculateGaussian(x , heigth ,FWHm ,LambdaMax ):

    return heigth * numpy.exp(-0.5 * ((x - LambdaMax)/FWHm)**2) 


def calculateLorentz(x , heigth ,FWHm ,LambdaMax ):

    return heigth/1+((x-LambdaMax)/FWHm)**2



def getGaussianModel(heigth:float,FWHm:float,LambdaMax:float,abscissa:list,polynome = (0,0)):
    """
    polynome : les coefficients pour la puissance de x qui augmente :
    2x² - 3x  -> (2,-3,0)
    """
    outputValues = []
    for x in abscissa :
        yValue = calculateGaussian(x , heigth ,FWHm ,LambdaMax)

        outputValues.append(yValue)

    return outputValues




def getLorenzModel(heigth:float,FWHM:float,LambdaMax:float,abscissa:list,polynome = ()):

    """
    polynome : les coefficients pour la puissance de x qui augmente :
    2x² - 3x  -> (2,-3,0)
    """


    outputValues = []
    for x in abscissa :
        yValue = calculateLorentz(x , heigth ,FWHM ,LambdaMax)

        for index,coef in enumerate(polynome) :
            yValue += coef * x**(len(polynome)-index-1)

        outputValues.append(yValue)

    return outputValues




