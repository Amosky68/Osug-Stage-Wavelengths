from symbol import parameters
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy 
import FitsFonctions as fts
import scipy.optimize as So  
import sys
import csv


def getTotalImageAmount(TreatmentFilesPaths):

    Imagecount = 0 
    for i in TreatmentFilesPaths :
        i =  fits.open(i)
        Imagecount += len(i)-1 
        i.close()
    
    return Imagecount



def getWaveLengths(TreatmentFilesPaths):

    WaveLengths = []

    for path in TreatmentFilesPaths :     # parcours les fichiers de scan/
        PrimaryFile = fits.open(path)
        for i in range(1,len(PrimaryFile)):   # iter toute les image du fichier

            WaveLengths.append( float(PrimaryFile[i].header['LAMBDA']) )   # on ajoute la longueur d'onde à la liste

    return WaveLengths





def SaveDataInFits(Dataset , SavingFilePath):
    
    print("Saving Datas")
    data = fits.PrimaryHDU(Dataset) 
    data.writeto(SavingFilePath)






def  CalculateModelFromRows(Rows:int  , BeginingRowIndex: int  , BeginingCoordinates: tuple ,  Totalshapes  ,  TreatmentFilesPaths  ,  MasqueArray  ):

    ParameterOutputArray = numpy.zeros((Rows , Totalshapes[1] , 3)) # y , x , z
    print(ParameterOutputArray)

    ImageAmount = getTotalImageAmount(TreatmentFilesPaths)
    print("ImageAmount : " , ImageAmount)

    WaveLengthList = getWaveLengths(TreatmentFilesPaths)


    
    PixelsStegths = numpy.zeros((Rows , Totalshapes[1] , ImageAmount))    # crée la variable suis sauvegardera toute les valeurs des pixels 
    imageindex = -1

    try :
        del ImageData   # on détruit ImageData pour libérer de l'espace ram 
    except :
        pass

    for path in TreatmentFilesPaths :     # parcours les fichiers de scan/
        PrimaryFile = fits.open(path)
        for i in range(1,len(PrimaryFile)):   # iter toute les image du fichier
            imageindex += 1

            ImageData = numpy.array(PrimaryFile[i].data)   # données de l'image 
            print(imageindex , end=' \r')


            for ArrayIndex in range(Rows):
                ImageRowIndex = ArrayIndex + BeginingRowIndex + BeginingCoordinates[1]
                for Xcoordinates in range(Totalshapes[1]) :   # parcours tout les pixel de l'image 
                    ImageColumnIndex = Xcoordinates + BeginingCoordinates[0] 

                    if MasqueArray[ImageRowIndex][ImageColumnIndex] :
                        PixelsStegths[ArrayIndex][Xcoordinates][imageindex] = ImageData[ImageRowIndex][ImageColumnIndex]



    FailureCount = [0,0,0]  # liste comptant le taux de reussite | [trouvé , non trouvé , total]

    for ArrayIndex in range(Rows):
        ImageRowIndex = ArrayIndex + BeginingRowIndex + BeginingCoordinates[1]
        for Xcoordinates in range(Totalshapes[1]) :   # parcours tout les pixel de l'image 
            ImageColumnIndex = Xcoordinates + BeginingCoordinates[0] 

            if MasqueArray[ImageRowIndex][ImageColumnIndex] :   # si le pixel se trouve sur une image

                PixelData = PixelsStegths[ArrayIndex][Xcoordinates].tolist()   # données du pixel
                FailureCount[2] += 1    # ajoute 1 au total de pixel à compter


                heigth = max(PixelData) # estimation de la hauteur maximal
                Lambdamax = WaveLengthList[PixelData.index(heigth)]    # esimation du x maximal
                factor = 0.05
                datasetX = [WaveLengthList[PixelData.index(i)] for i in PixelData if i > heigth*factor]
                datasetY = [i for i in PixelData if i > heigth*factor]    # crée deux nouveaux set de données basé sur la courbe la plus haute

                if datasetX == [] :
                    print('heigth',heigth)
                    print('Lambdamax',Lambdamax)
                FWHM =( datasetX[-1] - datasetX[0])/4
                initGuess = [heigth , FWHM , Lambdamax] # estimation initial des constantes de mesure 



                try :
                    constants , covariance = So.curve_fit(fts.calculateGaussian , datasetX , datasetY , p0=initGuess )   # , bounds=((0,Lambdamax*0.8,0.4),(heigth*1.2,Lambdamax*1.2,20))
                    # print(constants)
                    ParameterOutputArray[ArrayIndex][Xcoordinates] = constants    # on transcrit les valeurs trouvé 
                    FailureCount[0] += 1    #  ajoute 1 au pixel trouvé
                    
                except :

                    print(f"Le Model du pixel {(ImageRowIndex,ImageColumnIndex)} n'a pas pu être trouvé , {100*FailureCount[1]/FailureCount[0]} % d'érreur")
                    FailureCount[1] += 1    #  ajoute 1 au pixel non trouvé

    return ParameterOutputArray , FailureCount  # revoie les valeurs des paramètre et du taux de reussite
                


    

                    
    









def CalculateEveryPixel(MasqueFilePath:str , TreatmentFilesPaths : list , BorderPoints = ((0,0),(100,100))  , RowCalculQuantity=100 , SavingFilePath='Saving') :
    """
    MasqueFilePath : chemin d'accès vers le fichier du masque\n
    TreatmentFilesPaths : chemins d'accès au fichiers Fits à traiter    \n
    BorderPoints : Liste de deux points : en haut a gauche et en bas à droite qui délimiterons la zone 
    dans laquelle les calculs se feronts, le (0,0) étant en haut a gauche. | ((x,y),(x2,y2))
    """


    shape = (BorderPoints[1][1]-BorderPoints[0][1] , BorderPoints[1][0]-BorderPoints[0][0] )    #shape : (y,x)
    BeginingCoordinates = BorderPoints[0]   # coordonées du début du calcul | x,y


    masquefile = fits.open(MasqueFilePath)
    masquedata = masquefile[0].data

    # needToCalculate = numpy.ndarray(shape,dtype=bool)   # shape : (y,x)
    masque = numpy.where(masquedata > 0 , True , False )
    masque = numpy.flip(masque,axis=0)  # flip le tableau sur l'axe horizontal (inverse le haut et le bas)

    ParametersList = [[] for i in range(shape[0])]  # liste final des paramètre trouvé sur chaque pixels




    for SectionIndex in range(shape[0]//RowCalculQuantity+1):   # on répète l'opération pour chaque section de plusieurs lignes 
        ToCalculateRows = shape[0] - RowCalculQuantity * SectionIndex # quantité de ligne a calculer restante 
        ToCalculateNow  = min(ToCalculateRows,RowCalculQuantity)   # quantité de ligne a calculer maintenant
        BeginingRowIndex = RowCalculQuantity * SectionIndex + BorderPoints[0][1]  # index du début de la ligne a calculer

        print(ToCalculateNow,BeginingRowIndex,BeginingCoordinates ,SectionIndex )
        if ToCalculateNow != 0 :
            ParameterOutputArray , SuccesRate = CalculateModelFromRows(ToCalculateNow,BeginingRowIndex,BeginingCoordinates,shape,TreatmentFilesPaths,masque)
            for Index in range(ToCalculateNow): # on ajoute à la liste final les données calculer par la fonctions 
                RowIndex = Index + BeginingRowIndex
                ParametersList[RowIndex] = ParameterOutputArray[Index]
            
            del ParameterOutputArray    # on détruit ParameterOutputArray pour libérer de l'espace ram


    print(ParametersList , SuccesRate)
    ParametersList = numpy.array(ParametersList)

    SaveDataInFits(ParametersList,SavingFilePath)






MasqueFilePath  = '../Data/mapWavelength_gauss.fits'
TreatmentFilesPaths = ['../Data/scan/20220725_105601_scanZolix_00000.fits' , '../Data/scan/20220725_111456_scanZolix_00001.fits']
RowCalculQuantity = 50
BorderPoints = ((0,0),(2808,100)) 
SavingFilePath = '../Data/dataSaves/SavedData2808100.fits'
CalculateEveryPixel(MasqueFilePath= MasqueFilePath  , TreatmentFilesPaths=TreatmentFilesPaths , RowCalculQuantity=RowCalculQuantity , BorderPoints = BorderPoints , SavingFilePath=SavingFilePath)
# print(getWaveLengths(TreatmentFilesPaths))