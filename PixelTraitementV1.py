from astropy.io import fits
import matplotlib.pyplot as plt
import numpy 
import FitsFonctions as fts
import scipy.optimize as So  
import sys


def ShowAllImagesGraph(mode = "precalculated"):

    """ mode = precalculated/rowdata"""

    if mode == "precalculated" :    # donnée pré calculer -> permet de gagner du temps
        file = fits.open('../Data/dataSaves/withoutNoiseSave1.fits')
        withoutNoise = numpy.array(file[0].data)

    elif mode == "rowdata":  # re-calcul les données depuis celles de base

        lightDirectoryPath = '../Data/data_633nm'   # chemin d'accès des deux dossiers
        darkDirectoryPath = '../Data/dark_500cit'


        averageDark = fts.getAverageImages(darkDirectoryPath)   # charge la moyennes de tout ces types d'images
        print("Dark images Done !")
        averageLight = fts.getAverageImages(lightDirectoryPath)
        print("Light images Done !")
        withoutNoise = numpy.subtract(averageLight,averageDark)  # crée une image sans bruit de fond pour améliorer la qualité




    # crée les variables de plus importantes :

    Informations = fts.getTilesCoordinatesInformation((29,11))
    tilesCoordinates , midlePixel = Informations[0] , Informations[1]
    


    waveLengthList = fts.getwaveLength('../Data/mapWavelength_gauss.fits' , midlePixel)   # revoie la liste des 319 longueures d'onde associé à chaque tuile 


    # execute le programme

    
    Radius = 1
    MediumValues = []
    for pixel in range(319):    #  crée la liste de la valeur médiane de chaques tuile dans un certain rayon depuis le centre de celle-ci
        mediumValue = fts.getMediumValueAroundPixel(Radius, midlePixel[pixel],withoutNoise)
        MediumValues.append(mediumValue)


    # affichage du graphique donné :

    print(waveLengthList,len(waveLengthList))

    columns = 1
    rows = 2
    fig = plt.figure(figsize=(8, 8))
    # affichage de lu graphique de l'intensité en fonction de longueurs d'ondes
    fig.add_subplot(rows, columns, 1)
    plt.plot(waveLengthList,MediumValues,marker="P", linestyle='-')
    # affchage de l'image moyenne avec les pixels choisis 
    fig.add_subplot(rows, columns, 2)
    plt.imshow(withoutNoise)

    tempX , tempY = [] , []
    for value in midlePixel :   # crée des listes temporaire pour afficher l'emplacement de la prise d'information
        tempX.append(value[0])
        tempY.append(value[1])

    plt.plot(tempX , tempY, linestyle='None', marker='o', color='b')
    plt.show() 



    # x = [i for i in range(319)]
    # plt.plot(x,MediumValues)
    # plt.xlabel('Tiles Number')
    # plt.ylabel('Light Value')
    # plt.title("Light")
    # plt.show()





def ShowOnePixelGraph(PixelCoordinates):

    # file[1].header['LAMBDA']

    Files = ['../Data/scan/20220725_105601_scanZolix_00000.fits' , '../Data/scan/20220725_111456_scanZolix_00001.fits']
    Informations = fts.getTilesCoordinatesInformation((29,11))
    #tilesCoordinates , Midles = Informations[0] , Informations[1]


    #Center = Midles[tileNumber]
    waveLengths = []
    pixelValues = []
    for path in Files :     # parcours les fichiers de scan/
        multipleImage = fits.open(path)
        for i in range(1,len(multipleImage)):

            waveLengths.append( float(multipleImage[i].header['LAMBDA']) )   # on ajoute la longueur d'onde à la liste
            data = numpy.array(multipleImage[i].data)

            Radius = 4
            Stregth = data[PixelCoordinates[1]][PixelCoordinates[0]]  #fts.getMediumValueAroundPixel(Radius,Center,data)
            pixelValues.append(Stregth)

            print(f"{round(i*100  /len(multipleImage),2)} %" , end="\r")

        multipleImage.close()


            



          
    #gaussianValues = fts.getGaussianModel(heigth=786  ,  FWHm=4.31  ,  LambdaMax=619  ,  abscissa= waveLengths, polynome=())



    # affichage de lu graphique 

    heigth = max(pixelValues)
    Lambdamax = waveLengths[pixelValues.index(heigth)]

    print('heigth',heigth)
    print('Lambdamax',Lambdamax)

    factor = 0.05
    datasetX = [waveLengths[pixelValues.index(i)] for i in pixelValues if i > heigth*factor]
    datasetY = [i for i in pixelValues if i > heigth*factor]

    FWHM =( datasetX[-1] - datasetX[0])/4
    initGuess = [heigth , FWHM , Lambdamax]

    print("FWHM" , FWHM)
    print('datasetX',datasetX)
    print('datasetY',datasetY)

    constants , covariance = So.curve_fit(fts.calculateGaussian , datasetX , datasetY , p0=initGuess )   # , bounds=((0,Lambdamax*0.8,0.4),(heigth*1.2,Lambdamax*1.2,20))
    print(constants)
    print("waveLengths[0] : ",waveLengths[0])
    print("waveLengths[1] : ", waveLengths[1])
    
    
    gaussianValues = fts.getGaussianModel(heigth=constants[0] ,  FWHm=constants[1]  ,  LambdaMax=constants[2]  ,  abscissa= waveLengths, polynome=())
    Abssica = numpy.linspace(waveLengths[0],waveLengths[-1] , num=1500)
    ModelValues = fts.getGaussianModel(heigth=constants[0] ,  FWHm=constants[1]  ,  LambdaMax=constants[2]  ,  abscissa= Abssica, polynome=())


    
    columns = 1
    rows = 2
    fig = plt.figure(figsize=(8, 8))
    # affichage de lu graphique de l'intensité en fonction de longueurs d'ondes
    fig.add_subplot(rows, columns, 1)
    plt.plot(waveLengths,pixelValues,marker="P", linestyle='-')
    plt.plot(Abssica,ModelValues,marker="None", linestyle='-')
    # affchage de l'image moyenne avec les pixels choisis 

    
    fig.add_subplot(rows, columns, 2)
    plt.plot(Abssica,ModelValues,marker="None", linestyle='-')
    plt.show() 
    






























def CalculateRowOfPixel(firstIndex,length,shape,Imagecount,masque) :

    Files = ['../Data/scan/20220725_105601_scanZolix_00000.fits' , '../Data/scan/20220725_111456_scanZolix_00001.fits']

    waveLengths = []    # liste de toutes les longueurs d'ondes
    PixelsStegths = numpy.zeros((length,shape[1],Imagecount) , dtype=None)    # crée la variable sui sauvegardera toute les valeurs des pixels 
    imageindex = -1
    for path in Files :     # parcours les fichiers de scan/
        multipleImage = fits.open(path)
        for i in range(1,len(multipleImage)):   # iter toute les image du fichier

            imageindex += 1
            waveLengths.append( float(multipleImage[i].header['LAMBDA']) )   # on ajoute la longueur d'onde à la liste
            data = numpy.array(multipleImage[i].data)   # données de l'image 


            print(imageindex , length)
            for Ycoordinates in range(length) :
                for Xcoordinates in range(shape[1]) :   # parcours tout les pixel de l'image  

                    if masque[firstIndex+Ycoordinates][firstIndex+Xcoordinates] : # si le pixel se trouve dans le masque 
                        PixelsStegths[Ycoordinates][Xcoordinates][imageindex] = data[firstIndex+Ycoordinates][firstIndex+Xcoordinates]    # ajoute a la liste la valeur du pixel
                        # print([Ycoordinates,Xcoordinates,imageindex])
                        # print("len() ! ",len(PixelsStegths) , Ycoordinates )


            # print(f"{round(i*100  /len(multipleImage),2)} %" , end="\r")

        multipleImage.close()



    TemporaryModelArray = numpy.zeros((length,shape[1],3) , dtype=None)



    working = 0
    notworking = []

    for Ycoordinates in range(length):
        
        # print(round(Ycoordinates/shape[0]*100,2) , "%" , end='\r')
        for Xcoordinates in range(shape[1]) :   # parcours tout les pixel de l'image 
            if masque[Ycoordinates][Xcoordinates] :

                
                PixelData = PixelsStegths[Ycoordinates][Xcoordinates].tolist()   # données du pixel

                heigth = max(PixelData) # estimation de la hauteur maximal
                Lambdamax = waveLengths[PixelData.index(heigth)]    # esimation du x maximal

                # print('heigth',heigth)
                # print('Lambdamax',Lambdamax)

                factor = 0.05
                datasetX = [waveLengths[PixelData.index(i)] for i in PixelData if i > heigth*factor]
                datasetY = [i for i in PixelData if i > heigth*factor]    # crée deux nouveaux set de données basé sur la courbe la plus haute

                print('datasetX',datasetX)
                print('datasetY',datasetY)
                if datasetX == [] :
                    print('heigth',heigth)
                    print('Lambdamax',Lambdamax)


                FWHM =( datasetX[-1] - datasetX[0])/4
                initGuess = [heigth , FWHM , Lambdamax] # estimation initial des constantes de mesure 

                # print("FWHM" , FWHM)
                # print('datasetX',datasetX)
                # print('datasetY',datasetY)

                try :
                    constants , covariance = So.curve_fit(fts.calculateGaussian , datasetX , datasetY , p0=initGuess )   # , bounds=((0,Lambdamax*0.8,0.4),(heigth*1.2,Lambdamax*1.2,20))
                    # print(constants)
                    TemporaryModelArray[Ycoordinates][Xcoordinates] = constants    # on transcrit les valeurs trouvé 
                    working += 1
                    # print(Ycoordinates,Xcoordinates , "working")
                except :

                    notworking.append((Ycoordinates,Xcoordinates))
                    # print(Ycoordinates,Xcoordinates , "not working !")
                    # columns = 1
                    # rows = 1
                    # fig = plt.figure(figsize=(8, 8))
                    # # affichage de lu graphique de l'intensité en fonction de longueurs d'ondes
                    # fig.add_subplot(rows, columns, 1)
                    # plt.plot(waveLengths,PixelData,marker="P", linestyle='-')
                    # plt.show()

    

    print(working,len(notworking))
    print(notworking)

    return TemporaryModelArray












def CalculateEveryPixelModel(shape=(1096,2808),rowsCalculQuantity = 10):
    """
    shape : (y,x) 
    """

    # création du masque 

    Hzz = fits.open('../Data/mapWavelength_gauss.fits')
    HzzData = Hzz[0].data

    needToCalculate = numpy.ndarray(shape,dtype=bool)   # shape : (y,x)
    print(needToCalculate)
    masque = numpy.where(HzzData > 0 , True , False )
    masque = numpy.flip(masque,axis=0)  # flip le tableau sur l'axe horizontal (inverse le haut et le bas)


    file = fits.open('../Data/dataSaves/withoutNoiseSave1.fits')
    Files = ['../Data/scan/20220725_105601_scanZolix_00000.fits' , '../Data/scan/20220725_111456_scanZolix_00001.fits']

    
    # calcul du model de chaque pixel

    print(masque)
    FinalModelsArray = numpy.zeros((shape[0],shape[1],3) , dtype=None)
    print(FinalModelsArray)


    totalPixel = shape[1] *shape[0]

    Imagecount = 0 
    for i in Files :
        i =  fits.open(i)
        Imagecount += len(i)-1 
        i.close()


    
    for RowSection in range(0,shape[0],rowsCalculQuantity):
        length = min(rowsCalculQuantity,shape[0]-RowSection)
        print("length :",length , "RowSection" , RowSection)
        ReturnedArray = CalculateRowOfPixel(RowSection,length,shape,Imagecount,masque)

        for index,value in enumerate(ReturnedArray) :
            FinalModelsArray[RowSection+index] = value




    return FinalModelsArray



    # columns = 1
    # rows = 1
    # fig = plt.figure(figsize=(8, 8))
    # # affichage de lu graphique de l'intensité en fonction de longueurs d'ondes
    # ax = fig.add_subplot(rows, columns, 1)
    # ax.imshow(withoutNoise, aspect='auto', cmap="copper", interpolation='nearest')

    # plt.show()
    




print(CalculateEveryPixelModel(shape=(200,200),rowsCalculQuantity=100))
# fts.showOneImage('../Data/mapWavelength_gauss.fits')
# ShowOnePixelGraph((1236,430))