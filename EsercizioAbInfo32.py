from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os

Path0 = input("Inserisci il path completo della cartella .tgz da estrarre\n Esempio C:\\Users\\userName\\... \n")
NOME = input("\n Inserisci il nome della cartella .tgz da estrarre:\n")
Path=input("\nInserisci il path completo della destinazione dell'estrazione. \nEsempio C:\\Users\\userName\\... \n")
Path0= "cd "+Path0
Pathtgz="tar zxf "+ NOME + " -C " + Path
print("Estraendo...") 
os.system(Path0)
os.system(Pathtgz)
print("Estrazione completata..")
kk=1
##Path="C:/Users/mtcic/OneDrive/Desktop/perAbInf"
while (kk <= 3):
    test=kk
    Nbins=200
    NDati=100
    Nmisure=100
    measures1=[]
    measures2=[]
    measures3=[]
    #Inizio primo punto
    for i in np.arange(NDati)+1:
        j=0
        j1="XI0"
        j2="XI2"
        j3="XI4"
        fname = Path + f'/data/MockMeasures_2PCF_Test{test}/MockMeasures_2PCF_Correlation_MULTIPOLES_Test{test}_{i}.fits'
        file = fits.open(fname)
        table = file[1].data.copy()
        measures1.append(table[j1])
        measures2.append(table[j2])
        measures3.append(table[j3])
        if i==1:
            scale = table['SCALE']
        del table
        file.close()
    #Fine primo punto
    #Inizio secondo punto
    measures1=np.asarray(measures1).transpose()		#np.asarray converte measures in un array, transpose traspone da riga a colonna (o viceversa). Poiché numpy.cov richiede una colonna come argomento è necessario trasporre
    measures2=np.asarray(measures2).transpose()	
    measures3=np.asarray(measures3).transpose()	
    AVE1 = np.zeros((Nbins,),dtype=float)                                       #inizializza vettore con tanti 0 quanto Nbins
    AVE2 = np.zeros((Nbins,),dtype=float)                       
    AVE3 = np.zeros((Nbins,),dtype=float)                                 
    COV11 = np.zeros((Nbins,Nbins),dtype=float)                                  #inizializza matrice N x N con tanti 0 quanto Nbins x Nbins
    COV22 = np.zeros((Nbins,Nbins),dtype=float) 
    COV33 = np.zeros((Nbins,Nbins),dtype=float) 
    COV12 = np.zeros((Nbins,Nbins),dtype=float) 
    COV13 = np.zeros((Nbins,Nbins),dtype=float) 
    COV23 = np.zeros((Nbins,Nbins),dtype=float)  
    COV21 = np.zeros((Nbins,Nbins),dtype=float) 
    COV31 = np.zeros((Nbins,Nbins),dtype=float) 
    COV32 = np.zeros((Nbins,Nbins),dtype=float)  
    for i in range(Nmisure):                                                 #per calcolare la media si fa scorrere un indice da 0 a Nmeasures
        AVE1 += measures1[:,i]                                                   #la media sarà data dalla somma di tutti i valori contenuti in measures (il quale è un vettore che contiene)
    AVE1 /= Nmisure
    for i in range(Nmisure):   
        AVE2 += measures2[:,i]     
    AVE2 /= Nmisure
    for i in range(Nmisure):   
        AVE3 += measures3[:,i]     
    AVE3 /= Nmisure
    for i in range(Nbins):
        for j in range(Nbins):
            COV11[i,j] = (np.sum(measures1[i]*measures1[j]) - AVE1[i]*AVE1[j]*Nmisure) / (Nmisure-1)
            COV22[i,j] = (np.sum(measures2[i]*measures2[j]) - AVE2[i]*AVE2[j]*Nmisure) / (Nmisure-1)
            COV33[i,j] = (np.sum(measures3[i]*measures3[j]) - AVE3[i]*AVE3[j]*Nmisure) / (Nmisure-1)
            COV12[i,j] = (np.sum(measures1[i]*measures2[j]) - AVE1[i]*AVE2[j]*Nmisure) / (Nmisure-1)
            COV13[i,j] = (np.sum(measures1[i]*measures3[j]) - AVE1[i]*AVE3[j]*Nmisure) / (Nmisure-1)
            COV23[i,j] = (np.sum(measures2[i]*measures3[j]) - AVE2[i]*AVE3[j]*Nmisure) / (Nmisure-1)
    COV21 = COV12.transpose()
    COV31 = COV13.transpose()
    COV32 = COV23.transpose()
    # correlation matrix
    corr_xi11 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi11[i,j]=COV11[i,j]/(COV11[i,i]*COV11[j,j])**0.5
    corr_xi22 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi22[i,j]=COV22[i,j]/(COV22[i,i]*COV22[j,j])**0.5
    corr_xi33 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi33[i,j]=COV33[i,j]/(COV33[i,i]*COV33[j,j])**0.5
    corr_xi12 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi12[i,j]=COV12[i,j]/(COV12[i,i]*COV12[j,j])**0.5
    corr_xi13 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi13[i,j]=COV13[i,j]/(COV13[i,i]*COV13[j,j])**0.5
    corr_xi23 = np.zeros((Nbins,Nbins),dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            corr_xi23[i,j]=COV23[i,j]/(COV23[i,i]*COV23[j,j])**0.5
    TEST_COVARIANCE=True
    PLOTS=True
    if test==1:
        sigs = [0.02, 0.02, 0.02]
        ls = [25, 50, 75]
    elif test==2:
        sigs = [0.02, 0.01, 0.005]
        ls = [50, 50, 50]
    else:
        sigs = [0.02, 0.01, 0.005]
        ls = [5, 5, 5]
    #Fine secondo punto
    #Inizio terzo punto
    ## Definitions to build the covariance matrices based on Squared Exponential kernel
    def covf(x1, x2, sig, l):
        return sig**2.*np.exp(-(x1 - x2)**2./(2.*l**2.))
    def covf1f2(x1, x2, sig1, l1, sig2, l2):
        return (np.sqrt(2.*l1*l2)*np.exp(-(np.sqrt((x1 - x2)**2.)**2./(l1**2. + l2**2.)))*sig1*sig2)/np.sqrt(l1**2. + l2**2.)
    cov_th11 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th22 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th33 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th12 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th13 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th21 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th31 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th23 = np.zeros((Nbins,Nbins),dtype=float)
    cov_th32 = np.zeros((Nbins,Nbins),dtype=float)
    #for i in range(Nbins):
    #    for j in range(Nbins):
    #        cov_th[i,j] = covf(scale[i],scale[j],sigs[0],ls[0])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th11[i, j] = covf1f2(scale[i], scale[j], sigs[0], ls[0], sigs[0], ls[0])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th22[i, j] = covf1f2(scale[i], scale[j], sigs[1], ls[1], sigs[1], ls[1])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th33[i, j] = covf1f2(scale[i], scale[j], sigs[2], ls[2], sigs[2], ls[2])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th12[i, j] = covf1f2(scale[i], scale[j], sigs[0], ls[0], sigs[1], ls[1])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th13[i, j] = covf1f2(scale[i], scale[j], sigs[0], ls[0], sigs[2], ls[2])
    for i in range(Nbins):
        for j in range(Nbins):
            cov_th23[i, j] = covf1f2(scale[i], scale[j], sigs[1], ls[1], sigs[2], ls[2])
    cov_th21 = cov_th12.transpose()
    cov_th32 = cov_th23.transpose()
    cov_th31 = cov_th13.transpose()
    #Fine terzo punto
    #Inizio quinto punto
    FULL_MEASURED_MATRIX=np.zeros((3*Nbins,3*Nbins),dtype=float)
############################################## Creare matrice singola di sottomatrici per la cov misurata
    for i in range(Nbins):
        for j in range(Nbins):
            for t1 in range(3):
                for t2 in range (3):
                    I=i+t1*Nbins
                    J=j+t2*Nbins
                    if I<Nbins and J<Nbins:
                        FULL_MEASURED_MATRIX[I,J]= COV11[i,j]
                    if I<2*Nbins and I>Nbins and J<2*Nbins and J>Nbins:
                        FULL_MEASURED_MATRIX[I,J]= COV22[i,j]
                    if I>2*Nbins and J>2*Nbins:
                        FULL_MEASURED_MATRIX[I,J]= COV33[i,j]
                    if I<Nbins and J<2*Nbins and J>Nbins:
                        FULL_MEASURED_MATRIX[I,J]= COV12[i,j]
                    if I<Nbins and J>2*Nbins:
                        FULL_MEASURED_MATRIX[I,J]=COV13[i,j]
                    if I<2*Nbins and I>Nbins and J<Nbins:
                        FULL_MEASURED_MATRIX[I,J]=COV21[i,j]
                    if I<2*Nbins and I>Nbins and J>2*Nbins:
                        FULL_MEASURED_MATRIX[I,J]=COV23[i,j]
                    if I>2*Nbins and J<Nbins:
                        FULL_MEASURED_MATRIX[I,J]=COV31[i,j]
                    if I>2*Nbins and J<2*Nbins and J>Nbins:
                        FULL_MEASURED_MATRIX[I,J]=COV32[i,j]
    ############################################## Creare matrice singola di sottomatrici per la cov teorica
    FULL_THEORETICAL_MATRIX=np.zeros((3*Nbins,3*Nbins), dtype=float)
    for i in range(Nbins):
        for j in range(Nbins):
            for t1 in range(3):
                for t2 in range (3):
                    I=i+t1*Nbins
                    J=j+t2*Nbins
                    if I<Nbins and J<Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]= cov_th11[i,j]
                    if I<2*Nbins and I>Nbins and J<2*Nbins and J>Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]= cov_th22[i,j]
                    if I>2*Nbins and J>2*Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]= cov_th33[i,j]
                    if I<Nbins and J<2*Nbins and J>Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]= cov_th12[i,j]
                    if I<Nbins and J>2*Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]=cov_th13[i,j]
                    if I<2*Nbins and I>Nbins and J<Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]=cov_th21[i,j]
                    if I<2*Nbins and I>Nbins and J>2*Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]=cov_th23[i,j]
                    if I>2*Nbins and J<Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]=cov_th31[i,j]
                    if I>2*Nbins and J<2*Nbins and J>Nbins:
                        FULL_THEORETICAL_MATRIX[I,J]=cov_th32[i,j]
    if PLOTS:
        gratio = (1. + 5. ** 0.5) / 2.
        dpi = 300
        #climit=max(np.max(theoretical_covariance),np.max(measured_covariance))
        cmin = -np.max(cov_th11)*0.05
        cmax =  np.max(cov_th11)*1.05
        # Matrix plot of measured covariance matrix
        fig = plt.figure(figsize=(6,4))
        plt.title('measured covariance matrix')
        plt.imshow(FULL_MEASURED_MATRIX, vmin=cmin, vmax=cmax)
        cbar = plt.colorbar(orientation="vertical", pad=0.02)
        cbar.set_label(r'$ C^{\xi}_{N}$')
        # PLOTNAME = 'Test%s_Measured_Matrix.png'%test
        # plt.savefig(PLOTNAME,dpi = dpi)
        plt.show()
        # Matrix plot of theoretical covariance matrix
        fig = plt.figure(figsize=(6,4))
        plt.title('theoretical covariance matrix')
        plt.imshow(FULL_THEORETICAL_MATRIX, vmin=cmin, vmax=cmax)
        cbar = plt.colorbar(orientation="vertical", pad=0.02)
        cbar.set_label(r'$ C^{\xi}_{N}$')
        # PLOTNAME = 'Test%s_Measured_Matrix.png'%test
        # plt.savefig(PLOTNAME,dpi = dpi)
        plt.show()
        # Matrix plot of residuals
        fig = plt.figure(figsize=(6,4))
        plt.title('residuals')
        plt.imshow(FULL_THEORETICAL_MATRIX-FULL_MEASURED_MATRIX, vmin=cmin, vmax=-cmin)
        cbar = plt.colorbar(orientation="vertical", pad=0.02)
        cbar.set_label(r'$ C^{\xi}_{N}$')
        # PLOTNAME = 'Test%s_Measured_Matrix.png'%test
        # plt.savefig(PLOTNAME,dpi = dpi)
        plt.show()
    #Fine quinto punto
    #Inizo quarto punto
    covteo = []
    norm_residuals11 = np.zeros_like(cov_th11)
    norm_residuals22 = np.zeros_like(cov_th11)
    norm_residuals33 = np.zeros_like(cov_th11)
    norm_residuals12 = np.zeros_like(cov_th11)
    norm_residuals13 = np.zeros_like(cov_th11)
    norm_residuals23 = np.zeros_like(cov_th11)
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th11[i,j]**2./(np.sqrt(cov_th11[i,i]*cov_th11[j,j])**2.)
            norm_residuals11[i,j]=(cov_th11[i,j]-COV11[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th11[i,i]*cov_th11[j,j]))
    rms_deviation=np.std(norm_residuals11.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th22[i,j]**2./(np.sqrt(cov_th22[i,i]*cov_th22[j,j])**2.)
            norm_residuals22[i,j]=(cov_th22[i,j]-COV22[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th22[i,i]*cov_th22[j,j]))
    rms_deviation=np.std(norm_residuals22.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th33[i,j]**2./(np.sqrt(cov_th33[i,i]*cov_th33[j,j])**2.)
            norm_residuals33[i,j]=(cov_th33[i,j]-COV33[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th33[i,i]*cov_th33[j,j]))
    rms_deviation=np.std(norm_residuals33.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th12[i,j]**2./(np.sqrt(cov_th12[i,i]*cov_th12[j,j])**2.)
            norm_residuals12[i,j]=(cov_th12[i,j]-COV12[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th12[i,i]*cov_th12[j,j]))
    rms_deviation=np.std(norm_residuals12.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th13[i,j]**2./(np.sqrt(cov_th13[i,i]*cov_th13[j,j])**2.)
            norm_residuals13[i,j]=(cov_th13[i,j]-COV13[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th13[i,i]*cov_th13[j,j]))
    rms_deviation=np.std(norm_residuals13.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    for i in range(Nbins):
        for j in range(Nbins):
            rho2 = cov_th23[i,j]**2./(np.sqrt(cov_th23[i,i]*cov_th23[j,j])**2.)
            norm_residuals23[i,j]=(cov_th23[i,j]-COV23[i,j])*np.sqrt((NDati-1.)/((1.+rho2)*cov_th23[i,i]*cov_th23[j,j]))
    rms_deviation=np.std(norm_residuals23.reshape(Nbins**2))
    print(f"rms deviation of normalized residuals: {rms_deviation}")
    #Fine quarto punto
    if rms_deviation<1.0:
        print("Test root-mean-square:")
        print("PASSATO")
    else:
        print("Test root-mean-square")
        print("FALLITO")
    kk=kk+1

