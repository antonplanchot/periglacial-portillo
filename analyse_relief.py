import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

#import xdem une bibliotheque pour faire des calculs de terrain
#https://xdem.readthedocs.io/en/stable/advanced_examples/plot_slope_methods.html


#############################
#importation ALOS
pas=30.92207
path='C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/data_dem_index/'
fp1 = path+'output_DEM_ALOS_dcp_plus1pixel.csv'
fp2 = path+'output_DEM_ALOS_dcp.csv'

n,p=133+2,238+2
dem1=pd.read_csv(fp1,usecols=[2])
_1=pd.read_csv(fp2,usecols=[0,1])
dem2=np.array(dem1)
_2=np.array(_1)
dem=np.reshape(dem2,(n,p))

#############################
#création indices lat/long
Lnm=[238,133]
step=0.000277777778535437
Lccord=[[-69.6150054931641051,-69.5488815307617045],[-33.6458358764648011,-33.6088867187500000]]
for k in range(2):
    Lccord[k][0]+=step
    Lccord[k][1]-=step


step2=37
Llat_long=[0,0]
for k in range(2):
    a=Lccord[k][1]-Lccord[k][0]
    Llat_long[k]=[Lccord[k][0]+i*a/Lnm[k] for i in range(0,Lnm[k],step2)]
    Llat_long[k]=[str(np.round(x,2)) for x in Llat_long[k]]
print(Llat_long)

r=np.cos(33.64*np.pi/180)*6471000
xx=1000/(2*np.pi*r/360)
print(round(xx*100*step2))
############################
#importation cartographie final
map=pd.read_csv("C://Users//Utilisateur//Documents//_Internship_Anton//cartography_Portillo//treatment//attribut_intersection.csv",sep=',',usecols=[1,12,13,14])
map1=np.array(map)
nmap,mmap=map1.shape

L_index_int=[]
for k in range(nmap):
    ind=int(map1[k,1]*238+map1[k,2])
    if map1[k,0]!=0 and ind not in L_index_int :
        L_index_int.append(ind)
nb_int=len(L_index_int)
L_index_ext=[x for x in range(133*238) if x not in L_index_int]

def rognage(M,L_index_ext):
    M=np.reshape(M,(133*238,1))
    M[L_index_ext]=np.nan
    M=np.reshape(M,(133,238))


#############################
#DEM
ddem=pd.read_csv(fp2,usecols=[2])
ddem=np.array(ddem)

ddem=ddem.astype('float64')
ddem[L_index_ext]=np.nan
ddem=np.reshape(ddem,(133,238))

print(np.nanmin(ddem),np.nanmin(ddem))


def plot_map(step2,Llat_long,label1,title_image,M,cmap):
    size=16
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)
    plt.rc('axes', axisbelow=True)
    plt.figure(figsize=(17,7.7))
    plt.grid(zorder=1)
    plt.imshow(M,cmap=cmap,zorder=2)
    plt.colorbar(label=label1)
    plt.plot([184,223],[111,111],linewidth=7.0,color='k')
    plt.text(195,105,'1 km', fontsize=size)
    plt.xticks(range(0,238,step2),Llat_long[0])
    plt.yticks(range(0,133,step2),reversed(Llat_long[1]))
    plt.xlabel('longitude (°)')
    plt.ylabel('latitude (°)')
    plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/"+title_image+".png",bbox_inches="tight")
    plt.show()
    

plot_map(step2,Llat_long,'altitude (m)','demALOS',ddem,'viridis')

def dem2maat(ddem):
    z_station=3950
    maat_station=-1.285
    grad=-0.006
    return (ddem-z_station)*grad+maat_station

maat=dem2maat(ddem)

plot_map(step2,Llat_long,'mean annual air temperature (°C)','maat',maat,'viridis')


a=np.concatenate((_2,np.reshape(ddem,(133*238,1))),axis=1)
np.savetxt(path+'output_DEM_ALOS_dcp_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')




#############################
#PENTE

#Horn
M=np.zeros((n-2,p-2))
for i in range(1,n-1):
    for j in range(1,p-1):
        d1=((dem[i-1,j-1]+2*dem[i-1,j]+dem[i-1,j+1])-(dem[i+1,j-1]+2*dem[i+1,j]+dem[i+1,j+1]))/(8*pas)
        d2=((dem[i-1,j-1]+2*dem[i,j-1]+dem[i+1,j-1])-(dem[i-1,j+1]+2*dem[i,j+1]+dem[i+1,j+1]))/(8*pas)
        M[i-1,j-1]=np.sqrt(d1**2+d2**2)

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_slope_ALOS_dcp.csv', a, delimiter=',',header='x,y,slope', comments='')

rognage(M,L_index_ext)

plot_map(step2,Llat_long,'slope (m/m)','slope',M,'viridis')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_slope_ALOS_dcp_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')



if False:
    #Zevenbergen_Thorne  
    M=np.zeros((n-2,p-2))
    for i in range(1,n-1):
        for j in range(1,p-1):
            d1=(-dem[i,j-1]+dem[i,j+1])/(2*pas)
            d2=(dem[i-1,j]-dem[i+1,j])/(2*pas)
            M[i-1,j-1]=np.sqrt(d1**2+d2**2)
            
    rognage(M,L_index_ext)        
    plt.imshow(M)
    plt.colorbar(label='pente (m/m)')
    plt.show() 
  
#############################
#ORIENTATION

#ORIENTATION
#Zevenbergen_Thorne  /  Wilson and Gallant (2000)
M=np.zeros((n-2,p-2))
for i in range(1,n-1):
    for j in range(1,p-1):
        d1=(-dem[i,j-1]+dem[i,j+1])/(2*pas)
        d2=(dem[i-1,j]-dem[i+1,j])/(2*pas)
        if d1==0 or np.sqrt(d1**2+d2**2)<0.05:
            M[i-1,j-1]=np.nan
        else:
            M[i-1,j-1]=180-np.arctan(d2/d1)*(180/np.pi)+90*(d1/np.abs(d1))

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_orientation_ALOS_dcp.csv', a, delimiter=',',header='x,y,slope', comments='')

rognage(M,L_index_ext)        

plot_map(step2,Llat_long,'orientation (°)','orientation',M,'hsv')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_orientation_ALOS_dcp_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')

###############################
#COURBURE V ET H

#Zevenbergen_Thorne  
M1=np.zeros((n-2,p-2))
M2=np.zeros((n-2,p-2))
for i in range(1,n-1):
    for j in range(1,p-1):
        D=((dem[i,j-1]+dem[i,j+1])/2-dem[i,j])/pas**2
        E=((dem[i-1,j]+dem[i+1,j])/2-dem[i,j])/pas**2
        F=(-dem[i-1,j-1]+dem[i-1,j+1]+dem[i+1,j-1]-dem[i+1,j+1])/(4*pas**2)
        G=(-dem[i,j-1]+dem[i,j+1])/(2*pas)
        H=(dem[i-1,j]-dem[i+1,j])/(2*pas)
        if G**2+H**2>0:
            CV=2*(D+G**2+E*H**2+F*G*H)/(G**2+H**2)
            CH=-2*(D+H**2+E*G**2+F*G*H)/(G**2+H**2)
        else:
            CV,CH=np.nan,np.nan
        if -0.5<CV<2.5:
            M1[i-1,j-1]=CV
        else:
            M1[i-1,j-1]=np.nan
        if -2.5<CH<0.5:
            M2[i-1,j-1]=CH
        else:
            M2[i-1,j-1]=np.nan
        
rognage(M1,L_index_ext)        
plot_map(step2,Llat_long,'CurvatureV','curvature',M1,'binary')
plt.show()

if False:
    rognage(M2,L_index_ext)
    plt.imshow(M2,cmap='binary')
    plt.colorbar(label='CurvatureH')
    plt.show()




#############################
#TERRAIN RUGOSITY INDEX (TRI)

M=np.zeros((n-2,p-2))
for i in range(1,n-1):
    for j in range(1,p-1):
        neighbour=dem[i-1:i+2,j-1:j+2]
        neighbour=np.reshape(neighbour,(1,9))[0]
        S=0
        for x in neighbour:
            S+=(dem[i,j]-x)**2
        M[i-1,j-1]=np.sqrt(S)
        
a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_TRI_ALOS_dcp.csv', a, delimiter=',',header='x,y,slope', comments='')

rognage(M,L_index_ext)        


plot_map(step2,Llat_long,'TRI (m)','TRI',M,'viridis')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_TRI_ALOS_dcp_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')


#TRI par classes 
# M1=np.ones((n-2,p-2))*2
# M1[M<55]=1
# M1[M<15]=0
# rognage(M1,L_index_ext) 
deb=45
nrange=30
choix=57
for k in range(nrange):
    M1=np.ones((n-2,p-2))
    M1[M<deb+k]=0
    rognage(M1,L_index_ext) 
    
    if deb+k==choix:
        plot_map(step2,Llat_long,'range TRI','range_TRI',M1,'viridis')
        M1b=np.reshape(M1,(133*238,))
        L_TRI_bedrock=list(np.where(M1b==1)[0])
        
    
    a=np.concatenate((_2,np.reshape(M1,(133*238,1))),axis=1)
    np.savetxt(path+'output_rangeTRI'+str(deb+k)+'.csv', a, delimiter=',',header='x,y,slope', comments='')


'''remarque: cela discrimine bien le bedrock, ormis quelques irrégularités 
avec les debris slope'''

#############################
#SHADOW ROGNER

fp1 = path+'output_shadow_qgis_ALOS.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M[M!=0]#quelques 0 sur les bords
M=M.astype('float64')

M[L_index_ext]=np.nan
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'index shadow','shadow',M,'viridis')


a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_shadow_qgis_ALOS_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')

#INSOLATION ROGNER

fp1 = path+'output_total_insolation.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M=np.reshape(M,(135,240))
M=M[1:134,1:239]
M=np.reshape(M,(133*238,1))
M[L_index_ext]=np.nan
M=np.reshape(M,(133,238))


plot_map(step2,Llat_long,'Potential incoming solar radiation (kWh/m2)','PISR',M,'viridis')


a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_total_insolation_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')


#NDVI ROGNER

fp1 = path+'output_Pleiade_NDVI_2022_grilleALOS_dcp.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M[L_index_ext]=np.nan
M[M<0.025]=np.nan #enlever les ombres et l'eau
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'NDVI Pléiade 1st march 2022','NDVI',M,'viridis')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_Pleiade_NDVI_2022_grilleALOS_dcp_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')


#INDEX PERMAFROST ROGNER

fp1 = path+'output_index_permafrost.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M[L_index_ext]=np.nan
if False:
    M[L_TRI_bedrock]=np.nan
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'index permafrost field','permafrost',M,'copper')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_index_permafrost_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')


##########################################################
#BEDROCK

fp1 = path+'output_bedrock.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M[L_index_ext]=np.nan
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'bedrock area fraction (%)','bedrock',M*100,'viridis')

a=np.concatenate((_2,np.reshape(M,(133*238,1))),axis=1)
np.savetxt(path+'output_bedrock_rogne.csv', a, delimiter=',',header='x,y,slope', comments='')




#############################
#-------------------------------------------------------------------------
#############################
#ANALYSE CROISEE BEDROCK VS RANGE-TRI

n_fichier=31654
data=np.zeros((n_fichier,nrange+1))
for k in range(nrange):
    data[:,k]=np.array(pd.read_csv(path+'output_rangeTRI'+str(deb+k)+'.csv',usecols=[2]))[:,0]
data[:,-1]=np.array(pd.read_csv(path+'output_bedrock_rogne'+'.csv',usecols=[2]))[:,0]

data=data[L_index_int,:]

Lcoef=[]
Lnonspecif=[]
Ltropspecif=[]
for k in range(nrange):
    coef=np.corrcoef(data[:,k],data[:,-1])[0,1]
    Lcoef.append(coef**2)
    Mnonspecif=data[:,0][data[:,k]-data[:,-1]>0.5]#spécificité
    #on va compter le nombre de cellule sur lequel l'indice rangeTRI est supérieur de 0.5 au bedrock
    Mtropspecif=data[:,0][data[:,-1]-data[:,k]>0.5]
    Lnonspecif.append(100*Mnonspecif.shape[0]/nb_int)
    Ltropspecif.append(100*Mtropspecif.shape[0]/nb_int)
    
x=range(deb,deb+nrange)


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('threshold TRI')
ax1.set_ylabel('R^2', color=color)
ax1.plot(x, Lcoef, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('no specificity (%)', color=color)  # we already handled the x-label with ax1
ax2.plot(x, Lnonspecif, color=color,label='nonspé')
ax2.plot(x, Ltropspecif, linestyle='dashed',color=color,label='tropspé')
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.axvline(x=57, color='k')
plt.text(57.25, 22.5, 'x=57')
plt.legend()
plt.show()



#############################
#ANALYSE CROISEE

L_fichier_ctxt=['output_DEM_ALOS_dcp_rogne','output_slope_ALOS_dcp','output_orientation_ALOS_dcp','output_TRI_ALOS_dcp','output_total_insolation_rogne','output_Pleiade_NDVI_2022_grilleALOS_dcp_rogne','output_index_permafrost_rogne']
L_name_ctxt=['altitude (m)','slope (m/m)','orientation (°)','terrain rugosity index','insolation_sg22 (kWh/m2)','NDVI 2022','index permafrost field']
L_name_ctxt3=['altitude','slope','orientation','terrain rugosity index','insolation saga2022','NDVI 2022','index permafrost field']


n_fichier,m_fichier=31654,len(L_fichier_ctxt)
data=np.zeros((n_fichier,m_fichier))
for k in range(m_fichier):
    data[:,k]=np.array(pd.read_csv(path+L_fichier_ctxt[k]+'.csv',usecols=[2]))[:,0]

#specificite orientation
for k in range(n_fichier):
    if data[k,2]>250:
        data[k,2]=data[k,2]-360

Lcoef_perm=[]
Mcoef=np.zeros((m_fichier,m_fichier))
for i in range(m_fichier):
    for j in range(i+1,m_fichier):
        #plt.plot(data[:,i],data[:,j],'.')
        a,b=data[:,i],data[:,j]
        a,b=a[~np.isnan(a)],b[~np.isnan(a)]
        a,b=a[~np.isnan(b)],b[~np.isnan(b)]
        #plt.hist2d(a,b, bins=40, norm=colors.LogNorm())
        plt.hist2d(a,b,weights=np.ones((a.shape[0]))*1000/nb_int, bins=20,cmap='viridis',cmin=0.05,vmin=0)
        plt.colorbar(label="fraction of watershed area (‰)")
        plt.xlabel(L_name_ctxt[i])
        plt.ylabel(L_name_ctxt[j])
        coef=np.corrcoef(a,b)[0,1]
        if j==6:
            Lcoef_perm.append(np.round(coef,2))
        a1,b1=np.polyfit(a,b,1)
        Mcoef[i,j]=coef
        plt.title("r: "+str(np.around(coef,3))+" a: "+str(np.around(a1,3))+" b: "+str(np.around(b1,3)))
        plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/"+L_name_ctxt3[i]+"_X_"+L_name_ctxt3[j]+".png",bbox_inches="tight")
        plt.show()

plt.figure(figsize=(7,7))
plt.imshow(Mcoef,cmap='seismic',vmin=-1,vmax=1)
plt.yticks(range(7),L_name_ctxt3)
plt.xticks(range(7),L_name_ctxt3)
plt.xticks(rotation=90)
plt.colorbar(label='correlation coeficient R')
for k in range(6):
    plt.text(5.6,k+0.1,str(Lcoef_perm[k]), fontsize=15)
plt.grid()
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/graphR.png",bbox_inches="tight")
plt.show()


M2coef=Mcoef*Mcoef
plt.imshow(M2coef,cmap='plasma',vmin=0,vmax=1)
plt.yticks(range(7),L_name_ctxt3)
plt.xticks(range(7),L_name_ctxt3)
plt.xticks(rotation=90)
plt.colorbar(label='R^2')
plt.grid()
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/graphR2.png",bbox_inches="tight")
plt.show()

def around(LM):
    return [round(x,4) for x in LM]

print('R_permafrost',around(Mcoef[:,-1]))
print('R^2_permafrost',around(M2coef[:,-1]))

print('R_altitud',around(Mcoef[0,:]))
print('R^2_altitud',around(M2coef[0,:]))
    
''''Relation intéressantes
altitude-slope
altitude-insolation
altitude-NDVI
slope-insolation
slope-NDVI
orientation-insolation
orientation-NDVI
insolation-NDVI
'''
fig,axs=plt.subplots(1,3,figsize=(15,4))
nbb=data[:,i].shape[0]
for i in [0,1,2]:
    axs[i].hist(data[:,i],bins=20,weights=np.array([100/nbb]*nbb))
    axs[i].set_xlabel(L_name_ctxt[i])
    axs[i].set_ylabel("fraction of wtshd area (%)")
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/3var_dem.png",bbox_inches="tight")
plt.show()


#########################################################################
#QUANTIFIER L'AGROUPATION
if False:
    test=[data[:,0],data[:,1]]
    test=[test[k][L_index_int] for k in range(2)]
    L_dist=[]
    k=0
    nn=len(L_index_int)
    nb=0.5*nn*(nn-1)
    for i in range(n_fichier):
        for j in range(i+1,n_fichier):
            L_dist.append(np.sqrt((test[0][i]-test[0][i])**2+(test[1][i]-test[1][i])**2))
            sys.stdout.write('\r k= {} %'.format(round(100*k/nb,3)) )                         
            sys.stdout.flush()
            k+=1

#Puis calculer par exemple: 


##############################################################################
###CORRELATION SPATIALES
##carte de corrélations spatiales

n,p=133,238
d=20
Lfact=[[-1,0],[-1,1],[0,1],[1,1]]
L_nom_dir=['0°N','45°NE','0°E','45°SE']
ndir=len(Lfact)
MCS=np.zeros((m_fichier,m_fichier,n-(d-1)*2,p-(d-1)*2,d,ndir))

for a in range(m_fichier):
    for b in range(a,m_fichier):
        LX=[data[:,a],data[:,b]]
        LX=[(x-np.nanmean(x))/np.nanstd(x) for x in LX]
        LX=[np.reshape(x,(n,p)) for x in LX]

        for dire in range(ndir):
            fact=Lfact[dire]
            for i in range(d-1,n-d+1):
                for j in range(d-1,p-d+1):
                    for k in range(d):          
                        MCS[a,b,i-(d-1),j-(d-1),k,dire]=abs(LX[0][i+k*fact[0],j+k*fact[1]]-LX[1][i,j])


for a in range(m_fichier):
    for b in range(a,m_fichier):        
        if a!=b:
            plt.imshow(MCS[a,b,:,:,0,0],cmap='gray')
            plt.colorbar()
            plt.title(L_name_ctxt[a]+" et "+L_name_ctxt[b])
            plt.show()
 

####################  
 
MCS2=np.reshape(MCS,(m_fichier,m_fichier,n-(d-1)*2*p-(d-1)*2,d,ndir))
variog=0.5*np.nanmean(MCS2*MCS2,axis=2)
        

for a in range(m_fichier):
    for b in range(a+1,m_fichier):
        for dire in range(ndir):
            plt.plot(range(d),variog[a,b,:,dire], label=L_nom_dir[dire])
        plt.title(L_name_ctxt[a]+" et "+L_name_ctxt[b])
        plt.legend()
        plt.grid()
        plt.show()

#variogramme pour chaque variable seule
L_name_ctxt2=['altitude','slope','orientation','terrain_rugosity_index','insolation_sg22','NDVI_2022','index_permafrost_field']

fig,axs=plt.subplots(4,2, figsize=(10,13))
for a in range(m_fichier):
    for dire in range(ndir):
        axs[a//2,a%2].plot(range(d),variog[a,a,:,dire], label=L_nom_dir[dire])
    axs[a//2,a%2].title.set_text(L_name_ctxt2[a]+" X2")
    axs[a//2,a%2].legend()
    axs[a//2,a%2].grid()
    color='k'
    if a<5:
        color='w'
    axs[a//2,a%2].set_xticks(range(0,d,5),[str(k*30) for k in range(0,d,5)],color=color)
    axs[a//2,a%2].set_xlabel('distance (m)',color=color)
    if a%2==0:
        axs[a//2,a%2].set_ylabel('semivariance (A.U.)')
fig.delaxes(axs[3,1])
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/variog.png",bbox_inches="tight")
plt.show()









#############################################################################
###############################################################################
###############################################################################
#compare index permafrost et model Carla
list1=[]
path2='C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/permafrost_model/'
fp1 = path2+'output_PF_valley_aligned_dcp.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M[L_index_ext]=np.nan
M[M<0]=0
if False:
    M[L_TRI_bedrock]=np.nan
list1.append(M)
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'index permafrost model valley','permafrost_model',M,'copper')

    
    
############################################################################
path='C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/data_dem_index/'
fp1 = path+'output_index_permafrost.csv'

s1=pd.read_csv(fp1,usecols=[2])
M=np.array(s1)
M=M.astype('float64')

M[L_index_ext]=np.nan
if False:
    M[L_TRI_bedrock]=np.nan
list1.append(M)
M=np.reshape(M,(133,238))

plot_map(step2,Llat_long,'index permafrost field','permafrost',M,'copper')


################################################################
from matplotlib import colors

list2=[np.reshape(x,133*238) for x in list1]
Lsum=[np.nansum(x) for x in list2]
print(Lsum)
list2=[x[~np.isnan(x)] for x in list2]
plt.hist2d(list2[0],list2[1], bins=40, norm=colors.LogNorm())
#plt.hist2d(a,b, bins=40, norm=colors.LogNorm())
#plt.hist2d(a,b,weights=np.ones((a.shape[0]))*1000/nb_int, bins=20,cmap='viridis',cmin=0.05,vmin=0)
plt.colorbar(label="portion of total watershed area (‰)")
plt.xlabel('model')
plt.ylabel('field')
plt.show()



#correlation spatiales...

n,p=133,238
d=20
Lfact=[[-1,0],[-1,1],[0,1],[1,1]]
L_nom_dir=['0°N','45°NE','0°E','45°SE']
ndir=len(Lfact)
MCS=np.zeros((n-(d-1)*2,p-(d-1)*2,d,ndir))


LX=[list1[0],list1[1]]
#LX=[(x-np.nanmean(x))/np.nanstd(x) for x in LX] #centrer-réduire
LX=[np.reshape(x,(n,p)) for x in LX]

for dire in range(ndir):
    fact=Lfact[dire]
    for i in range(d-1,n-d+1):
        for j in range(d-1,p-d+1):
            for k in range(d):          
                MCS[i-(d-1),j-(d-1),k,dire]=abs(LX[0][i+k*fact[0],j+k*fact[1]]-LX[1][i,j])


plt.imshow(MCS[:,:,0,0],cmap='gray')
plt.colorbar()
#plt.title(L_name_ctxt[a]+" et "+L_name_ctxt[b])
plt.show()


####################  
#variogramme field VS model
MCS2=np.reshape(MCS,(n-(d-1)*2*p-(d-1)*2,d,ndir))
variog=0.5*np.nanmean(MCS2*MCS2,axis=0)
        

for dire in range(ndir):
    plt.plot(range(d),variog[:,dire], label=L_nom_dir[dire])
#plt.title(L_name_ctxt[a]+" et "+L_name_ctxt[b])
plt.legend()
plt.xticks(range(0,d,5),[str(k*30) for k in range(0,d,5)],color=color)
plt.xlabel('distance (m)',color=color)
plt.ylabel('semivariance (A.U.)')
plt.grid()
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/variog_permafrost.png",bbox_inches="tight")
plt.show()




#########################################"

d=1
MCS_simple=np.zeros((n-(d-1)*2,p-(d-1)*2))
for i in range(d-1,n-d+1):
    for j in range(d-1,p-d+1):         
            MCS_simple[i-(d-1),j-(d-1)]=LX[0][i,j]-LX[1][i,j]



def grid_ok():
    size=16
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)
    plt.rc('axes', axisbelow=True)
    plt.figure(figsize=(17,7.7))
    plt.grid(zorder=2)
    plt.plot([184,223],[111,111],linewidth=7.0,color='k')
    plt.text(195,105,'1 km', fontsize=size)
    plt.xticks(range(0,238,step2),Llat_long[0])
    plt.yticks(range(0,133,step2),reversed(Llat_long[1]))
    plt.xlabel('longitude (°)')
    plt.ylabel('latitude (°)')
    

from pylab import *#PiYG
grid_ok()
cmap=cm.get_cmap('Spectral', 8)
plt.imshow(MCS_simple,cmap=cmap,vmin=-1,vmax=1)
plt.colorbar(label="index model minus index field",ticks=[-1+0.25*k for k in range(9)])
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/diff_permafrost.png",bbox_inches="tight")
plt.show()

from pylab import *#PiYG
grid_ok()
cmap=cm.get_cmap('Spectral', 3)
plt.imshow(MCS_simple,cmap=cmap,vmin=-1.05,vmax=1.05)
plt.colorbar(label="index model minus index field",ticks=[-1.05,-0.35, 0.35,1.05])
#plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/diff_permafrost.png",bbox_inches="tight")
plt.show()

##################################




grid_ok()
plt.imshow(abs(MCS_simple),cmap='BuGn')
plt.colorbar(label="absolute value of index model minus index field")
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/abs_diff_permafrost.png",bbox_inches="tight")
plt.show()











