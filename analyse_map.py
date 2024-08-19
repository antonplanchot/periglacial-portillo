import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Lfichier=["output_DEM_ALOS_dcp","output_slope_ALOS_dcp","output_orientation_ALOS_dcp","output_TRI_ALOS_dcp"]

def settxtsize(size):
    size=size
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)

#importation cartographie final
map=pd.read_csv("C://Users//Utilisateur//Documents//_Internship_Anton//cartography_Portillo//treatment//attribut_intersection.csv",sep=',',usecols=[1,12,13,14])
map1=np.array(map)
nmap,mmap=map1.shape


#importation DEM ALOS et cie
Ldem=[]
for k in range(4):
    dem1=pd.read_csv('C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/data_dem_index/'+Lfichier[k]+'.csv',sep=',')
    Ldem.append(np.array(dem1))


ngrid,mgeof=31654,27  

#construction d'une matrice des aires de chaque geoforme pour chaque cellule
M=np.zeros((ngrid,mgeof))
for k in range(nmap):
    geoform=map1[k,0]
    if geoform>0:
        i=int(map1[k,1]*238+map1[k,2])
        j=int(geoform-1)
        M[i,j]+=map1[k,3]

Lsum=[np.sum(M[k,:]) for k in range(ngrid)]#aire de chaque cellule du watershed
Lsum_geof=[np.sum(M[:,k]) for k in range(mgeof)]#aire totale de chaque type de geoforme
M=np.c_[M,np.array(Lsum)]

L_index_exterior=[]
for k in range(ngrid):
    if Lsum[k]==0:
        L_index_exterior.append(k)

###########################################################################
Lmax=[]
for k in range(M.shape[0]):
    L=list(M[k,:27])
    ind=L.index(max(L))
    Lmax.append(ind)

Mmax=np.array([Lmax])
Mmax1=np.reshape(Mmax,(133,238))
plt.imshow(Mmax1)
plt.colorbar()
plt.show()


###########################################################################        
Lgeof_permafrost=[14,15,16,17,20,21,22]
Lgeof_maybe_permafrost=[13,23,24]
Lgeof_maybe2_permafrost=[9,12]
############################################################

Lcompt=[0]*3
LLgeof=[Lgeof_permafrost,Lgeof_maybe_permafrost,Lgeof_maybe2_permafrost]
for x in Lmax:
    for k in range(3):
        if x in LLgeof[k]:
            Lcompt[k]+=1
print(len(Lmax)-np.sum(Lcompt))

###############################################################"

Mpermafrost=np.sum(M[:,Lgeof_permafrost],axis=1)
Mmaybe_permafrost=np.sum(M[:,Lgeof_maybe_permafrost],axis=1)
Mmaybe2_permafrost=np.sum(M[:,Lgeof_maybe2_permafrost],axis=1)


##test des grosses valeurs anormales
if True:
    b=list(Mpermafrost) 
    for k in range(100):
        a=np.max(b)  
        print(a)
        b.remove(a)
    
    mp=Mpermafrost[Mpermafrost>1]
    plt.hist(mp,bins=3)
    plt.show()
    
    mmp=Mmaybe_permafrost[Mmaybe_permafrost>1]
    plt.hist(mmp,bins=3)
    plt.show()
    
    mm2p=Mmaybe_permafrost[Mmaybe2_permafrost>1]
    plt.hist(mm2p,bins=3)
    plt.show()


#indice parmafrost
seuil=956.1750847546101
seuil2=956.17508475461
seuil3=956.2
Mpermafrost[Mpermafrost==0]=0.01
Mpermafrost[Mpermafrost>seuil]=seuil2 #retirer valeur aberrantes
Mind_permafrost=Mpermafrost%seuil/seuil3+0.5*Mmaybe_permafrost%seuil/seuil3+0.25*Mmaybe2_permafrost%seuil/seuil3

Mind_permafrost=np.reshape(Mind_permafrost,(133,238))

plt.imshow(Mind_permafrost)
plt.colorbar()
plt.show()

path='C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/data_dem_index/'
a=np.concatenate((Ldem[0][:,:2],np.reshape(Mind_permafrost,(133*238,1))),axis=1)
np.savetxt(path+'output_index_permafrost.csv', a, delimiter=',',header='x,y,slope', comments='')



################################################################################
#camembert
name_geoform=['bedrock_massive','bedrock_thin','debris_undifferentiated','debris_slope','debris_slope_alim_transfert','debris_cone','debris_cone_non_recharging','landslide_debris_slope_cone','debris_mud_flow_cone','nival_avalanche_debris_cone','rock_avalanche','moraine','frozen_debris_slope_small','frozen_debris_slope_big','cryoplain','protalus_rampart','protalus_lobe','embryonic_rock_glacier','embryonic_rock_glacier_fossil','fossil_rock_glacier','transitioning_rock_glacier','active_rock_glacier','glacier_covered_with_rock_glacier','covered_glacier','glaciaret','current_fluvial_stream','vegas','watershed']
name_geoform_acro=['BM','BTh','DU','DS','DSAT','DC','DCnR','LDSC','DMFC','NADC','RAv','M','FrDSS','FrDSB','Cp','PR','PL','ERG','ERGF','RGF','RGT','RGA','CGwRG','CG','Gla','CFlS','Vg','WTSHD']
Lsum_geof_tri=list(np.copy(Lsum_geof))
name_geoform_tri=list(np.copy(name_geoform))[:-1]
if False:
    for _ in range(mgeof-1):
        for k in range(mgeof-1):
            if Lsum_geof_tri[k]<Lsum_geof_tri[k+1]:
                Lsum_geof_tri[k],Lsum_geof_tri[k+1]=Lsum_geof_tri[k+1],Lsum_geof_tri[k]
                name_geoform_tri[k],name_geoform_tri[k+1]=name_geoform_tri[k+1],name_geoform_tri[k]
                
    plt.figure(figsize=(10,10))
    
    labels = name_geoform_tri[:19]+['other']
    size=Lsum_geof_tri[:19]+[np.sum(Lsum_geof_tri[19:])]
else:
    plt.figure(figsize=(10,10))
    
    labels = name_geoform_tri
    size=Lsum_geof_tri

# Création du graphique circulaire
plt.pie(size, labels=labels, autopct="%1.1f%%")
plt.title("Proportion de chaque géoformes")
plt.show()
#idée triée en non périglaciaire // périglaciaire



#datas ING Tunuyan Norte
#Tipo de geoforma inventariada Área (km2)
L_ING=[['Uncovered \n glacier','Glaciarets','Covered glacier','Covered \n glacier with \n rock glacier','Rock glacier \n active and transitioning','Total'],[38.22,7.17,2.30,45.96,38.67,132.32]]
colors=['deepskyblue', 'cyan','sandybrown','peru', 'saddlebrown']

print(np.sum(Lsum_geof_tri[20:25])/(100*10**4))
print((1772047.611+56283.552+221032.402+54513.199+80575.375+84976.758+40886.894+186103.019+224658.45)/(100*10**4))

##icy geoforms
fig,axs=plt.subplots(1,2,figsize=(10,10))
settxtsize(10)
axs[0].pie(L_ING[1][:-1], labels=L_ING[0][:-1], autopct="%1.1f%%", colors=colors,radius=1,pctdistance=0.8,labeldistance=1.25)
axs[0].set_title('Watershed river Tunuyan norte \n Total area with icy geoforms: 132,32 km2', fontweight="bold")
axs[1].pie([Lsum_geof_tri[k] for k in [24,23,22]]+[np.sum(Lsum_geof_tri[20:22])], labels=L_ING[0][1:-1], autopct="%1.1f%%", colors=colors[1:],radius=1,pctdistance=0.8,labeldistance=1.2)
axs[1].set_title('Watershed river Portillo \n ING: 2,72 km2 / this work: 2,65 km2 (2.5% less)', fontweight="bold")
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/icy_geoforms.png",bbox_inches="tight")
plt.show()


##periglacial and glacial geoforms
colors2=['deepskyblue', 'cyan','sandybrown','peru', 'saddlebrown','sienna','darkgoldenrod','gold','orange','orangered','indianred','brown']
L_name_forms=['Glaciarets','Covered glacier','Covered \n glacier with \n rock glacier','Rock glacier \n active','Rock glacier \n transitioning','Embryonic \n rock glacier','Protalus lobe','Protalus rempart','Cryoplain','Frozen debris \n slope large','Frozen debris \n slope small']
plt.figure(figsize=(6,6))
settxtsize(10)
L_periglacial=[24,23,22,21,20,17,16,15,14,13,12]
plt.pie([Lsum_geof_tri[k] for k in L_periglacial], autopct="%1.1f%%", 
        labels=L_name_forms,colors=colors2[1:],radius=0.9,pctdistance=1.12,labeldistance=1.31) #, labels=[name_geoform_acro[k] for k in L_periglacial]
#draw circle
centre_circle = plt.Circle((0,0),0.42,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
L_perig=[np.sum(Lsum_geof_tri[20:25]),np.sum(Lsum_geof_tri[12:18])]
plt.pie(L_perig, autopct="%1.1f%%", labels=["ING \n geoforms",""],colors=['darkgrey','whitesmoke'],radius=0.42,pctdistance=0.52,labeldistance=0.2)
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/geoforms_periglacial_glacial.png",bbox_inches="tight")
plt.show()


##not glacial and periglacial deposits
colors3=['fuchsia','orchid','salmon','yellow','darkcyan','darkslategrey','cyan','lightcyan','k','silver','yellowgreen','olive','w']
L_name_forms2=['Bedrock massive','Bedrock sporadic','Debris slope','Debris slope \n alimentation transfert \n (and Landslide debris \n slope cone: 0.1%)','Debris cone','Debris cone unsupplied','Debris mud flow cone','Nival avalanche debris cone','Rock avalanche','Debris undifferentiated','Moraine','Fossil developed \n and embryonic \n rock glacier','Current riverbed \n and vegas']

fig,axs=plt.subplots(1,2,figsize=(12,10))
settxtsize(10)
L_other=[0,1,3,4,5,6,8,9,10,2,11]
L_areas=[Lsum_geof_tri[k] for k in L_other[:3]]+[Lsum_geof_tri[3]+Lsum_geof_tri[7]]+[Lsum_geof_tri[k] for k in L_other[4:]]+[Lsum_geof_tri[18]+Lsum_geof_tri[19]]+[Lsum_geof_tri[25]+Lsum_geof_tri[26]]
axs[0].pie(L_areas, autopct="%1.1f%%", colors=colors3,labels=L_name_forms2, 
        radius=1,pctdistance=1.1,labeldistance=1.25) 
L_areas2=[np.sum(L_areas[:2]),np.sum(L_areas[2:4]),np.sum(L_areas[4:8]),np.sum(L_areas[8:13])]
axs[0].pie(L_areas2,colors=['darkgrey','whitesmoke','darkgrey','whitesmoke','darkgrey','whitesmoke'], autopct="%1.1f%%",radius=0.45,pctdistance=0.5,labeldistance=0.2)
axs[0].set_title('Not glacial or periglacial deposits', fontweight="bold")
axs[1].set_title('Summary', fontweight="bold")
axs[1].pie(L_perig+L_areas2,colors=['darkgrey','whitesmoke','darkgrey','whitesmoke','darkgrey','whitesmoke'], autopct="%1.1f%%", 
        labels=["ING \n geoforms","Other \n periglacial geoforms", "Bedrock","Debris slope \n total","Debris cone total","Other"],
        radius=0.9,pctdistance=0.8,labeldistance=1.1) #, labels=[name_geoform_acro[k] for k in L_periglacial]
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/geoforms_other.png",bbox_inches="tight")
plt.show()




# Différence explicable par la relative faible altitude de notre vallée


###############################################################################
#plot et save bedrock
Mbedrock=np.reshape(M[:,0],(238,133),order='F').T

seuil=956.1750847546101
seuil2=956.17508475461
seuil3=956.2
Mbedrock[Mbedrock>seuil]=seuil2 #retirer valeur aberrantes
Mbedrock=Mbedrock%seuil/seuil3
       
plt.imshow(Mbedrock)
plt.colorbar(label="area (m2)")
plt.title('bedrock')
plt.show()

path='C:/Users/Utilisateur/Documents/_Internship_Anton/cartography_Portillo/data_dem_index/'
a=np.concatenate((Ldem[0][:,:2],np.reshape(Mbedrock,(133*238,1))),axis=1)
np.savetxt(path+'output_bedrock.csv', a, delimiter=',',header='x,y,slope', comments='')

 

#plot des cartes de présence des géoformes pour chaque pixel ALOS
for k in [3,5,20,21,22]:
    Mgeoform_map=np.reshape(M[:,k],(238,133),order='F').T
    plt.imshow(Mgeoform_map)
    plt.colorbar(label="area (m2)")
    plt.title(name_geoform[k])
    plt.show()

#plot sans somme de chaque occurence de bedrock
if False:
    plt.plot(Ldem[0][:,2],M[:,0],",")
    plt.show()





###############################################################################
###############################################################################
## ANALYSE avec data DEM et cie
pas_slope=0.01
len_variable=[int(max(Ldem[0][:,2])-min(Ldem[0][:,2]))+1,int(max(Ldem[1][:,2])//pas_slope)+1,360,int(max(Ldem[3][:,2])-min(Ldem[3][:,2]))+1]
Mdensite=[np.zeros((x,mgeof+1)) for x in len_variable]

for i in range(mgeof+1):
    for j in range(ngrid):
        L_index=[Ldem[0][j,2]-3223,Ldem[1][j,2]//pas_slope,Ldem[2][j,2]//1,Ldem[3][j,2]-2]
        if not np.isnan(L_index[2]):
            L_index=[int(x) for x in L_index]
            LL=range(4)
        else:
            L_index[2]=False
            L_index=[int(x) for x in L_index]
            LL=[0,1,3]
        for k in LL:
            Mdensite[k][L_index[k],i]+=M[j,i]


Mintegral=[np.zeros((x,mgeof+1)) for x in len_variable]
for i in range(mgeof+1):
    for k in range(4):
        S=0
        for j in range(len_variable[k]):
            S+=Mdensite[k][j,i]
            Mintegral[k][j,i]=S
        Mintegral[k][:,i]=100*Mintegral[k][:,i]/S
        
  
if not np.isnan(np.nan):
    print(1)
################################################################
###############################################################
#Diagramme moustache !!!


L_range_variable=[range(3223,4830+1), [pas_slope*i for i in range(len_variable[1])],range(360),range(2,162)]


Mdensite_stat=[np.zeros((2,mgeof+1)) for x in len_variable]
Mmodelage_violin1=[]
Mmodelage_violin2=[]
Mmodelage_violin3=[]
#moyenne, médiane, std, Q1, Q2, min, max
Lsum2=Lsum_geof+[np.sum(Lsum_geof)]
nb_stat=len(Mdensite_stat)
for k in range(nb_stat):
    for i in range(mgeof+1):
        #violinplot preparation
        if k==0:
            L=[]
            j=3223
            for x in Mdensite[0][:,i]:
                L+=[j]*int(x)
                j+=1
            Mmodelage_violin1.append(L)
        if k==1:
            L=[]
            j=0
            for x in Mdensite[1][:,i]:
                L+=[j]*int(x)
                j+=0.01
            Mmodelage_violin2.append(L)
        if k==2:
            L=[]
            j=0
            for x in Mdensite[2][:,i]: 
                L+=[j-360*(j>250)]*int(x)
                j+=1
            Mmodelage_violin3.append(L)
        #statistiques autres
        mw=np.sum([x*y for x,y in zip(Mdensite[k][:,i],L_range_variable[k])])/Lsum2[i]
        stdw=np.sqrt(np.sum([x*(y-mw)**2 for x,y in zip(Mdensite[k][:,i],L_range_variable[k])])/Lsum2[i])
        Mdensite_stat[k][:,i]=np.array([mw,stdw]).T

if False:
    for k in range(nb_stat):
        plt.figure(figsize=(10,3))
        plt.plot(range(mgeof+1),Mdensite_stat[k][0,:],'*')
        plt.xticks(range(mgeof+1),name_geoform_acro)
        plt.xticks(rotation=45)
        plt.show     

def dem2maat(ddem):
    z_station=3950
    maat_station=-1.285
    grad=-0.006
    return (ddem-z_station)*grad+maat_station

def maat2dem(mmaat):
    z_station=3950
    maat_station=-1.285
    grad=-0.006
    return (mmaat-maat_station)/grad+z_station

L_index_keep=[24,23,22,21,20,19,18,17,16,15,14,13,12,0,1,3,4,5,6,8,9,10,2,11,25,26,27]
L_name_forms_violin=['Glaciarets','Covered glacier','Covered glacier \n with ock glacier','Rock glacier active','Rock glacier \n transitioning','Rock glacier fosil',
                     'Embryonic \n rock glacier fossil','Embryonic RG','Protalus lobe','Protalus rempart','Cryoplain','Frozen debris \n slope large','Frozen debris \n slope small',
                     'Bedrock massive','Bedrock sporadic','Debris slope','Debris slope \n alimentation transfert','Debris cone','Debris cone unsupplied','Debris mud flow cone',
                     'Nival avalanche \n debris cone','Rock avalanche','Debris undifferentiated','Moraine','Current riverbed','Vegas','Watershed']


##violin for the altitude / maat
hauteur=6


ymin,ymax=3300,4850
med_w=np.median(Mmodelage_violin1[27])

fig, ax1 = plt.subplots(figsize=(15,hauteur))
settxtsize(14)
ax1.axhline(y=med_w,color='orange',linestyle='dashed')
ax1.axhline(y=maat2dem(0),color='red',linestyle='dashed')
vp=ax1.violinplot([Mmodelage_violin1[i] for i in L_index_keep], widths=0.65, showmeans=False,showextrema=False, points=15)
colors_violin=['blue']+['saddlebrown']*6+['orange']*4+['brown']*2+['fuchsia']*2+['salmon']*2+['darkslategrey']*4+['olive']*3+['green']*2+['grey']
for pc, color in zip(vp['bodies'], colors_violin):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
ax1.boxplot([Mmodelage_violin1[i] for i in L_index_keep],showfliers=False)
ax1.set_ylim(ymin,ymax)
ax1.set_ylabel('Altitude (m)')
plt.xticks(range(1,mgeof+1),L_name_forms_violin)
plt.xticks(rotation=90)
plt.xticks(color='w')
ax2 = ax1.twinx()
ax2.set_ylim(dem2maat(ymin),dem2maat(ymax))
ax2.set_ylabel('Mean annual air temperature (°C)')
ax2.set_yticks(list(range(2,-8,-1)))
ax2.grid()
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/violin_maat.png",bbox_inches="tight")
plt.show

##violin for the slope/TIR

def slope2tir(slope):
    a=74.069
    b=3.247
    return a*slope+b

ymin,ymax=-0.05,1.3
med_w=np.median(Mmodelage_violin2[27])

fig, ax1 = plt.subplots(figsize=(15,hauteur))
settxtsize(14)
ax1.axhline(y=med_w,color='orange',linestyle='dashed')
vp=ax1.violinplot([Mmodelage_violin2[i] for i in L_index_keep], widths=0.8, showmeans=False,showextrema=False, points=15)
colors_violin=['blue']+['saddlebrown']*6+['orange']*4+['brown']*2+['fuchsia']*2+['salmon']*2+['darkslategrey']*4+['olive']*3+['green']*2+['grey']
for pc, color in zip(vp['bodies'], colors_violin):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
ax1.boxplot([Mmodelage_violin2[i] for i in L_index_keep],showfliers=False)
ax1.set_ylim(ymin,ymax)
ax1.set_ylabel('Slope (m/m)')
ax1.grid()
plt.xticks(range(1,mgeof+1),L_name_forms_violin)
plt.xticks(rotation=90)
plt.xticks(color='w')
ax2 = ax1.twinx()
ax2.set_ylim(slope2tir(ymin),slope2tir(ymax))
ax2.set_ylabel('Terrain rugosity index')
#ax2.set_yticks(list(range(2,-8,-1)))
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/violin_slope.png",bbox_inches="tight")
plt.show


##essai violin orientation

ymin,ymax=-85,225
med_w=np.median(Mmodelage_violin3[27])

fig, ax1 = plt.subplots(figsize=(15,hauteur))
settxtsize(14)
#ax1.axhline(y=med_w,color='orange',linestyle='dashed')
ax1.boxplot([Mmodelage_violin3[i] for i in L_index_keep],showfliers=False)
vp=ax1.violinplot([Mmodelage_violin3[i] for i in L_index_keep], widths=0.8, showmeans=False,showextrema=False, points=25)
colors_violin=['blue']+['saddlebrown']*6+['orange']*4+['brown']*2+['fuchsia']*2+['salmon']*2+['darkslategrey']*4+['olive']*3+['green']*2+['grey']
for pc, color in zip(vp['bodies'], colors_violin):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
#ax1.boxplot([Mmodelage_violin3[i] for i in L_index_keep],showfliers=False)
ax1.set_ylim(ymin,ymax)
ax1.set_yticks(list(range(-45,226,45)))
ax1.set_ylabel('Orientation (°)')
ax1.grid()
plt.xticks(range(1,mgeof+1),L_name_forms_violin)
plt.xticks(rotation=90)
ax2 = ax1.twinx()
ax2.set_ylim(ymin,ymax)
#ax2.set_ylabel('Terrain rugosity index')
ax2.set_yticks(range(-45,226,45),['North-West','North','North-Est','Est','South-Est','South','South-West'])
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/violin_orientation.png",bbox_inches="tight")
plt.show


##orientation normalisée
Morient=np.copy(Mdensite[2])
Sorient=np.sum(Morient,axis=0)
for j in range(28):
    Morient[:,j]=Morient[:,j]/Sorient[j]
for i in range(360):
    if Morient[i,27]!=0:
        Morient[i,:-1]=Morient[i,:-1]/Morient[i,27]


Mmodelage_violin4=[]
for i in range(mgeof+1):
    L=[]
    j=0
    for x in Morient[:,i]: 
        L+=[j-360*(j>250)]*int(x*100)
        j+=1
    Mmodelage_violin4.append(L)


#plot
ymin,ymax=-90,250
med_w=np.median(Mmodelage_violin3[27])

fig, ax1 = plt.subplots(figsize=(15,hauteur))
settxtsize(14)
#ax1.axhline(y=med_w,color='orange',linestyle='dashed')
ax1.boxplot([Mmodelage_violin4[i] for i in L_index_keep],showfliers=False)
vp=ax1.violinplot([Mmodelage_violin4[i] for i in L_index_keep][:-1], widths=0.8, showmeans=False,showextrema=False, points=25)
colors_violin=['blue']+['saddlebrown']*6+['orange']*4+['brown']*2+['fuchsia']*2+['salmon']*2+['darkslategrey']*4+['olive']*3+['green']*2
for pc, color in zip(vp['bodies'], colors_violin):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
#ax1.boxplot([Mmodelage_violin3[i] for i in L_index_keep],showfliers=False)
ax1.set_ylim(ymin,ymax)
ax1.set_yticks(list(range(-90,226,45)))
ax1.set_ylabel('Orientation (°)')
ax1.grid()
#plt.xticks(range(1,mgeof),L_name_forms_violin[:-1])
plt.xticks(range(1,mgeof+1),L_name_forms_violin)
plt.xticks(rotation=90)
ax2 = ax1.twinx()
ax2.set_ylim(ymin,ymax)
#ax2.set_ylabel('Terrain rugosity index')
ax2.set_yticks(range(-90,226,45),['West','North-West','North','North-Est','Est','South-Est','South','South-West'])
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/violin_orientation_normal.png",bbox_inches="tight")
plt.show











##############################################################
###############################################################
###############################################################
##ex
if False:
    A=[[500,500,0,0],[250,250,250,250]]
    for L in A:
        S=np.nansum([x*4000 for x in L])/np.nansum(L)
        print(S)

##lissage
int_lissage=[100,10,30,8]
Mdensite_liss=[np.zeros((len_variable[k]+(1-int_lissage[k])*(k!=2),mgeof+1)) for k in range(4)]
for k in range(4):
    if k!=2:
        for i in range(mgeof+1):
            for j in range(len_variable[k]-int_lissage[k]):
                Mdensite_liss[k][j,i]=np.mean(Mdensite[k][j:j+int_lissage[k],i])
    else:
        Mdensite_orientation=np.concatenate((Mdensite[k],Mdensite[k][:int_lissage[k],:]),axis=0)
        for i in range(mgeof+1):
            for j in range(len_variable[k]):
                Mdensite_liss[k][j,i]=np.mean(Mdensite_orientation[j:j+int_lissage[k],i])




#plot de ces densités
L_range_variable=[range(3223,4830+1),range(0,pas_slope*len_variable[1],pas_slope),range(360),range(150)]
L_legend=['altitud (m)', 'slope', 'orientation (°)', 'rugosity']

if True:
    Lplot=[0,5,17,27]
    for x in range(4):
        for k in Lplot:
            plt.figure()
            plt.plot(L_range_variable[x],Mdensite[x][:,k],'.')
            if x!=2:
                plt.plot(L_range_variable[x][int_lissage[x]//2:-int_lissage[x]//2+1],Mdensite_liss[x][:,k])
            else:
                plt.plot(L_range_variable[x],Mdensite_liss[x][:,k])
            plt.title(name_geoform[k]+" // "+Lfichier[x])
            plt.xlabel(L_legend[x])
            plt.ylabel("area (m2)")
            plt.show()

#plot intégrals
Lplot=[0,3,4,5,8,9,27]
for x in range(4):
    plt.figure(figsize=(10,10))
    for k in Lplot:
        plt.plot(L_range_variable[x],Mintegral[x][:,k],label=name_geoform[k])
        plt.xlabel(L_legend[x])
        plt.ylabel("fraction area (%)")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(Lfichier[x])
    plt.show() 






####plot orientation

# positions angulaires
angles = 2*np.pi*np.linspace(0, 1, 9)
angles2 = 2*np.pi*np.linspace(0, 1, 360) # 360èmes de cercle
radDeg = 180/np.pi # conversion radians → degrés

geof1 = Mdensite_liss[2][:,0]
geof2 = Mdensite_liss[2][:,27]
geof3 = Mdensite_liss[2][:,17]
geof1 = geof1/geof1.max() # normalisation
geof2 = geof2/geof2.max()
geof3 = geof3/geof3.max()

# ***********
# * Traçage *
# ***********

courbes = plt.polar(angles2, geof1, angles2, geof2, angles2, geof3) 

# mise en forme
axes = plt.gca()
courbes[0].set(color="r", label="bedrock")
courbes[1].set(color="k", label="watershed")
courbes[2].set(color="g", label="embryonic_RG")
plt.legend(bbox_to_anchor=(1.45, 1.05), loc="upper right") # affichage à l'extérieur du tracé
axes.set_rticks([0.2, 0.4, 0.6, 0.8, 1]) # graduations radiales
axes.set_thetagrids(angles=radDeg*angles[0:-1], labels=("E", "NE", "N", "NO", "O", "SO", "S", "SE")) # graduations angulaires
plt.title("Orientation des terrains", va="bottom")

plt.show()















#########################################################################
#########################################################################
#importation MAP2, INDEX NDVI
Lpleiade=[]
Lfichier2=['Pleiade/output_test_Pleiade_NDVI_2022_projec_dcp','cartography_Portillo/treatment/output_portillo_raster_Pleiade_plusmarge']
for k in range(len(Lfichier2)):
    index1=pd.read_csv('C:/Users/Utilisateur/Documents/_Internship_Anton/'+Lfichier2[k]+'.csv',sep=',')
    Lpleiade.append(np.array(index1)[:,2])
threshold_shadow=0.025
#1864*3353
#+8+5=1877//+5=3358 ##vue sur raster_Pleiade
#1877*3359 et 1877*3358
Lpleiade[0]=np.reshape(Lpleiade[0],(1877,3359))
Lpleiade[0]=Lpleiade[0][:,1:]
Lpleiade[1]=np.reshape(Lpleiade[1],(1877,3358))

plt.imshow(Lpleiade[0])
plt.show()

plt.imshow(Lpleiade[1])
plt.show()

Lpleiade[0]=np.reshape(Lpleiade[0],(1877*3358))
Lpleiade[0][Lpleiade[0]<threshold_shadow]=np.nan
Lpleiade[1]=np.reshape(Lpleiade[1],(1877*3358))


Lndvi_geof=[Lpleiade[0][Lpleiade[1]==k] for k in range(1,28)]+[Lpleiade[0]]
Lndvi_geof=[Lndvi_geof[k][Lndvi_geof[k]>=threshold_shadow] for k in range(28)]


mgeof=27
if False:
    name_geoform_acro=['BM','BTh','DU','DS','DSAT','DC','DCnR','LDSC','DMFC','NADC','RAv','M','FrDSS','FrDSB','Cp','PR','PL','ERG','ERGF','RGF','RGT','RGA','CGwRG','CG','Gla','CFlS','Vg','WTSHD']
    plt.figure(figsize=(15,5))
    plt.violinplot(Lndvi_geof[1:], showmeans=True, showextrema=True, points=12)
    plt.xticks(range(1,mgeof+1),name_geoform_acro[:-1])
    plt.xticks(rotation=45)
    plt.grid()
    plt.show
    
    plt.figure(figsize=(15,5))
    plt.boxplot(Lndvi_geof[1:],showfliers=False)
    plt.xticks(range(1,mgeof+1),name_geoform_acro[:-1])
    plt.xticks(rotation=45)
    plt.grid()
    plt.show


'''The box extends from the first quartile (Q1) to the third quartile (Q3)
 of the data, with a line at the median. The whiskers extend from the box 
 to the farthest data point lying within 1.5x the inter-quartile range 
 (IQR) from the box. Flier points are those past the end of the whiskers.'''


L_index_keep=[24,23,22,21,20,19,18,17,16,15,14,13,12,0,1,3,4,5,6,8,9,10,2,11,25,26,27]
L_name_forms_violin=['Glaciarets','Covered glacier','Covered glacier \n with rock glacier','Rock glacier active','Rock glacier \n transitioning','Rock glacier fossil',
                     'Embryonic \n rock glacier fossil','Embryonic RG','Protalus lobe','Protalus rempart','Cryoplain','Frozen debris \n slope large','Frozen debris \n slope small',
                     'Bedrock massive','Bedrock sporadic','Debris slope','Debris slope \n alimentation transfer','Debris cone','Debris cone unsupplied','Debris mud flow cone',
                     'Nival avalanche \n debris cone','Rock avalanche','Debris undifferentiated','Moraine','Current riverbed','Vegas','Watershed']

if False:
    for st in L_name_forms_violin:
        print(st)

#L_index_keep=L_index_keep[:18]
##violin for ndvi

ymin,ymax=0,0.59
med_w=np.median(Lndvi_geof[27])

fig, ax1 = plt.subplots(figsize=(15,hauteur))
settxtsize(14)
ax1.axhline(y=med_w,color='orange',linestyle='dashed')
vp=ax1.violinplot([Lndvi_geof[i] for i in L_index_keep], showmeans=False,showextrema=False, points=15)
colors_violin=['blue']+['saddlebrown']*6+['orange']*4+['brown']*2+['fuchsia']*2+['salmon']*2+['darkslategrey']*4+['olive']*3+['green']*2+['grey']
for pc, color in zip(vp['bodies'], colors_violin):
    pc.set_facecolor(color)
    pc.set_edgecolor('black')
ax1.boxplot([Lndvi_geof[i] for i in L_index_keep],showfliers=False)
ax1.set_ylim(ymin,ymax)
ax1.set_yticks([i/20 for i in range(13)])
ax1.set_ylabel('NDVI')
plt.xticks(range(1,mgeof+1),L_name_forms_violin)
plt.xticks(rotation=90)
ax2 = ax1.twinx()
ax2.set_ylim(ymin,ymax)
ax2.set_yticks([i/20 for i in range(13)])
ax2.grid()
plt.savefig("C:/Users/Utilisateur/Documents/_Internship_Anton/$report/figures/violin_ndvi.png",bbox_inches="tight")
plt.show





###############################################################################
###############################################################################
###############################################################################
################################################################################





























if False:

    ##test
    
    a=np.array([[1,2,3,3,3,4,1,3,2,2],[1,2,3,3,3,4,1,3,2,2]]).T
    plt.violinplot(a)
    plt.show
    
    
    
    #######################################
    #cas particulier altitude
    nbclass=17
    multip=100
    Mdensite_class100=np.zeros((nbclass,mgeof+1))
    for i in range(mgeof+1):
        for j in range(nbclass):
            Mdensite_class100[j,i]=np.sum(Mdensite[0][j*multip:(j+1)*multip,i])
            
            
            
    L=range(3250,4851,100)          
    if False:        
          
        for k in range(28):
            plt.plot(L,Mdensite_class100[:,k])
            plt.title(name_geoform[k])
            plt.xlabel("altitud (m)")
            plt.ylabel("area (m2)")
            plt.show()
            
        for k in range(28):
            plt.plot(L,Mdensite_class100[:,k],label=name_geoform[k])
            plt.xlabel("altitud (m)")
            plt.ylabel("area (m2)")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show() 

