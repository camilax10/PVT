# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 13:25:32 2021

@author: Camila Espinosa
"""
import numpy
from scipy.linalg import solve
from scipy import linalg
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 12


#CONDUCTIVIDAD TERMICA MATERIALES [W/mK]import matplotlib.animation as animationt 
w=1
k=110
k_glass=1.8
k_eva=0.35
k_pv=148
k_ted=0.2


#DENSIDADES DEL PANEL POR CAPA [kg/m3]
rho=8530
rho_glass=3000
rho_eva=960
rho_pv=2330
rho_ted=1200


#CAPACIDAD DE CALOR ESPECIFICO [J/kg K] 
C=380
C_glass=500
C_eva=2090
C_pv=677
C_ted=1250

#PRODUCTO DE DENSIDAD Y CALOR ESPECIFICO
p=rho*C

p_glass=rho_glass*C_glass
p_eva=rho_eva*C_eva
p_pv=rho_pv*C_pv
p_ted=rho_ted*C_ted


#ESPESORES MATERIALES [mm]
t_glass=float(3)
t_eva=float(0.5)
t_pv=float(0.24)
t_ted=float(0.1)
ancho=t_glass+t_eva+t_ted 
l_silicio=float(828) #cambio de 78 a 786 (828-42) extension panel que es la mitad del panel, se estendio a casi completo 1740
l_ant=float(42)
largo=l_ant+l_silicio

#ESPESORES MATERIALES [m]
tglass=t_glass/1000 
teva=t_eva/1000 
tpv=t_pv/1000 
tted=t_ted/1000

ancho_m=ancho/1000
largo_m=largo/1000


ngxl=29 # cambio panel 29 nodos
ngxpv=87 # cambio nodos longitud 87 nodos
ngx=ngxl+ngxpv     
ngyg=6#6
ngyeva=21#21
ngypv=9#9
ngyted=5#5 
ngytot=ngyg+ngyeva+ngyted

dx=largo_m/(ngx-1)      #[m] delta x
dx2=dx*1000         #[mm]
dyg=(tglass)/(ngyg-1)   #[m] delta y del vidrio       
dye=(teva)/(ngyeva-1)
dypv=(tpv)/(ngypv-1)
dyted=(tted)/(ngyted-1)

nm=ngx*ngytot
dt=1

L=1         #[m] largo del panel

#PROPIEDADES OPTICAS DEL VIDRIO
e_glass=0.91     #emisividad dle vidrio
abs_glass=0.9   #absortividad del vidrio
trans_glass=0.9 #transmitividad del vidrio
e_ted=0.85     #emisividad de tedlar
e_eva=0.9

#PROPIEDADES AIRE A TEMP.
T_inf=25+273                #[K] temperatura ambente en celsius
T_sky=T_inf
S=850                  #[W/m2] Radiacion constante 
efi_pv=0.1491           #eficiencia del panel % 
A_pv=dx*dx            #area de celda [m2]
V_pv=dx*dypv*dx         #Volumen de celda [m3]
A_ted=dx*dx          #AreaDetedlar [m2]
V_ted=dx*dyted*dx       #Volumen de tedlar [m3]

V=2                      # [m/s] velocidad del aire
v_air=1.59*10**(-5)         #[m2/s]
Pr=0.707                    #prandtl aire en temperatura ambiente
k_air=0.0263                #[W/mK]
Re=V*L/v_air                #Reynolds
Nu_l=0.664*Re**(0.5)*Pr**(1/3)
h_conv=Nu_l*k_air/L          # W/m2K Coeficiente convectivo natural del aire
sigma=5.670*10**(-8)        #[W/m2*K4] constante Stefan-Boltzmann

rad_on=-1

#CONSTANTES EN ECUACUIONES:
ated=(dx/dyted+dyted/dx)
bted=((dx*dyted)/dt)
aglass=(dx/dyg+dyg/dx)
bglass=((dx*dyg)/dt)
aeva=(dx/dye+dye/dx)
beva=((dx*dye)/dt)
apv=(dx/dypv+dypv/dx)
bpv=((dx*dypv)/dt)

#GENERACION CALOR EN CELDA
Q_pv=S*A_pv*(1-efi_pv)/V_pv

#GENERACION CALOR EN TEDLAR
Q_ted=S*A_ted*1/V_ted


Ti=numpy.ones(nm)*(25+273)            #Temperatura inicial placa acero K
B_copy=numpy.zeros((ngytot,ngx))

#CANAL DE AIRE ACTUALIZACON DE TEMPERATURA DE CELDA 

# V_i=3                      # [m/s] velocidad del aire
# v_iair=1.59*10**(-5)         #[m2/s]
# Pr_i=0.707                    #prandtl aire en temperatura ambiente
# k_iair=0.0263                #[W/mK]
# Re_i=V*L/v_air                #Reynolds
# Nu_l_i=0.664*Re**(0.5)*Pr**(1/3)
# h_forced=Nu_l*k_air/L          # W/m2K Coeficiente convectivo forzado del aire
# m_air=0.0023                           #kg/s
# cp_air=1007                        #kJ/kg K
# alpha=h_forced*dx*dx/(m_air*cp_air)
# dy_aire=dx

# T_ia=numpy.ones(ngx)*(20+273)            #Temperatura inicial aire K
# #T_ia[0]=25+273                  #Temperatura inicial aire K
# vect_ngx=numpy.linspace(0, ngx, ngx)

# #CANAL DE AIRE ACTUALIZACON DE TEMPERATURA DE CELDA 

V_i=2                      # [m/s] velocidad del aire
v_iair=1.59*10**(-5)         #[m2/s]
Pr_i=0.707                    #prandtl aire en temperatura ambiente
k_iair=0.0263                #[W/mK]
rho_aire=1.184                  #kg/m3 densidad del aire a 25° C 
dy_aire=2/100
ancho_aire=0.3                 #m ancho canal de aire
A_aire=ancho_aire*dy_aire       #Area transversal del ducto del aire 

D_h=(4*dy_aire*ancho_aire)/(2*dy_aire+2*ancho_aire)
AR=ancho_aire/dy_aire
Re_i=V_i*D_h/v_air                #Reynolds
Nu_l_i=7.54+(0.03*(D_h/L)*Re_i*Pr)**(2/3)/(1+0.0016*((D_h/L)*Re_i*Pr)**(2/3))#7.541*(1-2.610*AR+4.970*AR**2-5.119*AR**3+2.702*AR**4-0.548*AR**5)#7.54+(0.03*(D_h/L)*Re_i*Pr)/(1+0.0016*((D_h/L)*Re_i*Pr)**(2/3))
h_forced=Nu_l*k_air/D_h          # W/m2K Coeficiente convectivo forzado del aire
m_air=0.0147                           #0.018 kg/s
cp_air=1007                        #J/kg K
alpha=h_forced*dx*dx/(m_air*cp_air)

T_ia=numpy.ones(ngx)*(25+273)            #Temperatura inicial aire K
#T_ia[0]=25+273                  #Temperatura inicial aire K
vect_ngx=numpy.linspace(0, ngx, ngx)


##GLAZING 
alpha_glass=k_glass/p_glass
Ti_g=numpy.ones(ngx)*(25+273)
delta=dt*alpha_glass/dx*dx
th_g=6/1000 #espesor cubierta de vidiro metros


# AIRE ESTANCADO EN GAP

ka_gap=0.02808    #W/mK
v_agap=1.896*10**(-5)
Pr_gap=0.7202    
th_gap=2/100    #m
beta=1/T_inf        #1/K
g=9.81          #m/s2
Ts_glass=60+273     #K
Lc=1            #m2
Ra_l=((g*beta*(Ts_glass-T_inf)*Lc**3)*Pr_gap)/(v_agap**2)
ter1=1-(1708/Ra_l)
ter2=(Ra_l**(1/3)/18)-1
print(ter1)
print(ter2)
Nu_gap=1+1.44*ter1+ter2
h_gap=(ka_gap/Lc)*Nu_gap

R_gap=th_gap/(ka_gap*dx)
# R_gap=1/(h_gap*dx)

#RESISTENCIAS MATERIAL 

R_conv_f=1/(h_forced*dx)
R_ted=tted/(k_ted*dx)
R_eva=teva/(k_eva*dx)
R_pv=tpv/(k_pv*dx)
R_glass=tglass/(k_glass*dx)
R_conv=1/(h_conv*dx)


R_tot1=R_glass+R_gap
R_tot=R_conv_f+R_ted+R_eva+R_pv+R_eva+R_glass+R_conv
alpha2=dx/(R_tot*m_air*cp_air)







def matriz (nm,ngytot,ngx,dx,h_conv,k_glass,k_pv,k_eva,k_ted,p,dt,dx2,R_conv,R_tot1,Ti_g):
    A=numpy.zeros((nm,nm))
    im=0
    jm=0
    for jg in range (0,ngytot):
        for ig in range (0,ngx):
            
            
            if ig==0 and jg==0: #Esquina inferior izquierda (TEDLAR)
                h_rad=rad_on*(e_ted*sigma*(Ti[ig]**2+T_sky**2)*(Ti[ig]+T_sky))
                A[jm,jm]=-((k_ted/2)*ated +(h_rad/2+h_conv/2)*(dyted+dx)+(p_ted/4)*bted)
                A[jm,jm+1]=(k_ted/2)*dyted/dx
                A[jm,jm+ngx]=(k_ted/2)*dx/dyted
                # print((A[jm,jm]+A[jm,jm+1]+A[jm,jm+ngx])*(25+273))

            elif ig==ngx-1 and jg==0: #Esquina inferior derecha (TEDLAR)
                h_rad=rad_on*(e_ted*sigma*(Ti[ig]**2+T_sky**2)*(Ti[ig]+T_sky))
                A[jm,jm]=-((h_rad+h_conv)*(dx)+k_ted*ated+(p_ted/2)*bted)
                A[jm,jm-1]=k_ted*dyted/dx
                A[jm,jm+ngx]=k_ted*dx/dyted
                # print((A[jm,jm]+A[jm,jm-1]+A[jm,jm+ngx])*(25+273))
                
            elif ig==ngx-1 and jg==ngytot-1: #Esquina superior derecha (VIDRIO)
                h_rad=rad_on*(e_glass*sigma*(Ti[-1]**2+T_sky**2)*(Ti[-1]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                A[jm,jm]=-(1/((1/(h_rad*dx)+1/R_conv)**(-1)+R_tot1)+k_glass*aglass+(p_glass/2)*bglass)
                A[jm,jm-1]= k_glass*dyg/dx
                A[jm,jm-ngx]=k_glass*dx/dyg
                # print((A[jm,jm]+A[jm,jm-1]+A[jm,jm-ngx])*(25+273))
                

            elif ig==0 and jg==ngytot-1:#Esquina superior Izquierda (VIDRIO)
                h_rad=rad_on*(e_glass*sigma*(Ti[-ngx]**2+T_sky**2)*(Ti[-ngx]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                A[jm,jm]=-((h_rad/2+h_conv/2)*(dyg)+1/((2/(h_rad*dx)+2/R_conv)**(-1)+R_tot1)+(k_glass/2)*aglass+(p_glass/4)*bglass)
                A[jm,jm+1]=(k_glass/2)*dyg/dx
                A[jm,jm-ngx]=(k_glass/2)*dx/dyg
                # print((A[jm,jm]+A[jm,jm+1]+A[jm,jm-ngx])*(25+273))
                
                
            elif ig==0 and jg==ngyted-1:#Esquina izquierda union tedlar y eva
                h_rad_ted=rad_on*(e_ted*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_ted/2)*dyted+(h_rad/2)*dye+(h_conv/2)*(dye+dyted)+(k_ted/2)*ated+(k_eva/2)*aeva+((p_eva)/4)*beva+((p_ted)/4)*bted)
                A[jm,jm-ngx]=(k_ted/2)*(dx/dyted)
                A[jm,jm+1]=(k_eva/2)*(dye/dx) +(k_ted/2)*(dyted/dx)
                A[jm,jm+ngx]=(k_eva/2)*(dx/dye)
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+1]+A[jm,jm+ngx])*(25+273))
                
            elif ig==ngx-1 and jg==ngyted-1: #Esquina derecha union tedlar y eva            
                A[jm,jm]=-((k_ted)*(dx/dyted+(dyted/(dx)))+(k_eva)*(dx/dye+(dye/(dx)))+(p_eva/2)*beva+(p_ted/2)*bted)
                A[jm,jm-ngx]=k_ted*dx/dyted
                A[jm,jm-1]=((k_eva)*dye/dx +k_ted*dyted/dx)
                A[jm,jm+ngx]=k_eva*dx/dye
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm-1]+A[jm,jm+ngx])*(25+273))
            
            elif ig==0 and jg==(ngytot-ngyg)-1:#Esquina Izquierda Eva y VIdrio (valor alto 2021)
                h_rad_g=rad_on*(e_glass*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_g/2)*dyg+(h_rad_e/2)*dye+(h_conv/2)*(dyg+dye)+(k_eva/2)*aeva+(k_glass/2)*aglass+(p_eva/4)*beva+(p_glass/4)*bglass)
                A[jm,jm+ngx]=(k_glass/2)*dx/dyg 
                A[jm,jm-ngx]=(k_eva/2)*dx/dye
                A[jm,jm+1]=(k_glass/2)*dyg/dx+(k_eva/2)*dye/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm+1])*(25+273))
                
            elif ig==ngx-1 and jg==(ngytot-ngyg)-1:#Esquina derecha Eva y VIdrio (valor alto 4042)
                A[jm,jm]=-(k_glass*aglass+k_eva*aeva+(p_eva/2)*beva+(p_glass/2)*bglass)
                A[jm,jm+ngx]=k_glass*dx/dyg 
                A[jm,jm-ngx]=k_eva*dx/dye
                A[jm,jm-1]=k_glass*dyg/dx+k_eva*dye/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm-1])*(25+273))
                
            elif ig==ngx-1 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Esquina inferior derecha EVA Y Celda
                A[jm,jm]=-(((k_eva)*aeva+(k_pv)*apv)+((p_eva)/2)*beva+(p_pv/2)*bpv)
                A[jm,jm-ngx]=k_eva*dx/dye
                A[jm,jm+ngx]=k_pv*dx/dypv
                A[jm,jm-1]=(k_pv)*dypv/dx +(k_eva)*dye/dx
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm-1])*(25+273))
                
            elif ig==ngx-1 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#esquina derecha superior celda y eva
                A[jm,jm]=-((k_pv)*apv+(k_eva)*aeva+((p_eva/2)*beva+(p_pv/2)*bpv))
                A[jm,jm+ngx]=k_eva*dx/dye
                A[jm,jm-ngx]=k_pv*dx/dypv
                A[jm,jm-1]=(k_pv)*(dypv/dx)+(k_eva)*(dye/dx)
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm-1])*(25+273))
            
            elif ig==0 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Esquina interna izquierda inferior celda y eva
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_pv+h_rad_e)*dypv/2+h_conv*dypv+((k_eva/2)*aeva+(k_pv/2)*apv)+((p_eva)/4)*beva+(p_pv/4)*bpv)
                A[jm,jm+ngx]=(k_pv/2)*dx/dypv
                A[jm,jm-ngx]=(k_eva/2)*dx/dye
                A[jm,jm+1]=(k_pv/2)*dypv/dx+(k_eva/2)*dye/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+ A[jm,jm+1])*(25+273))
                
            elif ig==0 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Esquina superior izquierda interna celda y eva
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_pv+h_rad_e)*dypv/2+h_conv*dypv+((k_eva/2)*aeva+(k_pv/2)*apv)+((p_eva)/4)*beva+(p_pv/4)*bpv)
                A[jm,jm+ngx]=(k_eva/2)*dx/dye
                A[jm,jm-ngx]=(k_pv/2)*dx/dypv
                A[jm,jm+1]=(k_eva/2)*dye/dx+(k_pv/2)*dypv/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+ A[jm,jm+1])*(25+273))
                
                

            elif jg==ngytot-1: #Borde Superior (VIDRIO)
                h_rad=rad_on*(e_glass*sigma*(Ti[-ngx+ig]**2+T_sky**2)*(Ti[-ngx+ig]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                A[jm,jm]=-(1/((1/(h_rad*dx)+1/R_conv)**(-1)+R_tot1)+k_glass*aglass+(p_glass/2)*bglass)
                A[jm,jm-1]=(k_glass/2)*(dyg/dx)
                A[jm,jm+1]=(k_glass/2)*(dyg/dx)
                A[jm,jm-ngx]=k_glass*(dx/dyg)
                # print((A[jm,jm]+A[jm,jm-1]+A[jm,jm+1]+A[jm,jm-ngx])*(25+273))
                
            elif jg==0: #Borde Inferior (TEDLAR)
                h_rad=rad_on*(e_ted*sigma*(Ti[ig]**2+T_sky**2)*(Ti[ig]+T_sky))
                A[jm,jm]=-((h_rad+h_conv)*dx+k_ted*ated+(p_ted/2)*bted)
                A[jm,jm-1]=(k_ted/2)*dyted/dx
                A[jm,jm+1]=(k_ted/2)*dyted/dx
                A[jm,jm+ngx]=k_ted*dx/dyted
                # print((A[jm,jm]+A[jm,jm-1]+A[jm,jm+1]+A[jm,jm+ngx])*(25+273))
            
            elif ig==0 and 0<=jg<=ngyted-1:# Borde izquierdo Tedlar
                h_rad_ted=rad_on*(e_ted*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_conv+h_rad_ted)*dyted+k_ted*ated+(p_ted/2)*bted)
                A[jm,jm+ngx]=(k_ted/2)*dx/dyted
                A[jm,jm-ngx]=(k_ted/2)*dx/dyted
                A[jm,jm+1]=k_ted*dyted/dx
                # print((A[jm,jm]+ A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm+1])*(25+273))
                
            elif ig==ngx-1 and 0<=jg<=ngyted-1:# Borde derecho Tedlar
                A[jm,jm]=-(2*k_ted*ated+p_ted*bted)
                A[jm,jm+ngx]=k_ted*dx/dyted
                A[jm,jm-ngx]=k_ted*dx/dyted
                A[jm,jm-1]=2*k_ted*dyted/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm-1])*(25+273))
            
            elif jg==ngyted-1: #Borde Tedlar y Eva
                A[jm,jm]=-((k_ted)*ated+(k_eva)*aeva+((p_ted/2)*bted+(p_eva/2)*beva))
                A[jm,jm-ngx]=k_ted*dx/dyted
                A[jm,jm+1]=(k_ted/2)*dyted/dx+(k_eva/2)*dye/dx
                A[jm,jm-1]=(k_ted/2)*dyted/dx+(k_eva/2)*dye/dx
                A[jm,jm+ngx]=k_eva*dx/dye
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+1]+A[jm,jm-1]+A[jm,jm+ngx])*(25+273))
                
            elif jg==(ngytot-ngyg)-1:#Borde EVA Y VIDRIO
                A[jm,jm]=-(k_glass*aglass+k_eva*aeva+(p_eva/2)*beva+(p_glass/2)*bglass)
                A[jm,jm-ngx]=k_eva*dx/dye
                A[jm,jm+ngx]=k_glass*dx/dyg 
                A[jm,jm+1]=((k_eva/2)*dye/dx+(k_glass/2)*dyg/dx)
                A[jm,jm-1]=((k_eva/2)*dye/dx+(k_glass/2)*dyg/dx)
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+1]+A[jm,jm-1]+A[jm,jm+ngx])*(25+273))
            
            elif ig==0 and (((ngytot-1)-(ngyg-1))-1)<jg<ngytot-1:#Borde izquierdo de vidrio
                h_rad_g=rad_on*(e_glass*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_g+h_conv)*dyg+k_glass*aglass+(p_glass/2)*bglass)
                A[jm,jm+ngx]=(k_glass/2)*dx/dyg 
                A[jm,jm-ngx]=(k_glass/2)*dx/dyg 
                A[jm,jm+1]=k_glass*dyg/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm+1])*(25+273))
            
            elif ig==ngx-1 and (((ngytot-1)-(ngyg-1))-1)<jg<ngytot-1:#Borde derecho de vidrio
                A[jm,jm]=-(2*k_glass*aglass+(p_glass)*bglass)
                A[jm,jm+ngx]=(k_glass*(dx/dyg)) 
                A[jm,jm-ngx]=(k_glass*(dx/dyg)) 
                A[jm,jm-1]=2*k_glass*(dyg/dx)
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm-1])*(25+273))
            
            
            elif ig==0 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Borde izquierdo celda
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_pv+h_conv)*dypv+(k_pv*apv)+(p_pv/2)*bpv)
                A[jm,jm+1]=k_pv*(dypv/dx)
                A[jm,jm+ngx]=(k_pv/2)*dx/dypv
                A[jm,jm-ngx]=(k_pv/2)*dx/dypv
                # print((A[jm,jm]+A[jm,jm+1]+A[jm,jm+ngx]+A[jm,jm-ngx])*(25+273))
            
            elif ig==ngx-1 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Borde derecho celda
                A[jm,jm]=-(2*k_pv*apv+p_pv*bpv)
                A[jm,jm-1]=(2*k_pv*dypv/dx)
                A[jm,jm+ngx]=k_pv*dx/dypv
                A[jm,jm-ngx]=k_pv*dx/dypv
                # print(( A[jm,jm]+A[jm,jm-1]+A[jm,jm+ngx]+A[jm,jm-ngx])*(25+273))
            
            elif ig==0 :# Borde izquierdo 
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                A[jm,jm]=-((h_rad_e+h_conv)*dye+k_eva*aeva+(p_eva/2)*beva)
                A[jm,jm+ngx]=(k_eva/2)*dx/dye
                A[jm,jm-ngx]=(k_eva/2)*dx/dye
                A[jm,jm+1]=k_eva*dye/dx
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm+1])*(25+273))
            
            elif ig==ngx-1 : #Borde derecho 
                A[jm,jm]=-(2*k_eva*aeva+p_eva*beva)
                A[jm,jm+ngx]=k_eva*dx/dye
                A[jm,jm-ngx]=k_eva*dx/dye
                A[jm,jm-1]=(2*k_eva)*(dye/dx)
                # print((A[jm,jm]+A[jm,jm+ngx]+A[jm,jm-ngx]+A[jm,jm-1])*(25+273))
                
            elif 0<ig<ngx-1 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Borde interior inferior celda y Eva
                A[jm,jm]=-((k_eva)*aeva+(k_pv)*apv+((p_eva/2)*beva+(p_pv/2)*bpv))
                A[jm,jm-ngx]=k_eva*dx/dye
                A[jm,jm+ngx]=k_pv*dx/dypv
                A[jm,jm+1]=(k_pv/2)*dypv/dx+(k_eva/2)*dye/dx
                A[jm,jm-1]=(k_pv/2)*dypv/dx+(k_eva/2)*dye/dx
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
            
            elif 0<ig<ngx-1 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1): #Borde interior superior celda y eva
                A[jm,jm]=-((k_eva)*aeva+(k_pv)*apv+((p_eva/2)*beva+(p_pv/2)*bpv))
                A[jm,jm+ngx]=k_eva*dx/dye
                A[jm,jm-ngx]=k_pv*dx/dypv
                A[jm,jm-1]=(k_pv/2)*dypv/dx+(k_eva/2)*dye/dx
                A[jm,jm+1]=(k_pv/2)*dypv/dx+(k_eva/2)*dye/dx
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
            
            
                
            elif 0<ig<ngx-1 and 0<jg<ngyted-1:# Nodos internos tedlar 
                A[jm,jm]=-(2*k_ted*ated+p_ted*bted)
                A[jm,jm+1]=k_ted*dyted/dx
                A[jm,jm-1]=k_ted*dyted/dx
                A[jm,jm+ngx]=k_ted*dx/dyted
                A[jm,jm-ngx]=k_ted*dx/dyted
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
            
            elif 0<ig<ngx-1 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1): #Nodos internos Celda
                A[jm,jm]=-(2*k_pv*apv+p_pv*bpv)
                A[jm,jm-1]=k_pv*dypv/dx
                A[jm,jm+1]=k_pv*dypv/dx
                A[jm,jm+ngx]=k_pv*dx/dypv
                A[jm,jm-ngx]=k_pv*dx/dypv
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
            
            elif 0<ig<ngx-1 and (ngytot-ngyg)-1<jg<ngytot-1: #Nodos internos vidrio
                A[jm,jm]=-(2*k_glass*aglass+p_glass*bglass)
                A[jm,jm-1]=k_glass*dyg/dx
                A[jm,jm+1]=k_glass*dyg/dx
                A[jm,jm+ngx]=k_glass*dx/dyg
                A[jm,jm-ngx]=k_glass*dx/dyg
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
                
            elif 0<ig<ngx-1 and 0<jg<ngytot-1:# Nodos internos 
                A[jm,jm]=-(2*k_eva*aeva+p_eva*beva)
                A[jm,jm+1]=k_eva*dye/dx
                A[jm,jm-1]=k_eva*dye/dx
                A[jm,jm+ngx]=k_eva*dx/dye
                A[jm,jm-ngx]=k_eva*dx/dye
                # print((A[jm,jm]+A[jm,jm-ngx]+A[jm,jm+ngx]+A[jm,jm+1]+A[jm,jm-1])*(25+273))
            

                
            jm+=1
    return A

def RHS(nm,ngytot,ngx,dx,h_conv,dx2,T_inf,Ti,T_air,R_conv,R_tot1,Ti_g):
       
    C=numpy.zeros(nm)
    
    jm=0
    for jg in range (0,ngytot):
        for ig in range (0,ngx):
            if ig==0 and jg==0: #Esquina inferior izquierda (TEDLAR)
                h_rad=rad_on*(e_ted*sigma*(Ti[0]**2+T_sky**2)*(Ti[0]+T_sky))
                C[jm]=-(((h_conv/2)*(dyted+dx)*T_inf+(h_rad/2)*(dyted+dx)*T_sky) +(p_ted/4)*bted*Ti[jm])
                # print(C[jm])
                                
            elif ig==ngx-1 and jg==0: #Esquina inferior derecha (TEDLAR)
                h_rad=rad_on*(e_ted*sigma*(Ti[ngx-1]**2+T_sky**2)*(Ti[ngx-1]+T_sky))
                C[jm]=-(((h_conv)*dx*T_inf+(h_rad)*dx*T_sky)+(p_ted/2)*bted*Ti[jm])
                # print(C[jm])
            
            elif ig==ngx-1 and jg==ngytot-1: #Esquina superior derecha (vidrio)
                h_rad=rad_on*(e_glass*sigma*(Ti[-1]**2+T_sky**2)*(Ti[-1]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                C[jm]=-((p_glass/2)*bglass*Ti[jm]+(1/((1/(h_rad*dx)+1/R_conv)**(-1)+R_tot1))*T_inf)
                # print(C[jm])
            
            elif ig==0 and jg==ngytot-1:#Esquina superior Izquierda (VIDRIO)
                h_rad=rad_on*(e_glass*sigma*(Ti[-ngx]**2+T_sky**2)*(Ti[-ngx]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                C[jm]=-((p_glass/4)*bglass*Ti[jm]+(h_conv/2)*(dyg)*T_inf+(h_rad/2)*(dyg)*T_sky+(1/((2/(h_rad*dx)+2/R_conv)**(-1)+R_tot1))*T_inf) 
                # print(C[jm])
            
            elif ig==0 and jg==ngyted-1:#Esquina izquierda union tedlar y eva
                h_rad_ted=rad_on*(e_ted*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-(((p_eva/4)*beva+(p_ted/4)*bted)*Ti[jm]+(h_conv/2)*(dye+dyted)*T_inf+((h_rad_ted/2)*dyted+(h_rad_e/2)*dye)*T_sky)
                # print(((p_eva/4)*beva+(p_ted/4)*bted)*Ti[jm]+(h_conv/2)*(dye+dyted)*T_inf+((dx/2)*(dyted/2)*Q_ted))
                # print(C[jm])
            
            elif ig==ngx-1 and jg==ngyted-1: #Esquina derecha union tedlar y eva 
                C[jm]=-(((p_eva/2)*beva+(p_ted/2)*bted)*Ti[jm])
                # print(C[jm])
            
            elif ig==0 and jg==(ngytot-ngyg)-1:#Esquina Izquierda Eva y VIdrio (valor grande 2021 comparado)
                h_rad_g=rad_on*(e_glass*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-(((p_eva/4)*beva+(p_glass/4)*bglass)*Ti[jm]+(h_conv/2)*T_inf*(dyg+dye)+((h_rad_g/2)*dyg+(h_rad_e/2)*dye)*T_sky)
                # print(C[jm])
            elif ig==ngx-1 and jg==(ngytot-ngyg)-1:#Esquina derecha Eva y VIdrio 
                C[jm]=-(((p_eva/2)*beva+(p_glass/2)*bglass)*Ti[jm])
                # print(C[jm])
            
            elif ig==ngx-1 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Esquina inferior derecha EVA Y Celda
                C[jm]=-(((p_eva/2)*beva+(p_pv/2)*bpv)*Ti[jm])
                # print(C[jm])
            
            elif ig==ngx-1 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#esquina derecha superior celda y eva
                C[jm]=-(((p_eva/2)*beva+(p_pv/2)*bpv)*Ti[jm]+(dx*(dypv/2)*Q_pv))
                # print(C[jm])
            
            elif ig==0 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Esquina interna izquierda inferior celda y eva
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-(((p_eva/4)*beva+(p_pv/4)*bpv)*Ti[jm]+h_conv*dypv*T_inf+h_rad_pv*dypv*T_sky)
                # print(C[jm])
            
            elif ig==0 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Esquina superior izquierda interna celda y eva
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-(((p_eva/4)*beva+(p_pv/4)*bpv)*Ti[jm]+h_conv*dypv*T_inf+h_rad_pv*dypv*T_sky+((dx/2)*(dypv/2)*Q_pv))
                # print(C[jm])
            
            elif jg==ngytot-1: #Borde Superior (VIDRIO)
                h_rad=rad_on*(e_glass*sigma*(Ti[-ngx+ig]**2+T_sky**2)*(Ti[-ngx+ig]+T_sky))
                h_rgap=h_rad/(((1-e_glass)/e_glass) +1)
                C[jm]=-((p_glass/2)*bglass*Ti[jm]+(1/((1/(h_rad*dx)+1/R_conv)**(-1)+R_tot1))*T_inf)
                # print(C[jm])
            
            elif jg==0: #Borde Inferior
                h_rad=rad_on*(e_ted*sigma*(Ti[ig]**2+T_sky**2)*(Ti[ig]+T_sky))
                C[jm]=-((p_ted/2)*bted*Ti[jm]+h_conv*dx*T_air[jm]+h_rad*dx*T_sky)
                # print(C[jm])
                
            elif ig==0 and 0<=jg<=ngyted-1:# Borde izquierdo Tedlar
                h_rad_ted=rad_on*(e_ted*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-((p_ted/2)*bted*Ti[jm]+h_conv*dyted*T_inf+h_rad_ted*dyted*T_sky)
                # print(C[jm])
            
            elif ig==ngx-1 and 0<=jg<=ngyted-1:# Borde derecho Tedlar
                C[jm]=-((p_ted)*bted*Ti[jm])
                # print(C[jm])
            
            # elif 0<ig<ngxl-1 and jg==ngyted-1: #Borde Tedlar y Eva seccion que genera calor
            #     C[jm]=-(((p_ted/2)*bted +(p_eva/2)*beva)*Ti[jm]+((dx)*(dyted/2)*Q_ted)) #
            #     # print(((p_ted/2)*bted +(p_eva/2)*beva)*Ti[jm]+((dx)*(dyted/2)*Q_ted))
                
            
            elif jg==ngyted-1: #Borde Tedlar y Eva
                C[jm]=-(((p_ted/2)*bted +(p_eva/2)*beva)*Ti[jm])
                # print(((p_ted/2)*bted +(p_eva/2)*beva)*Ti[jm])
                # print(C[jm])
                
            
            
            elif jg==(ngytot-ngyg)-1:#Borde EVA Y VIDRIO
                C[jm]=-(((p_eva/2)*beva+(p_glass/2)*bglass)*Ti[jm])
                # print(C[jm])
            
            elif ig==0 and (((ngytot-1)-(ngyg-1))-1)<jg<ngytot-1:#Borde izquierdo de vidrio
                h_rad_g=rad_on*(e_glass*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-((p_glass/2)*bglass*Ti[jm]+h_conv*dyg*T_inf+h_rad_g*dyg*T_sky)
                # print(C[jm])
            
            elif ig==ngx-1 and (((ngytot-1)-(ngyg-1))-1)<jg<ngytot-1:#Borde derecho de vidrio
                C[jm]=-((p_glass)*bglass*Ti[jm])
                # print(C[jm])
            
            elif ig==0 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Borde izquierdo celda
                h_rad_pv=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-(((p_pv/2)*bpv)*Ti[jm]+h_conv*dypv*T_inf+h_rad_pv*dypv*T_sky)
                # print(C[jm])
            
            elif ig==ngx-1 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1):#Borde derecho celda
                C[jm]=-(p_pv*bpv*Ti[jm])
                # print(C[jm])
                
            elif ig==0 :# Borde izquierdo
                h_rad_e=rad_on*(e_eva*sigma*(Ti[jg*ngx]**2+T_sky**2)*(Ti[jg*ngx]+T_sky))
                C[jm]=-((p_eva/2)*beva*Ti[jm]+h_conv*dye*T_inf+h_rad*dye*T_sky)
                # print(C[jm])
                
            
            elif ig==ngx-1 : #Borde derecho 
                C[jm]=-((p_eva)*beva*Ti[jm])
                # print(C[jm])
            
                
            elif 0<ig<ngx-1 and jg==((ngyted+(ngyeva-ngypv)/2)-1): #Borde interior inferior celda y Eva
                C[jm]=-(((p_eva/2)*beva+(p_pv/2)*bpv)*Ti[jm])
                # print(C[jm])
            
            elif 0<ig<ngx-1 and jg==((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1): #Borde interior superior celda y eva
                C[jm]=-(((p_eva/2)*beva+(p_pv/2)*bpv)*Ti[jm]+((dx*dypv/2)*Q_pv))
                # print(C[jm])
             
            elif 0<ig<ngx-1 and 0<jg<ngyted-1:# Nodos internos en tedlar
                C[jm]=-((p_ted)*bted*Ti[jm])
                # print(C[jm])
            
            elif 0<ig<ngx-1 and ((ngyted+(ngyeva-ngypv)/2)-1)<jg<((ngyted+((ngyeva-ngypv)/2)+ngypv-1)-1): #Nodos internos Celda
                C[jm]=-((p_pv)*bpv*Ti[jm])
                # print(C[jm])
                
            elif 0<ig<ngx-1 and (ngytot-ngyg)-1<jg<ngytot-1: #Nodos internos vidrio
                C[jm]=-((p_glass)*bglass*Ti[jm])
                # print(C[jm])
            
            elif 0<ig<ngx-1 and 0<jg<ngytot-1:# Nodos internos 
                C[jm]=-((p_eva)*beva*Ti[jm])
                # print(C[jm])

            jm+=1
    return C

def implicito(C, A):
    T = linalg.solve(A,C)
    
        
    return T

A=matriz(nm,ngytot,ngx,dx,h_conv,k_glass,k_pv,k_eva,k_ted,p,dt,dx2,R_conv,R_tot1,Ti_g)



###### GLAZING 

def GenMatriz(ngx, sigma, e_glass, Ti_g,T_sky,h_conv,k_glass,th_g,dx,p_glass):
    D=numpy.zeros((ngx,ngx))
    jm=0
    for ig in range (0,ngx):
        if ig==0:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[0]**2+T_sky**2)*(Ti_g[0]+T_sky))
            D[jm,jm]=-((h_conv)*(th_g+dx)+(h_rad)*(th_g+dx)+(p_glass)*th_g*dx/dt+k_glass*th_g/dx)
            D[jm,jm+1]=k_glass*th_g/dx
            # print((D[jm,jm]+D[jm,jm+1])*(25+273))
            
        elif ig==ngx-1:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[ig]**2+T_sky**2)*(Ti_g[ig]+T_sky))
            D[jm,jm]=-((h_conv)*(th_g+dx)+(h_rad)*(th_g+dx)+(p_glass)*th_g*dx/dt+k_glass*th_g/dx)
            D[jm,jm-1]=k_glass*th_g/dx
            # print((D[jm,jm]+D[jm,jm-1])*(25+273))
            
        else:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[ig]**2+T_sky**2)*(Ti_g[ig]+T_sky))
            D[jm,jm]=-((h_conv)*dx+(h_rad)*dx+(p_glass)*th_g*dx/dt+2*k_glass*th_g/dx)
            D[jm,jm-1]=k_glass*th_g/dx
            D[jm,jm+1]=k_glass*th_g/dx
            # print((D[jm,jm]+D[jm,jm+1]+D[jm,jm-1])*(25+273))
            
        jm+=1
            
    return D

# D=GenMatriz(ngx, sigma, e_glass, Ti_g,T_sky,h_conv,k_glass,th_g,dx,p_glass)


def lado_derecho(sigma,h_conv,T_sky,th_g,dx,Ti_g,T_inf,dt):
    rhs=numpy.zeros(ngx)
    jm=0
    
    for ig in range (0,ngx):
        if ig==0:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[0]**2+T_sky**2)*(Ti_g[0]+T_sky))
            rhs[jm]=-((h_conv)*(th_g+dx)*T_inf+(h_rad)*(th_g+dx)*T_sky+(p_glass)*th_g*dx/dt*Ti_g[ig])
            # print(rhs[jm])
                
            
        elif ig==ngx-1:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[ig]**2+T_sky**2)*(Ti_g[ig]+T_sky))
            rhs[jm]=-((h_conv)*(th_g+dx)*T_inf+(h_rad)*(th_g+dx)*T_sky+(p_glass)*th_g*dx/dt*Ti_g[ig])
            # print(rhs[jm])
               
            
        else:
            h_rad=rad_on*(e_glass*sigma*(Ti_g[ig]**2+T_sky**2)*(Ti_g[ig]+T_sky))
            rhs[jm]=-((h_conv)*(dx)*T_inf+(h_rad)*(+dx)*T_sky+(p_glass)*th_g*dx/dt*Ti_g[ig])
            # print(rhs[jm])
                
        jm+=1
    return rhs
  
# b=lado_derecho(sigma,h_conv,T_sky,th_g,dx,Ti_g,T_inf,dt)
# T_interior=solve(D,b)
# Ti_copy=T_interior.copy()            
    
    

nt=3600

vect_t=numpy.linspace(0, nt, nt)
t_pv=[]
t_ted=[]
t_vidrio=[]
t_eva=[]


##### CANAL DE AIRE INFERIOR 

def Temp_aire(ngx,T_ia,dx, h_forced,m_air,cp_air,T_s_ted):
    nx=len(T_air)
    for i in range(0,nx-2):
        # T_air[i+1]=(1/(1+alpha2))*T_air[i]+(alpha2/(1+alpha2))*T_s_ted[i]
        # #T_air[ngx-1]=T_air[ngx-2]+(T_air[ngx-2]-T_air[ngx-3])
        # T_air[ngx-1]=(1/(1+alpha2))*T_air[ngx-2]+(alpha2/(1+alpha2))*T_s_ted[i] 
        
        
        T_air[i+1]=(1/(1+alpha))*T_air[i]+(alpha/(1+alpha))*T_s_ted[i]
        #T_air[ngx-1]=T_air[ngx-2]+(T_air[ngx-2]-T_air[ngx-3])
        T_air[ngx-1]=(1/(1+alpha))*T_air[ngx-2]+(alpha/(1+alpha))*T_s_ted[i]
        
    return T_air

T_air=T_ia.copy()




######  INICIALIZACION FUNCIONES

for t in range(nt):
    
    C=RHS(nm,ngytot,ngx,dx,h_conv,dx2,T_inf,Ti,T_air,R_conv,R_tot1,Ti_g)
    T=implicito(C, A)
    Ti=T.copy()
    t_ted.append(T[35]-273)
    t_pv.append(T[185]-273)
    t_vidrio.append(T[275]-273)
    t_eva.append(T[195]-273)
    T_s_ted=T[:ngx]
    T_air=Temp_aire(ngx, T_ia, dx, h_forced, m_air, cp_air,T_s_ted)
    
    
    #   GLAZING
    D=GenMatriz(ngx, sigma, e_glass, Ti_g,T_sky,h_conv,k_glass,th_g,dx,p_glass)
    b=lado_derecho(sigma,h_conv,T_sky,th_g,dx,Ti_g,T_inf,dt)
    T_interior=solve(D,b)
    
    Ti_g=T_interior.copy()
    
    
    
    



vect_T_s=T_s_ted-273                # vector con las temperaturad del tedlar superficie
print("nt",nt)
print(T_s_ted-273)
        

#ORDEN DE TEMPERATURAS SEGUN MALLADO

B=numpy.zeros((ngytot,ngx))         #matriz para reordenar temperaturas segun posicion
# print(len(T))
for i in range (ngytot):
    B[(ngytot-i)-1,:]=T[i*ngx:(i+1)*(ngx)]-273
    B_copy=B
    

T_ted2=T[:(ngyted)*ngx]-273     #vector temperaturas de tedlar para sacar promedio de temperaturas
z=(sum(T_ted2)/len(T_ted2))     # promedio de temperaturas de tedlar
zz=(sum(T)/len(T))
print(zz-273)                   #Temperatura promedio panel
print(z)                        #temperatura promedio tedlar

print(T_air-273)
        
T_air_C=T_air-273               #TEMPERATURA DEL CANAL DE AIRE EN CELSIUS

    
    
            
        
plt.plot(vect_t,t_ted,label='tedlar')
plt.plot(vect_t,t_vidrio,label='vidrio')
plt.plot(vect_t,t_eva,label='eva')
plt.plot(vect_t,t_pv,label='celda')
plt.xlabel("tiempo s")
plt.ylabel("temperatura °C")
plt.legend()
plt.show()
plt.plot(vect_ngx,T_air_C,label='T aire en canal')
plt.xlabel("nodo x")
plt.ylabel("temperatura °C")
plt.legend()
plt.show()
plt.plot(vect_ngx,vect_T_s,label='T superficie tedlar')
plt.xlabel("nodo x")
plt.ylabel("temperatura °C")
plt.legend()
plt.show()
plt.plot(vect_ngx,T_air_C,label='T aire en canal')
plt.plot(vect_ngx,vect_T_s,label='T superficie tedlar')
plt.xlabel("nodo x")
plt.ylabel("temperatura °C")
plt.legend()
plt.show()