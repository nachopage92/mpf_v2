def SOLUCION_EXACTA(test,x,y,**kwargs):
    
    #''''''''''''''''''''''''''''''''''''''''
    # ... strain_stress_calculation ...
    #
    def strain_stress_calculation():
        # CAMINO CORTO
        factor = young/(1.0-poisson*poisson)
        sigmax = factor*(dudx+dvdy*poisson)
        sigmay = factor*(dvdy+dudx*poisson)
        tauxy  = factor*0.5*(1.0-poisson)*(dudy+dvdx)
        epsxx  = dudx
        epsyy  = dvdy
        epsxy  = dudy+dvdx
        ## CAMINO LARGO
        #import numpy as np
        #factor = young/(1.0-poisson*poisson)
        #mat = [ [   1.0   , poisson , 0.0 ] , \
        #                         [ poisson ,   1.0   , 0.0 ] ,  \
        #                         [   0.0   ,   0.0   , 0.5*(1.0-poisson)] ]
        #mat = np.array(mat)
        #mat = np.multiply(mat,factor)
        #voigt_strain = [ dudx , dvdy , dudy+dvdx ]
        #voigt_strain = np.array(voigt_strain)
        #voigt_stress = np.matmul(mat,voigt_strain)
        #sigmax,sigmay,tauxy = voigt_stress
        return epsxx,epsyy,epsxy,sigmax,sigmay,tauxy

    # ''''''''''''''''''''''''''''''''''''''''''

    if ( test == 'patchtest-0' ):
        young   = kwargs.get('E',None) # ojo aqui
        poisson = kwargs.get('v',None) # eventualmente produce error
        u = x+y
        dudx = 1.0
        dudy = 1.0
        v = x+y
        dvdx = 1.0
        dvdy = 1.0
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy

    elif (test=='patchtest-1'): # esfuerzo uniforme
        young   = kwargs.get('E',None) # ojo aqui
        poisson = kwargs.get('v',None) # eventualmente produce error
        aux = x*0.0+1.0 # 'artilugio' para graficar
        u    = x/young
        dudx = 1.0/young * aux
        dudy = 0.0 * aux
        v    = -y*poisson/young
        dvdx = 0.0 * aux
        dvdy = -poisson/young * aux
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy
        
    elif (test=='patchtest-2'): # esfuerzo lineal
        young   = kwargs.get('E',None) # ojo aqui
        poisson = kwargs.get('v',None) # eventualmente produce error
        u    = x*y/young
        dudx = y/young
        dudy = x/young
        v    = -0.5*(x*x+poisson*y*y)/young
        dvdx = -x/young
        dvdy = -poisson*y/young
        print(u,dudx,dudy)
        print(v,dvdx,dvdy)
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        print(v,dvdx,dvdy)
        print(sigmax,sigmay,tauxy)
        input()
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy
        
    elif (test=='cantilever'): # viga en voladizo
        P,young,poisson,I,c,G = constantes(test)
        Lx_i,Lx_f,Ly_i,Ly_f = constantes_geometricas(test)
        Lx = Lx_f-Lx_i ; Ly=Ly_f-Ly_i
        # una pequena correccion
        x_ = Lx-x
        y_ = y
        u = -P*x_*x_*y_/(2.0*young*I) - poisson*P*y_*y_*y_/(6.0*young*I) + P*y_*y_*y_/(6.0*I*G)\
                + ( P*Lx
                        *Lx/(2.0*young*I) - P*c*c/(2.0*G*I) ) * y_
        dudx = P*x_*y_/(young*I)
        dudy = -P*x_*x_/(2.0*young*I) - poisson*P*y_*y_/(2.0*young*I) + P*y_*y_/(2.0*I*G)\
                + P*Lx*Lx/(2.0*young*I) - P*c*c/(2.0*I*G)
        
        v = poisson*P*x_*y_*y_/(2.0*young*I) + P*x_*x_*x_/(6.0*young*I) - P*Lx*Lx*x_/(2.0*young*I)\
                + P*Lx*Lx*Lx/(3.0*young*I)
        dvdx = -poisson*P*y_*y_/(2.0*young*I) - P*x_*x_/(2.0*young*I) + P*Lx*Lx/(2.0*young*I)
        dvdy = poisson*P*x_*y_/(young*I)
        # una pequena correccion
        v    = - v
        dvdx = - dvdx
        dvdy = - dvdy
        
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy

    elif (test=='infinite-plate'):
        print("PENDIENTE") ; from sys import exit ; exit('stop')
        #epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        #return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy
    
    else:
        print('Error: test no implementado')
        return


"""
constantes_geometricas
    Entrega los limites en cada dimension (x,y) del dominio.
    Supone dominio rectangular. Para el caso 'meshgen = nurbs'
    los limites corresponderan al espacio parametrico (xi,eta)
    con valores predeterminados (Lx_i=0;Lx_f=1;Ly_i=0;Ly_f=1)
    Lx_i : limite izquierdo direccion x         ______        
    Lx_f : limite derecho direccion x          | Ly_f |       
    Ly_i : limite inferior direccion y     Lx_i|      |Lx_f   
    Ly_f : limite superior direccion y         |______|       
                                                 Ly_i       
"""
def constantes_geometricas(test):
    if (test=='patchtest-0') :
        Lx_i = 0.0 ; Lx_f = 2.0
        Ly_i = 0.0 ; Ly_f = 2.0
    if (test=='patchtest-1') :
        Lx_i = 0.0 ; Lx_f = 1.0
        Ly_i = 0.0 ; Ly_f = 1.0
    elif (test=='patchtest-2') :
        Lx_i = 0.0 ; Lx_f = 1.0
        Ly_i = 0.0 ; Ly_f = 1.0
    elif (test=='cantilever') :
        Lx_i = 0.0 ; Lx_f = 8.0
        Ly_i = -0.5 ; Ly_f = 0.5
    else:
        print('Test no ingresado')
        from sys import exit ; exit('Stopped')
    return Lx_i,Lx_f,Ly_i,Ly_f


"""
    MAS CONSTANTES
"""
def constantes(test,**kwargs):
    if ( test == 'patchtest-0' or\
         test == 'patchtest-1' or\
         test == 'patchtest-2' ):
        return self.young_modulus, self.poisson_coeff
    elif ( test == 'cantilever' ):
        Lx_i,Lx_f,Ly_i,Ly_f = constantes_geometricas(test)
        Lx=Lx_f-Lx_f ; Ly=Ly_f-Ly_i 
        young = self.young_modulus
        poisson = self.poisson_coeff
        P = 1.0
        #P = 1000.0
        I = 1.0*Ly*Ly*Ly/12.0
        c = Ly/2.0
        G = 2.0*E/(1.0+v)
        return P,E,v,I,c,G
    else:
        print('Test no ingresado!!')
        return
