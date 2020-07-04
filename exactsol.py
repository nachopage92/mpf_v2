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
    young   = kwargs.get('E',None) # ojo aqui
    poisson = kwargs.get('v',None) # eventualmente produce error

    if ( test == 'patchtest-0' ):
        u = x+y
        dudx = 1.0
        dudy = 1.0
        v = x+y
        dvdx = 1.0
        dvdy = 1.0
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy

    elif (test=='patchtest-1'): # esfuerzo uniforme
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
        u    = x*y/young
        dudx = y/young
        dudy = x/young
        v    = -0.5*(x*x+poisson*y*y)/young
        dvdx = -x/young
        dvdy = -poisson*y/young
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        # note: sigma y es identicamente cero, para evitar decimales indeseados
        sigmay = 0.0
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy
        
    elif (test=='cantilever'): # viga en voladizo
        Lx = 8.0   ; Ly = 1.0
        P = 1000.0 ; I = 1.0*Ly*Ly*Ly/12.0
        c = Ly/2.0 ; G = 2.0*young/(1.0+poisson)
        # cambio de variables
        x_ = Lx-x
        y_ = y
        # -------------------------
        # solucion analitica:
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
        # regreso a la variable inicial ----
        v    = - v
        dvdx = - dvdx
        dvdy = - dvdy
        # ---------------------------------
        
        epsxx,epsyy,epsxy,sigmax,sigmay,tauxy = strain_stress_calculation()
        # los esfuerzos normales en y son nulos
        sigmay = 0.0 # para evitar errores numericos:
        return u,v,epsxx,epsyy,epsxy,sigmax,sigmay,tauxy

    elif (test=='infinite-plate'):
        from numpy import sqrt,power,pi,arccos,arcsin,cos,sin
        # -------------- contantes ------------------
        a = 1.0                       # radio del orificio
        k = (3.-poisson)/(1.+poisson) # kolosov constant
        G = 0.5*young/(1.+poisson)    # modulo cortante
        # -------------- coordenadas polares --------
        r = sqrt( power(x,2) + power(x,2) )
        if ( x <= 0 ):
            theta = 0.5*pi
        if ( y <= 0 ):
            theta = 0.0
        else
            theta = arccos(x/r) # = arcsin(y,r)
        #-------------------------------------------------------#
        # se calculan los siguientes elementos infinitesimales: #
        #   du/dtheta ; du/dradius     dtheta/dx ; drradius/dx  #
        #   dv/dtheta ; dv/dradius     dtheta/dy ; drradius/dy  #  
        # se realiza la siguiente operacion                     # 
        #          __                            __             #
        #   u,x    |   u,t    u,r      0      0   |  t,x        #  
        #   v,x  = |   v,t    v,r      0      0   |  r,x        #    
        #   u,y    |    0      0      u,t    u,r  |  t,y        #    
        #   v,y    |_   0      0      v,t    v,r _|  r,y        #     
        #                                                       #   
        # o bien, si definimos la submatriz M, pvecx, pvecy,    #    
        # uvecx, uvecy como:                                    #    
        #                                                       #      
        #  M = | u,t  u,r | ; pvecx = | t,x | ; pvecy = | t,y | #
        #      | v,t  v,r |           | r,x |           | r,y | #
        #                                                       #
        #  uvecx = | u,x | ; uvecy = | u,y |                    #
        #          | v,x |           | v,y |                    #     
        #                                                       #       
        # entonces:                                             #           
        #                uvecx = matmul( M , pvecx )            #   
        #                uvecy = matmul( M , pvecy )            #   
        #                                                       #
        # los terminos dtdx,dtdy,drdx,drdy se determinar indi-  #
        # rectamente calculada la matrix inversa:               #       
        #                                                       #
        #  | x,t   y,t |  =  | t,x  t,y |                       #
        #  | x,r   y,r |     | r,x  r,y |                       #
        #-------------------------------------------------------#
        f1 = (10.*a)/(8.*G)
        f2 = (r/a)*(k+1.)
        f3 = 2.*a/r
        f4 = (1.+k)
        f5 = 2.*power(a,3)/power(r,3)

        u = f1*( f2*cos(theta) + f3*(f4*cos(theta)+cos(3.*theta)) - f5*cos(3.*theta) ) 
        v = f1*( f2*sin(theta) + f3*(f4*sin(theta)+sin(3.*theta)) - f5*sin(3.*theta) ) 

        dvdt = f1*( -f2*sin(theta) - f3*(f4*sin(theta)+3.*sin(3.*theta)) + 3.*f5*sin(3.*theta) ) 
        dudt = f1*( f2*cos(theta) + f3*(f4*cos(theta)+3.*cos(3.*theta)) - 3.*f5*cos(3.*theta) ) 

        f1 = (10.*a)/(8.*G)
        f2 = (1/a)*(k+1.)
        f3 = 2.*a/power(r,2)
        f4 = (1.+k)
        f5 = 6.*power(a,3)/power(r,4)

        dudr = f1*( f2*cos(theta) + f3*(f4*cos(theta)+cos(3.*theta)) - f5*cos(3.*theta) ) 
        dvdr = f1*( f2*sin(theta) + f3*(f4*sin(theta)+sin(3.*theta)) - f5*sin(3.*theta) ) 
    


    else:
        print('Error: test no implementado')
        return


