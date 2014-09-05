c        PROGRAM  construido a partir de retardo.f (20sep06)
c
c		Jose M. Albornoz
c		Grupo de Caos y Sistemas COmplejos
c		Universidad de Los Andes	
c
C==========================================================
C       Dos enzimas (alpha y beta) y dos substancias (A y B)
c	que conforman un ciclo futil
c	
c parametros:
c     ga y gb son los coeficientes de las tasas de captura de A y B respectivamente
c     tau_ap y tau_bp son los tiempos de procesamiento de las encimas alpha y beta 
c     tau_a y tau_b son los tiempos de ciclo completo de las encimas alpha y beta
c     
c REGISTROS
c	apm=ga*a*alpha1 ======== tasa de destruccion de A
c	bpm=gb*b*beta1 ======== tasa de destruccion de B
c
c VALORES RETAZADOS
c	r_apm_half=apm(t-tau_ap)
c       r_apm_one=apm(t-tau_a)
c	r_bpm_half=bpm(t-tau_ap)
c       r_bpm_one=bpm(t-tau_a)
C==========================================================
c VARIABLES
c                   y1==A
c                   y2==alpha1
c                   y3==alpha2
c                   y4==alpha3
c                   y5==B
c                   y6==beta1
c                   y7==beta2
c                   y8==beta3
c  8 ec.  dif. ord.
c                 no son estas
c                   dy1/dt=r_bpm_half-ga*y1*y2   
c                   dy2/dt=r_apm_one-ga*y1*y2   
c                   dy3/dt=ga*y1*y2-r_apm_half   
c                   dy4/dt=r_apm_half-r_apm_one
c
c                   dy5/dt=r_apm_half-gb*y5*y6
c                   dy6/dt=r_bpm_one-gb*y5*y6
c                   dy7/dt=gb*y5*y6-r_bpm_half
c                   dy8/dt=r_bpm_half-r_bpm_one   
c
c  CONSTANTES Y VARIABLES:
c  np=Numero de puntos sobre el eje del tiempo.
c  kmax=Numero de pasos que pueden ser guardados.
c  hmin=Paso minimo permitido (puede ser cero)
c  eps=Nivel de tolerancia
c  SUBRUTINAS:
c  Init:Se dan los valores iniciales para arrancar integracion
c  Odeint:Utiliza el metodo Runge-Kutta con control adaptador del paso.
c  deri5:Calcula las derivadas del integrando (derint=dy/dm)
c  hunt:Dado un arreglo x(i) de longitud N,nos devuelve un valor jlo tal
c       que  x se encuentra entre x(jlo) y x(jlo+1).
c
        implicit double precision (a-h,o-z)
        external deri5,rkqc
        
        parameter(nvar5=8,np=30000)
        common/evovar/a0,alpha10,alpha20,alpha30,b0,beta10,beta20,beta30
        common/histprev/a00,alpha100,b00,beta100
        common/evoret/jlo
        common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
        dimension ys5(nvar5)
        common/path5/kmax5,kount5,dxsav5,x5(2000002),y5(10,2000002)
        common/pfijo/As,Bs,alpha1,alpha2,alpha3,beta1,beta2,beta3
        common/signal/isig,iflag
        common/sust/S,xmg_a
        common/otrospars/Vsim,xMMa,xMMb,xNAv
        common/vels/vmaxa,vmaxb,rka,rkb,ca,cb,gga,ggb
        open(unit=3,file='delay.dat')
c        open(unit=3,file='ret3_0.res')
        open(unit=4,file='ret4.res')
        open(unit=7,file='retardo.res')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        jlo=1
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c
c       numero de avogadro        
        xNAv=6.022d23

c       volumen de simulacion, litros       
        Vsim=10.3d-18
c        Vsim=5.15d-18

        if(.false.)then
           xMMa=5.d4
           xMMb=5.d4
           xNA=200.
           xNB=200.
c          concentracion final de sustrato, microMolar        
           S = 400.d0
        endif
        if(.true.)read(5,*)xMMa,xMMb,xNA,xNB,S,x25

	icaso=1

        if(icaso.eq.0)then
           vmaxa=150.0
           rka=0.005
           vmaxb=80.0
           rkb=0.005
           ca=0.5
           cb=0.5
c
           ga=vmaxa/rka
           gb=vmaxb/rkb
           tau_a=1.0/vmaxa
           tau_b=1.0/vmaxb
           tau_ap=ca*tau_a
           tau_bp=cb*tau_b
        endif

c       Las velocidades de las enzimas (va,vb) estan expresadas en micromol/(min-mg)
c       Las afinidades (rka,rkb) estan expresadas en microMolar
c
        if(icaso.eq.1)then
           if(.false.)then
              va = 500.
              vb = 500.

c             las velocidades son expresadas en micromol/(seg-mg)
              rka=10000.
              rkb=10000.
              ca=0.5
              cb=0.5
           endif

           if(.true.)read(5,*)va,vb,rka,rkb,ca,cb
        
           vmaxa=va/60.
           vmaxb=vb/60.
           ga=vmaxa/rka
           gb=vmaxb/rkb
           gga=ga
           ggb=gb

c          tiempos de procesamiento para UNA enzima           
           tau_a=60.0/(va*xMMa*1.d-3)
           tau_b=60.0/(vb*xMMb*1.d-3)
           tau_ap=ca*tau_a
           tau_bp=cb*tau_b
        endif

	if(icaso.eq.2)then
	   vmaxa=100.0
	   vmaxb=100.0
	   rka=10.d4
           rkb=10.d4
	   ca=0.5
	   cb=0.5
c
	   gga=vmaxa/rka
	   ggb=vmaxb/rkb
	   ttau_a=1.0/vmaxa
	   ttau_b=1.0/vmaxb
c          parametros con el tiempo normalizado a ttau_a 	
	   ga=1.0/rka
	   gb=1.0/rkb
           tau_a=1.0
	   tau_b=ttau_b/ttau_a
	   tau_ap=ca
	   tau_bp=cb*tau_b
	endif

c       parametros
	if(icaso.eq.3)then
   	   facg=0.05
	   factau=1.0
           ga=1.0*facg
	   gb=1.0*facg
	   tau_ap=0.5*factau
	   tau_bp=0.5*factau
	   tau_a=1.0*factau
	   tau_b=1.0*factau
	endif

c       masa de cada enzima, miligramos
c        xmg_a = 200.d0
c        xmg_b = 200.d0
c       calculo del numero de enzimas
c        xNA = xmg_a*xNAv/(xMMa*1.d3)  
c        xNB = xmg_b*xNAv/(xMMb*1.d3)  
c        print *, 'xNA=',xNA,' xNb=',xNB

c       calculo de la masa de cada enzima, miligramos
        xmg_a = xNA*xMMa*1.d3/xNAv
        xmg_b = xNB*xMMb*1.d3/xNAv
        print *,'A_mg, B_mg',xmg_a,xmg_b
       
c       calculo del punto fijo
        call pfix(S, xmg_a, xmg_b)

        print *, 'ga=',ga,' gb=',gb
        
        print *, 'mgA (mg) =',xmg_a,' mgB (mg) =',xmg_b
        
        print *,'Bs=',Bs
        print *,'m_beta1=',beta1
        print *,'m_beta2=',beta2
        print *,'m_beta3=',beta3

        print *,'As=',As
        print *,'m_alpha1 =',alpha1
        print *,'m_alpha2 =',alpha2
        print *,'m_alpha3 =',alpha3

        print *,'iflag = ',iflag

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        x15=0.0
c        x25=1000.0*tau_a
        x25=x25*tau_a
        print *,'tauA es ',tau_a
        print *,'El tiempo final de integracion es ',x25
        read(5,*)iparada
        if(iparada.eq.1) pause
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c-----------CONDICIONES INICIALES------------------------

      isig = 0

c     concentraciones iniciales cero
      if(isig.eq.0)then 
         ys5(1)=0.0
         ys5(2)=xmg_a
         ys5(3)=0.0
         ys5(4)=0.0
         ys5(5)=0.0
         ys5(6)=xmg_b
         ys5(7)=0.0
         ys5(8)=0.0
      endif

c     caso estacionario
      if(isig.ge.1)then 
         ys5(1)=As*(1.0d0 + 0.000d0)
         ys5(2)=alpha1
         ys5(3)=alpha2
         ys5(4)=alpha3
         ys5(5)=Bs*(1.0d0 - 0.000d0)
         ys5(6)=beta1
         ys5(7)=beta2
         ys5(8)=beta3
      endif

c--------------- ESCRITURA ENCABESADOS------------------------------
        write(6,*)'     --  PARAMETROS Y CONDICIONES INICIALES  --'
c        write(3,*)'     --  PARAMETROS Y CONDICIONES INICIALES  --'
        write(4,*)'     --  PARAMETROS Y CONDICIONES INICIALES  --'
        write(7,*)'     --  PARAMETROS Y CONDICIONES INICIALES  --'

        write (6,'(/1x,a,t30,6(1PE12.5,3x))') 
     #'ga,gb,tau_ap,tau_bp,tau_a,tau_b: ',
     #ga,gb,tau_ap,tau_bp,tau_a,tau_b
        write (6,'(/1x,a,t30,4(1PE12.5,3x),/)') 
     #'A  alpha1  alpha2   alpha3 (t=0): ',
     #ys5(1),ys5(2),ys5(3),ys5(4)
        write (6,'(/1x,a,t30,4(1PE12.5,3x),/)')
     #'B  beta1  beta2   beta3 (t=0): ',
     #ys5(5),ys5(6),ys5(7),ys5(8)

c        write (3,'(/1x,a,t30,6(1PE12.5,3x))')
c     #'ga,gb,tau_ap,tau_bp,tau_a,tau_b: ',
c     #ga,gb,tau_ap,tau_bp,tau_a,tau_b
c        write (3,'(/1x,a,t30,4(1PE12.5,3x),/)')
c     #'A  alpha1  alpha2   alpha3 (t=0): ',
c     #ys5(1),ys5(2),ys5(3),ys5(4)
c        write (3,'(/1x,a,t30,4(1PE12.5,3x),/)')
c     #'B  beta1  beta2   beta3 (t=0): ',
c     #ys5(5),ys5(6),ys5(7),ys5(8)

        write (4,'(/1x,a,t30,6(1PE12.5,3x))')
     #'ga,gb,tau_ap,tau_bp,tau_a,tau_b: ',
     #a,gb,tau_ap,tau_bp,tau_a,tau_b
        write (4,'(/1x,a,t30,4(1PE12.5,3x),/)')
     #'A  alpha1  alpha2   alpha3 (t=0): ',
     #ys5(1),ys5(2),ys5(3),ys5(4)
        write (4,'(/1x,a,t30,4(1PE12.5,3x),/)')
     #'B  beta1  beta2   beta3 (t=0): ',
     #ys5(5),ys5(6),ys5(7),ys5(8)

        write (7,'(/1x,a,t30,6(1PE12.5,3x))')
     #'ga,gb,tau_ap,tau_bp,tau_a,tau_b: ',
     #a,gb,tau_ap,tau_bp,tau_a,tau_b
        write (7,'(/1x,a,t30,4(1PE12.5,3x),/)')
     #'A  alpha1  alpha2   alpha3 (t=0): ',
     #ys5(1),ys5(2),ys5(3),ys5(4)
        write (7,'(/1x,a,t30,4(1PE12.5,3x),/)')
     #'B  beta1  beta2   beta3 (t=0): ',
     #ys5(5),ys5(6),ys5(7),ys5(8)


c-------------------------------------------------------------------
c-------------- Estados de equilibrio ----------------
        write(6,26)
c        write(3,26)
        write(4,26)
        write(7,26)
 26     format(t27,'ESTADO DE EQUILIBRIO :',/,5x,
     * 'ss        sfr    equilms   equilmc   equilmw   equilmh',
     *'   equilmse')
c        write(6,226)sskmy,sfr0kmy,emskmy,emckmy,emwkmy,emhkmy,emsekmy
c        write(3,226)sskmy,sfr0kmy,emskmy,emckmy,emwkmy,emhkmy,emsekmy
c        write(4,226)sskmy,sfr0kmy,emskmy,emckmy,emwkmy,emhkmy,emsekmy
c        write(7,226)sskmy,sfr0kmy,emskmy,emckmy,emwkmy,emhkmy,emsekmy
c 226    format(t1,7(1x,1pd9.2),/)

C......................................................
        a0=ys5(1)
        a00 = ys5(1) 
        alpha10=ys5(2)
        alpha100=alpha10
        alpha20=ys5(3)
        alpha30=ys5(4)
        b0=ys5(5)
        b00=ys5(5)
        beta10=ys5(6)
        beta100=beta10 
        beta20=ys5(7)
        beta30=ys5(8)
c--------------- ESCRITURA ENCABEZADOS------------------------------
      write(6,67)
c      write(3,67)
67	format(t8,
     #'k        t         A      alpha1   alpha2    alpha3'
     #,'        B      beta1   beta2    beta3    S')
c-------------------------------------------------------------------
c ojo        eps5=1.d-6
        eps5=1.d-6
c        h15=(x25-x15)/20000000.
        h15=(x25-x15)/100000.
        hmin5=0.0
        kmax5=20000002
c        dxsav5=(x25-x15)/20000000.
        dxsav5=tau_a/500.d0
        call ode5(ys5,nvar5,x15,x25,eps5,h15,hmin5,
     *  nok5,nbad5,deri5,rkqc)
        write(*,'(/1x,a,t30,i7)') 'pasos sucesivos: ',nok5
        write(*,'(1x,a,t30,i7)') 'pasos malos: ',nbad5
        write(*,'(/1x,a,t30,i10)') 'valores intermedios almacenados: ',
     3  kount5
c        write(3,'(/1x,a,t30,i6)') 'pasos sucesivos: ',nok5
c        write(3,'(1x,a,t30,i6)') 'pasos malos: ',nbad5
c        write(3,'(/1x,a,t30,i6)') 'valores intermedios almacenados: ',
     3  kount5
c       write (6,'(/1x,t9,a,t20,a,t33,a)') 't','Ms','Mc'
c       write (3,'(/1x,t9,a,t20,a,t33,a)') 't','Ms','Mc'
c       do 11 i=1,kount5
c       write(6,'(1x,f10.4,2x,2(1pe12.5,2x))') x5(i),y5(1,i),y5(2,i)
c       write(3,'(1x,f10.4,2x,2(1pe12.5,2x))') x5(i),y5(1,i),y5(2,i)
c  11   continue

C calculo de valores maximos y minimos de las variables
	xmax=0.0
	y1max=0.0
	y2max=0.0
	y3max=0.0
	y4max=0.0
	y5max=0.0
	y6max=0.0
	y7max=0.0
	y8max=0.0
	y1min=0.0
	y2min=0.0
	y3min=0.0
	y4min=0.0
	y5min=0.0
	y6min=0.0
	y7min=0.0
	y8min=0.0
	smax=0.0
	do n=1,kount5-1
	xmax=max(xmax,x5(n))
	y1max=max(y1max,y5(1,n))
	y2max=max(y2max,y5(2,n))
	y3max=max(y3max,y5(3,n))
	y4max=max(y4max,y5(4,n))
	y5max=max(y5max,y5(5,n))
	y6max=max(y6max,y5(6,n))
	y7max=max(y7max,y5(7,n))
	y8max=max(y8max,y5(8,n))
	smax=max(smax,y5(1,n)+y5(3,n)+y5(5,n)+y5(7,n))
	y1min=min(y1min,y5(1,n))
	y2min=min(y2min,y5(2,n))
	y3min=min(y3min,y5(3,n))
	y4min=min(y4min,y5(4,n))
	y5min=min(y5min,y5(5,n))
	y6min=min(y6min,y5(6,n))
	y7min=min(y7min,y5(7,n))
	y8min=min(y8min,y5(8,n))
	enddo
c	write (3,112)xmax
 112	format(t2,'tiempo maximo',f8.3)
c	write (3,113)y1min,y2min,y3min,y4min,y5min,y6min,y7min,y8min
 113	format(t2,'minimos',8f8.3)
c	write (3,114)y1max,y2max,y3max,y4max,y5max,y6max,y7max,y8max
 114	format(t2,'maximos',8f8.3)

c----------------- ESCRITURA DE ARCHIVO GRAFICADOR --------------------
	linefin=kount5-1+13
        open(unit=17,file='jm0.mon')
        write (17,215)
 215    format('erase',/)
        write (17,115)5.0*tau_a, 0.0,smax*1.1
 115    format('lw 2',/,'window 2 2 3',/,
     #'LIMITS 0 ',3(f10.3,1x),/,'BOX',/)
c encabezado
        write (17,150)smax*1.2,ga,gb,tau_ap
 150  format('relocate 7 ',f10.2,/,
     #'putlabel 4 ga=',f9.2,' gb=',f9.2,
     #' \\gt\\\\d\\gap\\\\u=',f6.3,'   ')
        write (17,151)smax*1.2,tau_a,tau_bp,tau_b
 151  format('relocate 7 ',f10.2,/,
     #'putlabel 6   \\gt\\\\d\\ga\\\\u=',f6.2,
     #' \\gt\\\\d\\gbp\\\\u=',f6.2,' \\gt\\\\d\\gb\\\\u=',f6.2)

        ltype=0
        nx=2
        ny=11
	linedown=13
	lineup=5000
	do lll=1,10
	if(linefin.lt.lineup)lineup=linefin
        write (17,116)ltype,linedown,lineup
        write (17,117)nx,ny
        write (17,118)
	if(lineup.eq.linefin)goto 333
	linedown=lineup+1
	lineup=lineup+5000-1
	enddo
 116  format('ltype ',i1,/,'lw 1',/,'DATA ret3.res',/,'LINES ',i5,1x,i5)
 117  format('XCOLUMN ',i2,/,'YCOLUMN ',i2)
 118  format('CONNECT',/)

 333	continue

c ejes
        write (17,119)
 119    format('ylabel S',/,'xlabel t')
        write (17,120)
 120    format('!.................')

c--------
        write (17,125)xmax, y5min,max(y3max,y7max)*1.2
 125    format('lw 2',/,'window 2 2 4',/,
     #'LIMITS 0 ',3(f10.3,1x),/,'BOX',/)

       ltype=0
        nx=2
        ny=5
        linedown=13
        lineup=5000
        do lll=1,10
        if(linefin.lt.lineup)lineup=linefin
        write (17,116)ltype,linedown,lineup
        write (17,117)nx,ny
        write (17,118)
        if(lineup.eq.linefin)goto 334
        linedown=lineup+1
        lineup=lineup+5000-1
        enddo
 334	continue

        ltype=2
        nx=2
        ny=9
        linedown=13
        lineup=5000
        do lll=1,10
        if(linefin.lt.lineup)lineup=linefin
        write (17,116)ltype,linedown,lineup
        write (17,117)nx,ny
        write (17,118)
        if(lineup.eq.linefin)goto 335
        linedown=lineup+1
        lineup=lineup+5000-1
        enddo
 335	continue


c ejes
        write (17,129)
 129	format('ylabel n\\\\d\\ga2\\\\u , n\\\\d\\gb2\\\\u',/,'xlabel t')
        write (17,120)

c--------

        write (17,135)y3max*1.2,0.0,y7max*1.2
 135    format('lw 2',/,'window 2 2 1',/,
     #'LIMITS 0 ',3(f10.3,1x),/,'BOX',/)
        ltype=0
        nx=5
        ny=9
        linedown=13
        lineup=5000
        do lll=1,10
        if(linefin.lt.lineup)lineup=linefin
        write (17,116)ltype,linedown,lineup
        write (17,117)nx,ny
        write (17,118)
        if(lineup.eq.linefin)goto 336
        linedown=lineup+1
        lineup=lineup+5000-1
        enddo
 336	continue

c ejes
        write (17,139)
 139    format('xlabel n\\\\d\\ga2\\\\u',/,'ylabel n\\\\d\\gb2\\\\u')
        write (17,120)
c-------
        write (17,145)y1max*1.2,0.0,y5max*1.2
 145    format('lw 2',/,'window 2 2 2',/,
     #'LIMITS 0 ',3(f10.3,1x),/,'BOX',/)
        ltype=0
        nx=3
        ny=7
        linedown=13
        lineup=5000
        do lll=1,10
        if(linefin.lt.lineup)lineup=linefin
        write (17,116)ltype,linedown,lineup
        write (17,117)nx,ny
        write (17,118)
        if(lineup.eq.linefin)goto 337
        linedown=lineup+1
        lineup=lineup+5000-1
        enddo
 337	continue

c ejes
        write (17,149)
 149    format('ylabel B',/,'xlabel A')
        write (17,120)

c-------

        write (17,121)
 121    format('expand 1.0',/,'end')
	close(17)
c-----------------------------------------------------------------------

        stop
        end
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
c			SUBRUTINAS	
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
c---------------------------------------------------------------
        subroutine retardo(t,r_apm_half,r_apm_one,r_bpm_half,r_bpm_one)
        implicit double precision (a-h,o-z)
        common/evovar/a0,alpha10,alpha20,alpha30,b0,beta10,beta20,beta30
        common/histprev/a00,alpha100,b00,beta100
        common/evoret/jlo
        common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
        common/path5/kmax5,kount5,dxsav5,x5(2000002),y5(10,2000002)
        common/pfijo/As,Bs,alpha1,alpha2,alpha3,beta1,beta2,beta3
c       r_apm_half=apm(t-tau_ap)
c       r_apm_one=apm(t-tau_a)
c       r_bpm_half=bpm(t-tau_ap)
c       r_bpm_one=bpm(t-tau_a)

c	print *,' -------------  en retardo t=',t,' -------------'
	call retpm(t-tau_ap,1,r_apm_half,dy)
c        print *,' t-tau_ap=',t-tau_ap,' r_apm_half',r_apm_half,' dy=',dy
	call retpm(t-tau_a,1,r_apm_one,dy)
c        print *,'    t-tau_a=',t-tau_a,' r_apm_one',r_apm_one
	call retpm(t-tau_bp,5,r_bpm_half,dy)
c        print *,'    t-tau_bp=',t-tau_bp,' r_bpm_half',r_bpm_half
	call retpm(t-tau_b,5,r_bpm_one,dy)
c        print *,'    t-tau_b=',t-tau_b,' r_bpm_half',r_bpm_one

        return
        end
c---------------------------------------------------------------
        subroutine retpm(tfi,kvar,fit,dy)
        implicit double precision (a-h,o-z)
        common/evovar/a0,alpha10,alpha20,alpha30,b0,beta10,beta20,beta30
        common/histprev/a00,alpha100,b00,beta100
        common/evoret/jlo
        common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
        common/path5/kmax5,kount5,dxsav5,x5(2000002),y5(10,2000002)
        common/pfijo/As,Bs,alpha1,alpha2,alpha3,beta1,beta2,beta3
        common/otrospars/Vsim,xMMa,xMMb,xNAv
        common/sust/S,xmg_a
        real*8 coorx(4),coory(4)
c       print *,'jlo antes de hunt=',jlo

	if(kvar.eq.1) gg=ga
	if(kvar.eq.5) gg=gb

        call hunt(x5,kount5,tfi,jlo)
c        print *,'kount5, jlo DESPUES de hunt=',kount5,jlo
cc        if(kount5.lt.4)fit=Ms0
cc        if(kount5.lt.4)return
	if(tfi.le.0.0.and.kvar.eq.1)then
        fit=ga*a00*alpha100
	dy=0.0
	return
	endif

        if(tfi.le.0.0.and.kvar.eq.5)then
        fit=gb*b00*beta100
	dy=0.0
        return
        endif

        if(kount5.lt.4)print *,'problema en retapmhalf'
        if(kount5.lt.4)pause
        j=jlo
        np=kount5

	ret1=gg*y5(kvar,1)*y5(kvar+1,1)
	ret2=gg*y5(kvar,2)*y5(kvar+1,2)
	ret3=gg*y5(kvar,3)*y5(kvar+1,3)
	ret4=gg*y5(kvar,4)*y5(kvar+1,4)
        if (j .eq. 0)  then
           print *,'en retma extrapolation    j=0'
           coorx(1)=x5(1)
           coorx(2)=x5(2)
           coorx(3)=x5(3)
           coorx(4)=x5(4)
           coory(1)=ret1
           coory(2)=ret2
           coory(3)=ret3
           coory(4)=ret4
        else if (j .le. 2.) then
           coorx(1)=x5(1)
           coorx(2)=x5(2)
           coorx(3)=x5(3)
           coorx(4)=x5(4)
           coory(1)=ret1
           coory(2)=ret2
           coory(3)=ret3
           coory(4)=ret4
        else if(j .eq. np) then
           print *,'en retma extrapolation    j=np'
           coorx(1)=x5(np-3)
           coorx(2)=x5(np-2)
           coorx(3)=x5(np-1)
           coorx(4)=x5(np)
           coory(1)=gg*y5(kvar,np-3)*y5(kvar+1,np-3)
           coory(2)=gg*y5(kvar,np-2)*y5(kvar+1,np-2)
           coory(3)=gg*y5(kvar,np-1)*y5(kvar+1,np-1)
           coory(4)=gg*y5(kvar,np)*y5(kvar+1,np)
        else if (j .ge. np-1) then
           coorx(1)=x5(np-3)
           coorx(2)=x5(np-2)
           coorx(3)=x5(np-1)
           coorx(4)=x5(np)
           coory(1)=gg*y5(kvar,np-3)*y5(kvar+1,np-3)
           coory(2)=gg*y5(kvar,np-2)*y5(kvar+1,np-2)
           coory(3)=gg*y5(kvar,np-1)*y5(kvar+1,np-1)
           coory(4)=gg*y5(kvar,np)*y5(kvar+1,np)
        else
           coorx(1)=x5(j-1)
           coorx(2)=x5(j)
           coorx(3)=x5(j+1)
           coorx(4)=x5(j+2)
           coory(1)=gg*y5(kvar,j-1)*y5(kvar+1,j-1)
           coory(2)=gg*y5(kvar,j)*y5(kvar+1,j)
           coory(3)=gg*y5(kvar,j+1)*y5(kvar+1,j+1)
           coory(4)=gg*y5(kvar,j+2)*y5(kvar+1,j+2)
        end if
c       print *,'coorx',coorx
c       print *,'coory',coory
        call polint(coorx,coory,4,tfi,fit,dy)
c        if(tfi.gt.1.3*tau_a.and.tfi.lt.1.31*tau_a) then
c        if(kvar.eq.1) then
c           print *,'jlo, tfi, coorx/tau_a, coory, fit', jlo, tfi/tau_a, 
c     #coorx(1)/tau_a, 
c     #coorx(2)/tau_a, coorx(3)/tau_a, coorx(4)/tau_a, coory, fit
c           print *,'y5(kvar,jlo), y5(kvar+1,jlo), gg , producto',  
c     #y5(kvar,jlo), y5(kvar+1,jlo),gg,gg*y5(kvar,jlo)*y5(kvar+1,jlo) 
c        pause
c        if(jlo.gt.400.and.jlo.lt.500) pause
c        endif
c        endif
       if(dabs(dy/fit).gt.0.05)print *,'en retma tfi,fit,dy=',tfi,fit,dy
       if(fit.lt.0.0)print *,'//en retma/// INTERPOLACION < 0  fit=',fit
        if(fit.lt.0.0)print *,'   INTERPOLACION < 0  fit puesta a 0 en
     #kvar =',kvar
        if(fit.lt.0.0)fit=0.0
        return
        end
c---------------------------------------------------------------
        subroutine deri5(x,y,dydx)
        implicit double precision (a-h,o-z)
        common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
        common/histprev/a00,alpha100,b00,beta100
        dimension y(8),dydx(8)
        common/signal/isig,iflag
        common/pfijo/As,Bs,alpha1,alpha2,alpha3,beta1,beta2,beta3
        common/otrospars/Vsim,xMMa,xMMb,xNAv
        common/sust/S,xmg_a

C==========================================================
c VARIABLES
c                   y1==A
c                   y2==alpha1
c                   y3==alpha2
c                   y4==alpha3
c                   y5==B
c                   y6==beta1
c                   y7==beta2
c                   y8==beta3
c  8 ec.  dif. ord.
c                 no son estas
c                   dy1/dt=r_bpm_half-ga*y1*y2
c                   dy2/dt=r_apm_one-ga*y1*y2
c                   dy3/dt=ga*y1*y2-r_apm_half
c                   dy4/dt=r_apm_half-r_apm_one
c
c                   dy5/dt=r_apm_half-gb*y5*y6
c                   dy6/dt=r_bpm_one-gb*y5*y6
c                   dy7/dt=gb*y5*y6-r_bpm_half
c                   dy8/dt=r_bpm_half-r_bpm_one

crrrrrrrrrrrrrrrrrrrrr  ret= retraso
        call retardo(x,r_apm_half,r_apm_one,r_bpm_half,r_bpm_one)
c-----------------------------------------------------------
c
	y1=y(1)
	y2=y(2)
	y3=y(3)
	y4=y(4)
	y5=y(5)
	y6=y(6)
	y7=y(7)
	y8=y(8)
c	if(y1.lt.0.0)y1=0.0
c	if(y2.lt.0.0)y2=0.0
c	if(y3.lt.0.0)y3=0.0
c	if(y4.lt.0.0)y4=0.0
c	if(y5.lt.0.0)y5=0.0
c	if(y6.lt.0.0)y6=0.0
c	if(y7.lt.0.0)y7=0.0
c	if(y8.lt.0.0)y8=0.0
c se agrega A durante un lapso de tiempo
 	x0=2.0*tau_a
	sigma=0.1*tau_a**2.0d0
c	r=(800.0/112.0)*alpha100/(tau_a)

c	r=3.211425*alpha100/(tau_a)
	r=(200.d0/xmg_a)*(S/180.)*3.211425*alpha100/(tau_a)
c	r=(200.d0/xmg_a)*(S/180.)*3.223925*alpha100/(tau_a)
   
        if(isig.eq.0) dadt_externo=r*dexp((-(x-x0)**2.0d0)/sigma)
        if(isig.eq.1) dadt_externo = 0.0
        if(isig.eq.2) dadt_externo = 10.*dexp((-(x-x0)**2.0d0)/sigma)
        
        dydx(1)=(r_bpm_half/Vsim)-(ga/Vsim)*y1*y2+0.5*dadt_externo
        dydx(2)=r_apm_one-ga*y1*y2
        dydx(3)=ga*y1*y2-r_apm_half
        dydx(4)=r_apm_half-r_apm_one
        dydx(5)=(r_apm_half/Vsim)-(gb/Vsim)*y5*y6
        dydx(6)=r_bpm_one-gb*y5*y6
        dydx(7)=gb*y5*y6-r_bpm_half
        dydx(8)=r_bpm_half-r_bpm_one

        dydx(2)=dydx(2)*1.d-3*xMMa
        dydx(3)=dydx(3)*1.d-3*xMMa
        dydx(4)=dydx(4)*1.d-3*xMMa

        dydx(6)=dydx(6)*1.d-3*xMMb
        dydx(7)=dydx(7)*1.d-3*xMMb
        dydx(8)=dydx(8)*1.d-3*xMMb

c	if(y(1).le.0.0.and.dydx(1).lt.0.0)dydx(1)=0.0
c	if(y(2).le.0.0.and.dydx(2).lt.0.0)dydx(2)=0.0
c	if(y(3).le.0.0.and.dydx(3).lt.0.0)dydx(3)=0.0
c	if(y(4).le.0.0.and.dydx(4).lt.0.0)dydx(4)=0.0
c	if(y(5).le.0.0.and.dydx(5).lt.0.0)dydx(5)=0.0
c	if(y(6).le.0.0.and.dydx(6).lt.0.0)dydx(6)=0.0
c	if(y(7).le.0.0.and.dydx(7).lt.0.0)dydx(7)=0.0
c	if(y(8).le.0.0.and.dydx(8).lt.0.0)dydx(8)=0.0
c
c-------------------------------------------------------------------
c        if(x.gt.20.) then
c	print *,'en derivs x=',x,'  y=',y
c	print *,'en derivs ga=',ga,'  gb=',gb
c	print *,'en derivs tau_a=',tau_a,'  tau_b=',tau_b
c        pause
c        endif
c        print *,' dydx=',dydx
c	print *,' r_apm_half,r_apm_one,r_bpm_half,r_bpm_one',
c     #r_apm_half,r_apm_one,r_bpm_half,r_bpm_one
c-------------------------------------------------------------------
        return
        end
c-------------------------------------------------
      SUBROUTINE HUNT(XX,N,X,JLO)
      implicit double precision (a-h,o-z)
      DIMENSION XX(N)
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)goto 4
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
 4      return
      END
C----------------------------------------------------------------
      SUBROUTINE ODE5(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RK
     *QC)
      implicit double precision (a-h,o-z)
c      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.D-30)
      PARAMETER (MAXSTP=20000002,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.D-20)
      EXTERNAL DERIVS
      COMMON /PATH5/ KMAX5,KOUNT5,DXSAV5,XP5(2000002),YP5(10,2000002)
        common/evovar/a0,alpha10,alpha20,alpha30,b0,beta10,beta20,beta30
        common/histprev/a00,alpha100,b00,beta100
        common/evoret/jlo
        common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
        common/sust/S,xmg_a
        common/otrospars/Vsim,xMMa,xMMb,xNAv

      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT5=0

      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE

      XSAV=X-DXSAV5*TWO
      
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        
        DO 12 I=1,NVAR
          YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
c          YSCAL(I)=S
12      CONTINUE

c       Guarda resultados en xp5, yp5
        IF(KMAX5.GT.0)THEN
          IF(DABS(X-XSAV).GT.DABS(DXSAV5)) THEN
            IF(KOUNT5.LT.KMAX5-1)THEN
              KOUNT5=KOUNT5+1
              XP5(KOUNT5)=X
              DO 13 I=1,NVAR
                YP5(I,KOUNT5)=Y(I)
13            CONTINUE
              XSAV=X
c-------------------------- ESCRITURA ------------------------------
c              ss=y(1)+y(5)+(y(3)+y(7))*1.d3/(xMMa*Vsim)
              ss1=y(1)+y(3)*1.d3/(xMMa*Vsim)
              ss2=y(5)+y(7)*1.d3/(xMMa*Vsim)
              tot=ss1+ss2
              
              fa1=y(2)/xmg_a
              fa2=y(3)/xmg_a
              fa3=y(4)/xmg_a
              fb1=y(6)/xmg_a
              fb2=y(7)/xmg_a
              fb3=y(8)/xmg_a

              pa=a1*(y(1)/400.)
              pb=y(1)/400.
              bt=y(5)/400.
              dym1=DYDX(2)/(1000.*xmg_a)
c               ss=y(2)+y(3)+y(4)
c              if(X.gt.X2/2.)then
              write(6,'(5x,i7,f10.4,1x,9(1pe9.2,1x))') kount5,x/tau_a
     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),ss
c              write(3,'(5x,i7,f10.4,1x,9(1pd20.14,1x))') kount5,x/tau_a
c     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),ss
              write(3,'(5x,i7,f10.4,1x,20(1pd20.13,1x))') kount5,x/tau_a
     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),fa1,fa2,fa3,fb1,fb2,fb3
     1,pa,ss1,pb,bt,dym1,tot
              write(13,'(5x,i7,f10.4,1x,8(1pd20.13,1x))') kount5,x/tau_a
     1,DYDX(1),DYDX(2),DYDX(3),DYDX(4),DYDX(5),DYDX(6),DYDX(7),DYDX(8)
c              endif
              if(x/tau_a.gt.1000..and.x/tau_a.lt.1050.) then
       write(35,'(5x,i7,f10.4,1x,9(1pd20.14,1x))') kount5,x/tau_a
     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),ss
              endif 
              if(x/tau_a.gt.2000..and.x/tau_a.lt.2050.) then
       write(36,'(5x,i7,f10.4,1x,9(1pd20.14,1x))') kount5,x/tau_a
     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),ss
              endif
c-------------------------------------------------------------------
            ENDIF
          ENDIF
        ENDIF
        
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
c              write(44,'(5x,i5,f10.4,1x,9(1pe15.7,1x))') kount5,x
c     1,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),ss
        ELSE
          NBAD=NBAD+1
        ENDIF
        
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX5.NE.0)THEN
            KOUNT5=KOUNT5+1
            XP5(KOUNT5)=X
            DO 15 I=1,NVAR
              YP5(I,KOUNT5)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        
        IF(DABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum.'
        H=HNEXT
c	print *,'H en ODE5  ',H
16    CONTINUE
      PAUSE 'Too many steps.'
      write(6,*) kount5
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=8,FCOR=.0666666667,
     *    ONE=1.,SAFETY=0.9,ERRCON=6.D-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV)PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,DABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
c	print *,'XSAV,H,XSAV+HH,XSAV+H,ERRMAX',XSAV,H,XSAV+HH,XSAV+H,ERRMAX
c	print *,'YSCAL(I)',YSCAL
c        print *,'YTEMP(I)',YTEMP
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C------------------------------------------------------------------      
      subroutine pfix(S, xMA, xMB)
      implicit double precision (a-h,o-z)
      complex ci,s1,s2,s3,minus_one,w,y1,y2,y3
      common/param/ga,gb,tau_ap,tau_bp,tau_a,tau_b
      common/pfijo/As,Bs,alpha1,alpha2,alpha3,beta1,beta2,beta3
      common/signal/isig,iflag
      common/otrospars/Vsim,xMMa,xMMb,xNAv
      common/vels/vmaxa,vmaxb,rka,rkb,ca,cb,gga,ggb

       minus_one=-1.0  
       ci=csqrt(minus_one)
       pi = 4*atan(1.)

c      Determinacion de parametros auxiliares
       va = vmaxa
       vb = vmaxb

       tau_as = 1.d0/vmaxa
       tau_aps = ca*tau_as
       tau_bs = 1.d0/vmaxb
       tau_bps = cb*tau_bs

c      Determinacion de los coeficientes de la cubica
       if(va.ge.vb) then
c         solucion para B       
          a3=Vsim*xMMb*1.d-6*tau_bs*(xMMa*tau_as-(xMA/xMB)*xMMb*tau_bs)
          a2=1.d-3*((Vsim/ggb)*xMMa*tau_as+(xMA/xMB)*Vsim*(xMMb**2.d0)*
     #1.d-3*(tau_bs**2.d0)*S+(tau_aps+tau_bps)*(xMMa*tau_as*xMB-xMMb*
     #tau_bs*xMA)-2*(xMA*Vsim/(xMB*ggb))*xMMb*tau_bs-Vsim*xMMb*tau_bs*
     #(1./gga+xMMa*1.d-3*tau_as*S))
          a1=2*(xMA*Vsim/(xMB*ggb))*xMMb*1.d-3*tau_bs*S-(Vsim/ggb)*
     #(1./gga+xMMa*1.d-3*tau_as*S)-(xMA/ggb)*(tau_aps+tau_bps+Vsim/(xMB
     #*ggb))
          a0=xMA*Vsim*S/(xMB*(ggb**2.d0))
       else if(va.lt.vb) then
c         solucion para A
          a3=Vsim*xMMa*1.d-6*tau_as*(xMMb*tau_bs-(xMB/xMA)*xMMa*tau_as)
          a2=1.d-3*((Vsim/gga)*xMMb*tau_bs+(xMB/xMA)*Vsim*(xMMa**2.d0)*
     #1.d-3*(tau_as**2.d0)*S+(tau_aps+tau_bps)*(xMMb*tau_bs*xMA-xMMa*
     #tau_as*xMB)-2*(xMB*Vsim/(xMA*gga))*xMMa*tau_as-Vsim*xMMa*tau_as*
     #(1./ggb+xMMb*1.d-3*tau_bs*S))
          a1=2*(xMB*Vsim/(xMA*gga))*xMMa*1.d-3*tau_as*S-(Vsim/gga)*
     #(1./ggb+xMMb*1.d-3*tau_bs*S)-(xMB/gga)*(tau_aps+tau_bps+Vsim/(xMA
     #*gga))
          a0=xMB*Vsim*S/(xMA*(gga**2.d0))
       endif    

       print *,'a3,a2,a1,a0',a3,a2,a1,a0

c      Calculo del discriminante
       delta = 4.0*(a1**3.d0)*a3-(a1**2.d0)*(a2**2.d0)+4.0*a0*(a2**3.d0)
     #-18.0*a0*a1*a2*a3+27.0*(a0**2.d0)*(a3**2.d0) 

c      Las raices son reales
       if(delta.lt.0.0) then
       
c        Hay dos soluciones
          if(a3.eq.0.0) then

             x1 = (-a1 + sqrt((a1**2.d0) - 4.0*a2*a0))/(2.0*a2)
             x2 = (-a1 - sqrt((a1**2.d0) - 4.0*a2*a0))/(2.0*a2);

	     if(x1.eq.0.0.and.x1.lt.S) then
	        Bs = x1;
	     else if(x2.ge.0.0.and.x2.lt.S) then
	        Bs = x2;
             endif   

             print *,'x1,x2',x1,x2

	     beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	     beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3*tau_bs)
	     beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)

             As=Bs
             alpha1=beta1
             alpha2=beta2
             alpha3=beta3

c         Hay tres soluciones
          else if(a3.ne.0.0) then

	     e = (a1-a2*a2/(3.*a3))/a3
	     f = (a0+2.*(a2*a2*a2)/(27.*(a3*a3))-(a2*a1/(3.*a3)))/a3

	     sc = -e/3.
	     e = -(e*e*e)/27.

	     w = 0.5*(-f+csqrt(cmplx(f**2.d0-4*e)));

             s1r = real(w)
             s1i = aimag(w)
             s1mag = dsqrt(s1r**2.d0+s1i**2.d0)
             s1mag = s1mag**(1.d0/3.d0)
             s1phas = datan(s1i/s1r)
             
             s1 = cmplx(s1mag*dcos(s1phas/3.0),s1mag*dsin(s1phas/3.0))
             s2 = cmplx(s1mag*dcos((s1phas+2.*pi)/3.),s1mag*dsin((s1phas
     #+2.*pi)/3.))
             s3 = cmplx(s1mag*dcos((s1phas+4.*pi)/3.),s1mag*dsin((s1phas
     #+4.*pi)/3.))

             y1 = s1 + cmplx(sc)/s1
             y2 = s2 + cmplx(sc)/s2
             y3 = s3 + cmplx(sc)/s3
             
             x1 = real(y1 - cmplx(a2/(3.*a3)))
             x2 = real(y2 - cmplx(a2/(3.*a3)))
             x3 = real(y3 - cmplx(a2/(3.*a3)))

             print *,'x1,x2,x3',x1,x2,x3

c        Selecciona la solucion con sentido fisico
	     iflag = 0

             if(x1.gt.0.0.and.x1.lt.S) then
                if(va.ge.vb) then
	           As = S - x1 - (xMB/Vsim)*(tau_aps + tau_bps)/(1./(ggb*x1)
     #+xMMb*1.d-3*tau_bs)
	           if(As.gt.0.0) then
	              Bs = x1

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        else
	           Bs = S - x1 - (xMA/Vsim)*(tau_aps + tau_bps)/(1./(gga*x1)
     #+xMMa*1.d-3*tau_as)
	           if(Bs.gt.0.0) then
	              As = x1

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        endif   
	     endif   
	 
             if(x2.gt.0.0.and.x2.lt.S) then
                if(va.ge.vb) then
	           As = S - x2 - (xMB/Vsim)*(tau_aps + tau_bps)/(1./(ggb*x2)
     #+xMMb*1.d-3*tau_bs)
	           if(As.gt.0.0) then
	              Bs = x2

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        else
	           Bs = S - x2 - (xMA/Vsim)*(tau_aps + tau_bps)/(1./(gga*x2)
     #+xMMa*1.d-3*tau_as)
	           if(Bs.gt.0.0) then
	              As = x2

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        endif   
	     endif   
             
             if(x3.gt.0.0.and.x3.lt.S) then
                if(va.ge.vb) then
	           As = S - x3 - (xMB/Vsim)*(tau_aps + tau_bps)/(1./(ggb*x3)
     #+xMMb*1.d-3*tau_bs)
	           if(As.gt.0.0) then
	              Bs = x3

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        else
	           Bs = S - x3 - (xMA/Vsim)*(tau_aps + tau_bps)/(1./(gga*x3)
     #+xMMa*1.d-3*tau_as)
	           if(Bs.gt.0.0) then
	              As = x3

	              beta1=xMB/(1. + xMMb*1.d-3*Bs/rkb)
	              beta2=xMMb*1.d-3*tau_bps*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              beta3=xMMb*1.d-3*(1.-cb)*tau_bs*xMB/(1./(ggb*Bs)+xMMb*1.d-3
     #*tau_bs)
	              alpha1=xMA/(1. + xMMa*1.d-3*As/rka)
	              alpha2=xMMa*1.d-3*tau_aps*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              alpha3=xMMa*1.d-3*(1.-ca)*tau_as*xMA/(1./(gga*As)+xMMa*1.d-3
     #*tau_as)
	              delta = As + Bs + (1.d3/Vsim)*(alpha2/xMMa + beta2/xMMb) 

	              if(abs(delta - S).lt.1.d-3) iflag = 1
	           endif   
	        endif   
	     endif   
          endif
       endif
       end

