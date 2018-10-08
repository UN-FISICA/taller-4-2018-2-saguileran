#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#from sympy import *
from scipy import optimize
import numpy as np
import math
#from math import *

class Derivada:
    def __init__(self, f, metodo, dx=1e-1):
        self.f=f
        self.dx=float(dx)
        self.metodo=metodo
    #def __call__(self, x):
     #   f, dx = self.f, self.dx      # make short forms
      #  return (f(x+dx) - f(x))/dx
    def calc(self,x):        
        if self.metodo=='adelante':
            return((-float(self.f(x))+float(self.f(x+self.dx)))/self.dx)
        elif self.metodo=='central':
            return((float(self.f(x+(self.dx)/2))-float(self.f(x-(self.dx)/2)))/(self.dx))
        elif self.metodo=='extrapolada':
            dca=(float(self.f(x+self.dx/4))-float(self.f(x-self.dx/4)))/(self.dx/2)
            dcn=(float(self.f(x+self.dx/2))-float(self.f(x-self.dx/2)))/self.dx
            return((4*dca-dcn)/3)
        elif self.metodo=='segunda':
            #Segunda derivada
            return(float(self.f(x+self.dx)+self.f(x-self.dx)-2*self.f(x))/(self.dx**2))
        else:return('Método no válido')
  ##Núero de ceros en dx, a izquierda, igual a numero de digitos aproximados
        """ Método que retorna el valor numérico de la derivada de la
    función f evaluada en el punto x"""
#Mejorar metodos
#x = Symbol('x')

class Zeros:
    def __init__(self, f, metodo, error=1e-6, max_iter=100):
        self.f=f
        self.metodo=metodo
        self.error=float(error)
        self.max_inter=int(max_iter)

    def zero(self,vi): 
        #Newton
        def Newton(vi):
            i,vi=0,float(vi)
            while abs(float(self.f(vi)))>self.error:
                F=Derivada(self.f,'extrapolada',1e-6)
                vi=vi-(float(self.f(vi))/F.calc(vi))
                i+=1
                if i==self.max_inter: break
                elif float(self.f(vi))==0:break
            return(vi)
        #Bisectriz
        def Bisection(x):
            a,b=x[0],x[1]
            if float(self.f(b))*float(self.f(a))<0:
                a,b=float(vi[0]),float(vi[1])
                c,i=(a+b)/2,0
                while abs(float(self.f(c)))>self.error:
                    if float(self.f(c))==0:break
                    elif i==self.max_inter: break
                    elif float(self.f(a))*float(self.f(c))<0:b=c
                    else: a=c
                    c=(a+b)/2
                    i+=1
                return(c)
            else: 'Datos iniciales no apropiados'
        #Interpolación lineal
        def Interpol(vi):
                a,b=float(vi[0]),float(vi[1])
                c,i=(float(self.f(b))*a-float(self.f(a))*b)/(float(self.f(b))-float(self.f(a))),0
                while abs(float(self.f(c)))>self.error:
                    if float(self.f(c))==0:break
                    elif i==self.max_inter: break
                    elif float(self.f(a))*float(self.f(c))<0:b=c
                    else: a=c
                    c=(float(self.f(b))*a-float(self.f(a))*b)/(float(self.f(b))-float(self.f(a)))
                    i+=1
                return(c)
        #b=float(vi[0]),float(vi[1])
        if self.metodo=='newton':return((Newton(vi)))
        elif self.metodo=='bisectriz':return((Bisection(vi)))
        elif self.metodo=='interpolacion':return((Interpol(vi)))
        elif self.metodo=='newton-sp':
            return((optimize.newton(self.f,vi,fprime=None, tol=self.error, maxiter=int(self.max_inter), fprime2=None))) #Calcula la derivada con el método de la secante
        elif self.metodo=='fsolve-sp':
            F= lambda x : x**3
            return(optimize.fsolve(F,vi,xtol=self.error,maxfev=self.max_inter))
            #return(print(optimize.fsolve(self.f, np.array(vi), args=(), fprime=None, full_output=0, col_deriv=0, xtol=self.error, maxfev=int(self.max_inter))))
        elif self.metodo=='brentq-sp':
            if self.f(vi[0])*self.f(vi[1])<0:
                return(optimize.brentq(self.f,vi[0],vi[1],xtol=self.error,maxiter=self.max_inter))
            else: print('Valores iniciales no válidos')
        else: print('Método no válido')
    
if __name__ ==" __main__ ": pass
else:
    print(Derivada(exp,'segunda').calc(1))
    print(Derivada(exp,'central').calc(1))
    print(Derivada(exp,'extrapolada').calc(1))   
    print(Derivada(exp,'segunda').calc(1))
    print()
    print(Zeros(cos,'newton').zero(0.22))
    print(Zeros(cos,'bisectriz').zero((0,3.22)))   
    print(Zeros(cos,'interpolacion').zero((0.22,1.22)))
    print(Zeros(cos,'newton-sp').zero(5.23))
    print(Zeros(sin,'brentq-sp').zero((-1.2,2)))
    print(Zeros(sin,'fsolve-sp').zero((-1.2,2)))