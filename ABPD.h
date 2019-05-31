double am(double v)
{
	return(     ((127.0/105.0*v)+(201.0/7.0))  / ( 10.0-(10.0*exp((-201.0/70.0)-(127.0/1050.0*v))  )    )         ) ;
}

double bm(double v)
{
	return(4.0*  exp( (-188.0/63.0)-(127.0/1890.0*v) )         );
}

double mbar(double v)
{
	return( am(v)/(am(v)+bm(v) ) );
}


double ah(double v)
{
        return(7.0/100.0 *  exp( (-94.0/35.0)-(127.0/2100.0*v) )         );
}

double bh(double v)
{
        return(1.0/ (1.0+ exp( (-83.0/35.0)-(127.0/1050.0*v))   )         );
}
double hbar(double v)
{
	return( ah(v)/(ah(v)+bh(v) ) );
}


double an(double v)
{
        return(     ((127.0/105.0*v)+(166.0/7.0))  / ( 100.0-(100.0*exp((-83.0/35.0)-(127.0/1050.0*v))  )    )         ) ;
}

double bn(double v)
{
        return(1.0/8.0 *  exp( (-59.0/140.0)-(127.0/8400.0*v) )         );
}

double nbar(double v)
{
	return( an(v)/(an(v)+bn(v) ) );
}


double ma(double v)
{
        return(1.0/ (1.0+ exp( (v-va)/sa  )  )         );
}


double hai(double v)
{
        return(1.0/ (1.0+ exp( (v-vb)/sb  )  )         );
}



double zv(double v)
{
        return(1.0/ (1.0+ exp(-slope_ca_act/100.0* (v-zb) )  )         );
}



double psbar(double v)
{
        return ( 1 / ( exp(-2.0*(v+45.0)) + 1 ) ); //note its v not vs//0.15//-50
} 
double pstau(double v)
{
	return ((tps/(1+exp (-(v+50)/psst) )  ) + tpsmin ) ;
}

double bbar(double v)
{
        return ( 1 / ( exp(-2.0*(v+42.0)) + 1 ) ); //note its v not vs//0.15//-50
} 

double d_c (double v, double c, double z)
{
	//fprintf(rate_ca,"%lf\t%lf\n",c,1000*(rho*  (      ( (kca*z*(vca-v))/(1.0+2.0*c) ) - c      )    ));
	return(rho*  (      ( (kca*z*(vca-v))/(1.0+2.0*c) ) - c      )    );
}


double d_n (double v, double n )
{
	return(lambdan * (an(v) *(1.0-n) - bn(v) * n)       );
}

double d_h (double v, double h )
{
        return(lambdah * (ah(v) *(1.0-h) - bh(v) * h)       );
}


double d_z (double v, double z)
{
	return(  (zv(v)-z )/tauz );
}


double d_ha (double v, double ha)
{
        return(  (hai(v)-ha )*ka );
}

double f_ps(double v, double ps)
{
        return ( (psbar(v) - ps)/pstau(v));
}
double f_kf(double v, double b)
{
        return ( (bbar(v) - b)/btau);
}


double ina(double v, double h)
{
	return(gna * pow(mbar(v),3)*h*(v-vna) );
}


double ica(double v, double c,double z )
{
	tm_ica=gca * z/(ca_inact_var + c) *(v-vca);
        return(gca * z/(ca_inact_var + c) *(v-vca) ); //0.5
}

double ik(double v, double n )
{
	return(gk * pow(n,4)*(v-vk) );
}



double ikca(double v, double c)
{
	 tm_ikca=gkca * c/(0.5+c) * (v-vk) ;
	return (    gkca * c/(0.5+c) * (v-vk)     );
}

double ia(double v, double ha )
{
        return(ga * pow(ma(v),3)*ha * (v-vk) );
}

double i_ps(v,ps)
double v,ps;
{
	//t_ips=fopen("tm_ips.xls","w");
	//ips= G_PS * ps * (v-vk);
	//fprintf(t_ips,"%lf\t%lf\n",time_,ips);
        return (G_PS * ps * (v-vk));
}

double i_kf(v,b)
double v,b;
{
	return (G_kf * b * (v-vk));
}

double isyn(double v )
{
        return (gsyn * (v-vsyn) );
}


double il(double v )
{
        return (gl * (v-vl) );
}
double il_s(double v )
{
        return (gl_s * (v-vl) );
}
double il_pd(double v )
{
        return (gl_pd * (v-vl) );
}
double il_a(double v )
{
        return (gl_a * (v-vl) );
}

double f_vs(double vs, double vpd, double Va)
{
	
	return (iext_s - (isyn(vs) +  il_s(vs) + Gspd * (vs-vpd) + Gsa * (vs-Va)) /CMs);
}

double f_vpd(double vpd, double vs,double vd,double Va,double n, double h,double ps,double mk)
{
	
	return ( - ( il_pd(vpd) + i_ps(vpd,ps) - Gdpd * (vd - vpd) - Gapd* (Va-vpd) - Gspd * (vs-vpd)) /CMd);
}

double f_vd(double vpd, double vd,double c, double z,double ha,double Va, double b)
{
	//imk= i_mk(vd,mk);
	return( iext_d -( i_kf(vd,b) + ica(vd,c,z) + ikca(vd,c) + ia(vd,ha) + il(vd) + Gdpd * (vd - vpd) + Gda*(vd-Va) ) /CMd);
}

double f_va(double vpd, double Va, double n, double h, double vd, double vs)
{
	return( iext - ( ina(Va,h) + ik(Va,n) + il_a(Va) + Gapd* (Va-vpd) - Gda*(vd-Va)-Gsa * (vs-Va) ) /CMa);
}

void deriv_()
{
        // first neuron
        deriv[0] = f_vs(arg[0],arg[9],arg[7]);
	deriv[9] = f_vpd(arg[9],arg[0],arg[6],arg[7],arg[2],arg[3],arg[8],arg[10]);
	deriv[6] = f_vd(arg[9],arg[6],arg[1],arg[4],arg[5],arg[7], arg[10]);
	deriv[7] = f_va(arg[9],arg[7],arg[2],arg[3],arg[6],arg[0]);
	
	deriv[2] = d_n(arg[7],arg[2]);
	deriv[3] = d_h(arg[7],arg[3]);
        
	deriv[1] = d_c(arg[6],arg[1],arg[4]);
        deriv[4] = d_z(arg[6],arg[4]);
	deriv[5] = d_ha(arg[6],arg[5]);
	
	deriv[8] = f_ps(arg[9],arg[8]);
	deriv[10] = f_kf(arg[6],arg[10]);
}

void rkm()
{
        double old[N], epsilon[N];
        double k1[N], k2[N], k3[N], k4[N],k5[N];
        int i;

        start:
        for (i=0; i < N; i++) {
                arg[i]=state[i];
                //printf("%f\t",state[i]);
                }
        deriv_();
        for (i=0 ;i < N; i++)
        {
                old[i] = state[i];
                k1[i] = step*deriv[i]/3;
                arg[i] = state[i] + k1[i];
        }

        deriv_();
        for (i=0 ;i < N; i++)
        {
                k2[i] = step*deriv[i]/3;
                arg[i] = state[i] + k1[i]/2 +k2[i]/2;
        }

        deriv_();
        for (i=0; i < N; i++)
        {
                k3[i] = step*deriv[i]/3;
                arg[i] = state[i] + 3*k1[i]/8 +9*k3[i]/8;
        }

        deriv_();
        for (i=0; i < N; i++)
        {
                k4[i] = step*deriv[i]/3;
                arg[i] = state[i] + 3*k1[i]/2 - 9*k3[i]/2 + 6*k4[i];
        }

        deriv_();
        for (i=0; i < N; i++)
        {
                k5[i] = step*deriv[i]/3;
                state[i] = state[i] + (k1[i] +4*k4[i] +k5[i])/2 ;
                epsilon[i] = (k1[i] - 9*k3[i]/2 + 4*k4[i] - k5[i]/2)/5;
        }
 if(fabs(epsilon[0]) >= E_MAX || fabs(epsilon[4]) >= E_MAX)
        {
                for(i=0; i < N; i++)
                        state[i] = old[i];
                step = step/2;
                goto start;
        }
        time_ = time_ + step;
        if (fabs(epsilon[0]) <= E_MIN || fabs(epsilon[4]) <= E_MIN)
                step = step*2;
}

void rkm_()
{
        double k1[N],k2[N],k3[N],k4[N];
        int i;

        for(i=0;i<N;i++)
        {
        arg[i] = state[i];
        //printf("%d\t%lf\n",i,state[i]);
        }

deriv_();

        for(i=0;i<N;i++)
        {
        k1[i] = step *(deriv[i]);
        arg[i] = state[i]+k1[i]/2;
        }

deriv_();

        for(i=0;i<N;i++)
        {
        k2[i] = step * deriv[i];
        arg[i] = state[i]+ k2[i]/2;
        }


deriv_();

        for(i=0;i<N;i++)
        {
        k3[i] = step * deriv[i];
        arg[i] = state[i]+ k3[i];
        }

deriv_();

        for(i=0;i<N;i++)
        {
        k4[i] = step * deriv[i];
        state[i] = state[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
        }

        time_=time_+step;


}


cur()
{
FILE *psinf,*ps_tau,*ps, *vnull,*vaxonnull, *it, *ip, *itot,*zvbar,*ca,*var , *hinf , *nhrel; 
double v,p,i,q,m,k=0.43,a=0.5,s,t,u,iapd=0;
double n,h;
//vnull = fopen("vnull.xls","w");
//vaxonnull = fopen("vaxonnull.xls","w");
//
//ca = fopen("canull.xls","w"); 
//var = fopen("var.xls","w");
//zvbar = fopen("zvbar.xls","w");   
//
//psinf = fopen("psbar.xls","w");
//ps_tau = fopen("pstau.xls","w");
//
//hinf = fopen("hbar.xls","w");
//
//nhrel = fopen("nhrelation.xls","w");
	for(v=-80.0;v<=60;v=v+0.1)
	{
		
		//fprintf(psinf,"%lf\t%lf\n",v,  psbar(v) );
		//fprintf(ps_tau,"%lf\t%lf\n",v,  pstau(v) );
	
		m=kca*zv(v)*(vca-v) ;		
		//fprintf(ca,"%lf\t%lf\n",  (-1 + sqrt(1+8*m) )/4 ,v );
		//fprintf(ca,"%lf\t%lf\n",v,  m );
		//fprintf(zvbar,"%lf\t%lf\n",v,  zv(v) );
		//i = (ina(v,hbar(v))+ ik(v,nbar(v)) + ia(v,hai(v))+ il(v) );
		i = ( ia(v,hai(v))+ i_kf(v,bbar(v) ) + il(v) );
		q = gca * zv(v) * (v-vca);
		p = gkca * (v-vk);
		
		s = p+i;
		t = p*k + q + i*(a+k);
		u = a*(q+i*k);
		
		//fprintf(vnull,"%lf\t%lf\n",v,(-i*a-q)/(i+p));
		//fprintf(vnull,"%lf\t%lf\n",v,-(a*(i+q-w))/(-w+i+q+p));
		//fprintf(vnull, "%lf\t%lf\n",(-t + sqrt( pow(t,2) -4*s*u )) / (2*s) ,v      ); //(2*s) very important 5/2*6 is not the same as 5/(2*6)
		//fprintf(var,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",v,i,q,p,s,t,u);
	
		//fprintf(vaxonnull, "%lf\t%lf\n",v, (iext - gk* pow(nbar(v),4)*(v-vk) - il_a(v)) /(gna* pow(mbar(v),3)*(v-vna)));
		//fprintf(hinf,"%lf\t%lf\n",v,hbar(v));
		
	}

		for(v=-74.0;v<=60;v=v+0.001)
		{
			iapd = 0.227924;
			//fprintf(vaxonnull, "%lf\t%lf\n",v, (iext  - iapd - gk* pow(nbar(v),4)*(v-vk) - il_a(v)) /(gna* pow(mbar(v),3)*(v-vna)));//-0.274053
		}


	for(n=0;n<0.6;n=n+0.01)
	{
		h=-1.7253*n+0.9828; //n=(0.9822-h)/1.7253
		//fprintf(nhrel,"%lf\t%lf\n",n,h);
	}

}
