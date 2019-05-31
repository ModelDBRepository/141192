#include <stdio.h>
#include <math.h>
#include "ABPDcon.h"
#include "ABPD.h"
#include <stdlib.h>

main(int argc, char *argv[])
{

//-------------------------------
	/* variables */
//-------------------------------
//
int bt_count=0,sp_on=1,i,j=0;
double in_on=3650, in_off, stim_dur,increment=0.00,gsyn_strength,phase;//0.0025
FILE *tm_voltd,*tm_volts,*tm_volts_on,*tm_volts_off,*v_ca, *v_ca_on,*v_ca_off,*tm_voltpd,*tm_volta, *tm_c,*tm_ha,*tm_z,*tm_c2, *z_c,*ha_z,*inact_ca_act_kca, *prc1,*prc2,*tm_zc2, *soma_cur, *dend_cur,*sdr_cur,*t_dend_cur,*t_soma_cur, *t_soma_sum_cur, *t_dend_sum_cur, *t_soma_cur0, *t_dend_cur0,*tm_stim,*t_dend_cur_on,*t_soma_cur_on, *tm_Idpd;

FILE *va_h,*va_h_on,*va_h_off,*vpd_ps,*vpd_ps_off,*vpd_ps_on;
FILE *t_ids_cur0,*ds_sum,*t_ias_cur0,*as_sum,*t_ids,*ids_on,*t_ias,*ias_on,*volt_d_a;
double sp_time_old, sp_time, ISI_old, ISI, bt_on_time_old=0.0,  bt_on_time=0.0, period0, period1=0, period2, period1_st=0 ;
double isoma_sum=0,idend_sum=0,ids_sum=0,isoma0[5000],idend0[5000],dummy[5000],ids0[5000],ias_sum=0;

//scanf("%lf",&increment);

tm_volts = fopen("tm_volts.xls","w");
//tm_volts_off = fopen("tm_volts_off.xls","w");
//tm_volts_on = fopen("tm_volts_on.xls","w");

//tm_voltd = fopen("tm_voltd.xls","w");tm_volta = fopen("tm_volta.xls","w");tm_voltpd = fopen("tm_voltpd.xls","w");

//prc1=fopen("prc1.xls","a");
//prc2=fopen("prc2.xls","a");
//
//v_ca = fopen("v_ca.xls","w");
//v_ca_off = fopen("v_ca_off.xls","w");
//v_ca_on = fopen("v_ca_on.xls","w");
//va_h = fopen("va_h.xls","w");
//va_h_off = fopen("va_h_off.xls","w");
//va_h_on = fopen("va_h_on.xls","w");
//
//tm_stim = fopen("tm_stim.xls","w");
//t_ica=fopen("tm_ica.xls","w");
//t_ikca=fopen("tm_ikca.xls","w");
//tm_c = fopen("tm_c.xls","w");
//rate_ca= fopen("rate_ca.xls","w");
//volt_d_a= fopen("volt_d_a.xls","w");
//tm_Idpd = fopen("tm_Idpd.xls","w");
//
//vpd_ps = fopen("vpd_ps.xls","w");
//vpd_ps_off = fopen("vpd_ps_off.xls","w");
//vpd_ps_on = fopen("vpd_ps_on.xls","w");



//------------------------------------
	/*setting initial state*/
//------------------------------------

time_ = START_TIME;
step = STEP;
 
state[0] = -20;//v
state[6] = -20;
state[7] = -20;

state[8] = psbar(state[0]);
state[2] = an(state[0])/(an(state[0])+bn(state[0]));//n 
state[3] = ah(state[0])/(ah(state[0])+bh(state[0]));//h
state[4] = zv(state[0]);//z
state[5] = hai(state[0]);//ha
state[1] = 0.6 ;//c


	
//	gsyn_strength = strtod(argv[1],NULL); //0.1
//	phase = strtod(argv[2],NULL); //0.5
//	stim_dur = strtod(argv[3],NULL); //50000 
	gsyn_strength = 0; //0.1
	phase = 0; //0.5
	stim_dur = 0; //50000 

 
/*
for(i=0;i<5000;i++)
	{
	fscanf(t_soma_cur0, "%lf",&dummy[i]);	
	fscanf(t_soma_cur0, "%lf",&isoma0[i]);	
	fscanf(t_dend_cur0, "%lf",&dummy[i]);	
	fscanf(t_dend_cur0, "%lf",&idend0[i]);	
	//printf("%lf\t%lf\n",dummy[i],isoma0[i]);
	fscanf(t_ids_cur0, "%lf",&dummy[i]);	
	fscanf(t_ids_cur0, "%lf",&ids0[i]);
	}
*/
while(time_ <= END_TIME)
{
		fprintf(tm_volts, "%lf\t%lf\n",time_,state[0]);
		if(state[0]> -38 && sp_on==0 )
		{
			sp_time_old=sp_time;
			sp_time = time_;				
			ISI_old = ISI;
			ISI = sp_time - sp_time_old;
			
			if(ISI > 150)
			{
				bt_count = bt_count + 1;
				bt_on_time_old = bt_on_time;
				bt_on_time = time_;
				//printf("%lf\t%lf\n",time_,bt_on_time-bt_on_time_old);	//--------------
			

				if(bt_count == 4)
				{
					period1_st = time_;
					period0 = bt_on_time-bt_on_time_old;		
				}
				if(bt_count == 5)
				{	
					period1 = bt_on_time-bt_on_time_old;	
					//printf("%lf\t%lf\n",phase,(period1-period0)/period0);		
					
				}

				 if(bt_count == 6)
				{	
					period2 = bt_on_time-bt_on_time_old;	
					//printf("%lf\t%lf\n",phase,(period2-period0)/period0);
				} 		
				sp_on =1;
			}
		}

	if(sp_on == 1 && state[0]<-45)
		sp_on=0;


//printf( "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",time_,state[0],state[1],state[2],state[3],state[4],state[5]);
        if ((time_>1000)&&(time_<END_TIME))
        {
/*
		if(time_> 7711.360000 && time_< 7711.360000+972.72)
		{
			fprintf(tm_Idpd,"%f\t%f\t\n",-(7711.360000-time_)/972.720,Gdpd*(state[6]-state[9]) );
			fprintf(t_ica,"%f\t%f\t\n",-(7711.360000-time_)/972.720,tm_ica );
			fprintf(t_ikca,"%f\t%f\t\n",-(7711.360000-time_)/972.720,tm_ikca);
		}
*/
//			fprintf(tm_Idpd,"%f\t%f\t\n",time_,Gdpd*(state[6]-state[9]) );
//			fprintf(t_ica,"%f\t%f\t\n",time_,tm_ica );
//			fprintf(t_ikca,"%f\t%f\t\n",time_,tm_ikca);		
//		
//		fprintf(volt_d_a, "%lf\t%lf\n",state[7],state[6]);
//		fprintf(tm_voltpd, "%lf\t%lf\n",time_,state[9]);
//		fprintf(tm_voltd, "%lf\t%lf\n",time_,state[6]);
//		fprintf(tm_c, "%lf\t%lf\n",time_,state[1]); //c/(0.5+c)
//		fprintf(tm_volta, "%lf\t%lf\n",time_,state[7]);
		//fprintf(tm_stim, "%f\t%f\t\n",time_,-60.0 + (20*gsyn/(gsyn+0.01) ) );	//-------------------------------
		//fprintf(t_ica,"%lf\t%lf\n",time_,tm_ica);
		//fprintf(t_ikca,"%lf\t%lf\n",time_,tm_ikca);
		
		/*
		fprintf(t_soma_cur, "%lf\t%lf\n",time_,isoma);	
		fprintf(t_dend_cur, "%lf\t%lf\n",time_,idend);	
		fprintf(t_ids,"%lf\t%lf\t",time_,i_ds);
		fprintf(t_ias,"%lf\t%lf\t",time_,i_as);
		
		

		
		
		
		fprintf(tm_c2, "%lf\t%lf\n",time_,1/(ca_inact_var+state[1]));
		fprintf(tm_z, "%lf\t%lf\n",time_,state[4]);
		fprintf(tm_ha, "%lf\t%lf\n",time_,state[5]);
		fprintf(tm_zc2, "%lf\t%lf\n",time_,state[4] * ( 1/(ca_inact_var+state[1]) )    );
		fprintf(inact_ca_act_kca, "%lf\t%lf\n",1/(ca_inact_var+state[1]),state[1]/(0.5+state[1]) );
		
		fprintf(ha_z, "%lf\t%lf\n",state[5],state[4]);
		fprintf(z_c, "%lf\t%lf\n",state[4],state[1]);
		
		
               fprintf(tm_volt, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",time_,state[0],state[1],state[2],state[3],state[4],state[5]); 

		*/
	}

		if(period1_st> 0.0 && time_ > period1_st + phase * (period0)  && time_ < period1_st + phase * (period0) + stim_dur)
		{
			gsyn= gsyn_strength;
		}
	else
		gsyn=0;  
	
	if (time_>1000 && period1_st> 0.0 && time_ >= period1_st + phase * (period0) + stim_dur)
        	{ 
//          		fprintf(va_h_off,"%f\t%f\t\n",state[7],state[3]);
//			fprintf(vpd_ps_off,"%f\t%f\t\n",state[8],state[9]); 
//			fprintf(v_ca_off,"%f\t%f\t\n",state[1],state[6]); //----------------------------------
			//fprintf(tm_volts_off, "%lf\t%lf\n",time_,state[0]);
		}        
	else if(time_ >1000 && period1_st> 0.0 && time_ > period1_st + phase * (period0)  && time_ < period1_st + phase * (period0) + stim_dur)
		{
//			fprintf(va_h_on,"%f\t%f\t\n",state[7],state[3]);
//			fprintf(vpd_ps_on,"%f\t%f\t\n",state[8],state[9]);
//			fprintf(v_ca_on,"%f\t%f\t\n",state[1],state[6]); //-----------------------
			//fprintf(tm_volts_on, "%lf\t%lf\n",time_,state[0]);
		/* 
		isoma=isoma+(-isoma0[j]);
		idend=idend+(-idend0[j]);
		if(isoma>0)
			isoma_sum=isoma+isoma_sum;
		if(isoma<0)
			isoma_sum=-isoma+isoma_sum;

		
		if(idend>0)
			idend_sum=idend+idend_sum;
		if(idend<0)
			idend_sum=-idend+idend_sum;

		fprintf(t_soma_cur, "%lf\t%lf\n",time_,isoma);		
		fprintf(t_dend_cur, "%lf\t%lf\n",time_,idend);
		fprintf(t_soma_sum_cur, "%lf\t%lf\n",time_,isoma_sum);		
		fprintf(t_dend_sum_cur, "%lf\t%lf\n",time_,idend_sum);
		//printf("%lf\t%lf\t%lf\n",time_,isoma,isoma0[j]);	
		
		

		//i_ds=i_ds+(-ids0[j]);
		ids_sum=i_ds+ids_sum;
		
		j++;
		
		fprintf(t_dend_cur_on, "%lf\t%lf\n",time_,idend);
		fprintf(t_soma_cur_on, "%lf\t%lf\n",time_,isoma);
		fprintf(ids_on,"%lf\t%lf\n",time_,i_ds);
		fprintf(ias_on,"%lf\t%lf\n",time_,i_as);*/
		}
	else
		{
			if (time_>1000)
			{
//				fprintf(va_h,"%f\t%f\t\n",state[7],state[3]);
//				fprintf(v_ca,"%f\t%f\t\n",state[1],state[6]); //---------------------------------
				//fprintf(tm_volts, "%lf\t%lf\n",time_,state[0]);
//				fprintf(vpd_ps,"%f\t%f\t\n",state[8],state[9]);
			}	 
		}

	rkm_();

		
}


//cur();
fclose(tm_volts);
/*
fprintf(ds_sum, "%lf\t%lf\n",gsyn_strength,ids_sum);
fprintf(as_sum, "%lf\t%lf\n",gsyn_strength,ias_sum);
fprintf(soma_cur, "%lf\t%lf\n",gsyn_strength,isoma_sum);		
fprintf(dend_cur, "%lf\t%lf\n",gsyn_strength,idend_sum);
if(gsyn_strength>0)
fprintf(sdr_cur, "%lf\t%lf\n",gsyn_strength,idend_sum/isoma_sum);
*/
}
