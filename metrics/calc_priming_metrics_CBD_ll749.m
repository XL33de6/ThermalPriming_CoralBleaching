% Compute the priming metrics in the cells at the time when there is a
% coral bleaching report
% Author: Xinru Li; Date: Jun. 2022

clear


% load data
HS_ts=ncread('ts_HS_CoralBleachingDatabase.nc','ts_hs');
coor=ncread('ts_HS_CoralBleachingDatabase.nc','coor_cc');
date_DHDmax_ph41=ncread('DHDmax_MMMct5km_CBD_ph42.nc','Date_ph41_DHDmax');
date_DHDmax=date_DHDmax_ph41+41;


%%
tic

% priming metrics with the limit of at least 3 days which returns 11294 non zero outputs
nc=size(HS_ts,1);
Dre_hsbegin_c3=zeros(nc,1);         
Dtr_c3_v2=zeros(nc,1);   
Atr_c3_v2=zeros(nc,1);
Dc_c3=zeros(nc,1);
HSpeak_c3=zeros(nc,1);
HRc_c3=zeros(nc,1);


for n=1:nc    
 % find end date of recovery period
      HS_rec_e=0;
  for r=date_DHDmax(n):-1:7  
    if (HS_ts(n,r)==0 && HS_ts(n,r-1)==0 && HS_ts(n,r-2)==0)  % the end of recovery period with control of at least three continuous days of 0 HS
         HS_rec_e=r;  
    end
    if (HS_rec_e~=0)
        break
    end
  end
      dummy_rec_s=0;
  for x=HS_rec_e:-1:7
    if (dummy_rec_s==0)
         hs=HS_ts(n,x);
         hs_1=HS_ts(n,x-1);
         hs_2=HS_ts(n,x-2);
      if (hs>0 && hs_1>0 && hs_2>0)  % 3 days error control
           dummy_rec_s=1;
           HS_rec_s_v2=x+1;
           Dre_hsbegin_c3(n)=HS_rec_e-HS_rec_s_v2+1;
      end
    end
    if (dummy_rec_s==1)
         break
    end
    if (r==7)
         Dre_hsbegin_c3(n)=0;
    end
  end 
  
% find the start of pulsed training period (i.e.,V2 approach)
  if (Dre_hsbegin_c3(n) == 0)
       Dtr_c3_v2(n)=0;
       Atr_c3_v2(n)=0;
  else
        dummy_trc_s=0;
        HS_trc_e=HS_rec_s_v2-1;
    for r=HS_trc_e:-1:3
      if (dummy_trc_s==0)
           hs=HS_ts(n,r);
           hs_1=HS_ts(n,r-1);
           hs_2=HS_ts(n,r-2);           
        if (hs==0 && hs_1==0 && hs_2==0)
             HS_trc_s=r+1;
             Dtrc=HS_trc_e-HS_trc_s+1;    
          if (Dtrc >= 3)              % control of Dtr length to be at least 3 days
               dummy_trc_s=1;
               Dtr_c3_v2(n)=Dtrc;
               Atr_c3_v2(n)=sum(HS_ts(n,HS_trc_s:HS_trc_e));
          else
               dummy_trc_s=1;
               Dtr_c3_v2(n)=0;
               Atr_c3_v2(n)=0;
          end      
        end 
        if (dummy_trc_s==1)
             break
        end        
        if (r==3)
             Dtr_c3_v2(n)=0;
             Atr_c3_v2(n)=0;
        end 
      end
    end
  end
  
end

toc


%%
% priming metrics with the limit of at least 7 days which return 5835 non zero outputs 
nc=size(HS_ts,1);
Dre_hsbegin_c7=zeros(nc,1);         
Dtr_c7_v2=zeros(nc,1);  
Atr_c7_v2=zeros(nc,1);

for n=1:nc    
 % find end date of recovery period
      HS_rec_e=0;
  for r=data_DHDmax(n):-1:15  
    if (HS_ts(n,r)==0 && HS_ts(n,r-1)==0 && HS_ts(n,r-2)==0 && HS_ts(n,r-3)==0 && HS_ts(n,r-4)==0 && HS_ts(n,r-5)==0 && HS_ts(n,r-6)==0)  % the end of recovery period with control of at least three continuous days of 0 HS
         HS_rec_e=r;  
    end
    if (HS_rec_e~=0)
        break
    end
  end
      dummy_rec_s=0;
  for x=HS_rec_e:-1:15
    if (dummy_rec_s==0)
         hs=HS_ts(n,x);
         hs_1=HS_ts(n,x-1);
         hs_2=HS_ts(n,x-2);
      if (hs>0 && hs_1>0 && hs_2>0)  % 3 days error control
           dummy_rec_s=1;
           HS_rec_s_v2=x+1;
           Dre_hsbegin_c7(n)=HS_rec_e-HS_rec_s_v2+1;
      end
    end
    if (dummy_rec_s==1)
         break
    end
    if (r==15)
         Dre_hsbegin_c7(n)=0;
    end
  end 
  
% find the start of pulsed training period (i.e.,V2 approach)
  if (Dre_hsbegin_c7(n) == 0)
       Dtr_c7_v2(n)=0;
       Atr_c7_v2(n)=0;
  else
        dummy_trc_s=0;
        HS_trc_e=HS_rec_s_v2-1;
    for r=HS_trc_e:-1:7
      if (dummy_trc_s==0)
           hs=HS_ts(n,r);
           hs_1=HS_ts(n,r-1);
           hs_2=HS_ts(n,r-2);           
        if (hs==0 && hs_1==0 && hs_2==0)
             HS_trc_s=r+1;
             Dtrc=HS_trc_e-HS_trc_s+1;    
          if (Dtrc >= 7)              % control of Dtr length to be at least 3 days
               dummy_trc_s=1;
               Dtr_c7_v2(n)=Dtrc;
               Atr_c7_v2(n)=sum(HS_ts(n,HS_trc_s:HS_trc_e));
          else
               dummy_trc_s=1;
               Dtr_c7_v2(n)=0;
               Atr_c7_v2(n)=0;
          end      
        end 
        if (dummy_trc_s==1)
             break
        end        
        if (r==7)
             Dtr_c7_v2(n)=0;
             Atr_c7_v2(n)=0;
        end 
      end
    end
  end
  
end
           


%% write out to the spreadsheet
metrics=[Dre_hsbegin_c7,Dtr_c7_v2,Atr_c7_v2];
T=array2table(metrics,"VariableNames",["D_rec7","D_tr7","A_tr7"]);
writetable(T,'ph42_DHDmax_pulsedPriming_results.xlsx');

    