% Compute the priming metrics in the cells at the time when there is a
% coral bleaching report
% Author: Xinru Li; Date: Jun. 2022

clear

% load data
HS_ts=ncread('ts_HS_CoralBleachingDatabase.nc','ts_hs');
coor=ncread('ts_HS_CoralBleachingDatabase.nc','coor_cc');


%%
% priming metrics with the limit of at least 7 days which returns 2889 non-zero outputs     
nc=size(HS_ts,1);
Drec_c7=zeros(nc,1);         
Dtr_c7=zeros(nc,1);   
Atr_c7=zeros(nc,1);
Dc=zeros(nc,1);
Ac=zeros(nc,1);
HSpeak=zeros(nc,1);
hsmax_day=zeros(nc,1);
HRc=zeros(nc,1);

for n=1:nc    
  if (max(HS_ts(n,:))>0)
        hs_max=0;   % the SSTpeak within a HSY
    for r=1:(365-2)                
          hs=HS_ts(n,r);
          hs_1=HS_ts(n,r+1);
          hs_2=HS_ts(n,r+2);
          hs_3d=[hs, hs_1, hs_2]; 
          hs_3davg = mean(hs_3d);
          hs_max=max(hs_max,hs_3davg); 
    end
  
        HSpeak(n)=hs_max;
      
        dummy_hp=0;   % indicate that reach the the date of first SSTpeak
    for r=1:(365-2)
          hs=HS_ts(n,r);
          hs_1=HS_ts(n,r+1);
          hs_2=HS_ts(n,r+2);
          hs_3d=[hs, hs_1, hs_2]; 
          hs_3davg = mean(hs_3d);
      if (dummy_hp == 0)
        if (hs_3davg == HSpeak(n))
             hsmax_day(n)=r+1; % the date of first SSTpeak 
             dummy_hp=1;
        end
      end
    end 
    
 % find the date of the begining and end of the continuous HS period including HSpeak 
        dummy_c=0;  % a variable used to indicate the date of continuous sst>baseline
    for r=hsmax_day(n):-1:1  % if hsmax_day = 0, the loop will only run one relization
          hs=HS_ts(n,r);
      if (r >= 3)
           hs_1=HS_ts(n,r-1);
           hs_2=HS_ts(n,r-2);
        if (dummy_c == 0)
          if (hs == 0 && hs_1 == 0 && hs_2 == 0)
               hsbegin_day_c=r+1;   % ? the date of first continuous sst>baseline
               dummy_c=1;
          end
        end 
      else
        if (dummy_c == 0)
          if (hs == 0)
              hsbegin_day_c=r;
              dummy_c=1;
          end
        end
        if (dummy_c == 0)
          if (r==1)
               hsbegin_day_c=r;
               dummy_c=1;
          end
        end
      end
    end    
  
        HRc(n)=HSpeak(n)/(hsmax_day(n)-hsbegin_day_c+1);  
      
% calculate duration & accumulated HS during the continuous HS>0 period that includes the peak HS value
        dummy_e=0;
    for r=hsmax_day(n):365
          hs=HS_ts(n,r);
      if (r <= (365-2))
            hs_1=HS_ts(n,r+1);
            hs_2=HS_ts(n,r+2);
        if (dummy_e == 0)
          if (hs == 0 && hs_1 == 0 && hs_2 == 0)
               hslast_day_c=r-1;  % potential errors in singular 0 in a period of postive HS e.g. lon:224, lat=90, y=22
               dummy_e=1;
          end
        end
      else
        if (dummy_e == 0)
          if (hs == 0)
              hslast_day_c=r;
              dummy_e=1;
          end
        end
        if (dummy_e == 0)
          if (r == 365)
               hslast_day_c=r;
               dummy_e=1;
          end
        end   
      end
    end  
  
        Dc(n)=hslast_day_c-hsbegin_day_c+1; 
        Ac(n)=sum(HS_ts(n,hsbegin_day_c:hslast_day_c));
        
% detect priming        
    if (hsmax_day(n) < 16)
         Drec_c7(n)=0;         
         Dtr_c7(n)=0;   
         Atr_c7(n)=0; 
    else
 % find end date of recovery period
          HS_rec_e=0;
      for r=hsmax_day(n):-1:14  
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
               Drec_c7(n)=HS_rec_e-HS_rec_s_v2+1;
          end
        end
        if (dummy_rec_s==1)
             break
        end
        if (r==7)
             Drec_c7(n)=0;
        end
      end  
  
% find the start of pulsed training period (i.e.,V2 approach)
      if (Drec_c7(n) == 0)
           Dtr_c7(n)=0;
           Atr_c7(n)=0;
      else
           dummy_trc_s=0;
           HS_trc_e=HS_rec_s_v2-1;
        for r=HS_trc_e:-1:3
              hs=HS_ts(n,r);
              hs_1=HS_ts(n,r-1);
              hs_2=HS_ts(n,r-2);           
          if (hs>=0.5 && hs_1>=0.5 && hs_2>=0.5)
                dummy_trc_s=1;
                Dtr_c7(n)=0;
                Atr_c7(n)=0;  
          elseif (hs==0 && hs_1==0 && hs_2==0)
                HS_trc_s=r+1;
                Dtrc=HS_trc_e-HS_trc_s+1;
            if (Dtrc >= 7)              % control of Dtr length to be at least 3 days
                 dummy_trc_s=1;
                 Dtr_c7(n)=Dtrc;
                 Atr_c7(n)=sum(HS_ts(n,HS_trc_s:HS_trc_e));
            else
                 dummy_trc_s=1;
                 Dtr_c7(n)=0;
                 Atr_c7(n)=0;
            end              
          end        
          if (dummy_trc_s==1)
               break
          end  
          if (r==3)
               Dtr_c7(n)=HS_trc_e;
               Atr_c7(n)=sum(HS_ts(n,1:HS_trc_e));
               dummy_trc_s=1;
          end          
        end
      end 
    end
  else
        Drec_c7(n)=0;         
        Dtr_c7(n)=0;   
        Atr_c7(n)=0;
        Dc(n)=0;
        Ac(n)=0;
        HSpeak(n)=0;
        hsmax_day(n)=0;
        HRc(n)=0;
  end
end


%% write out to the spreadsheet
metrics=[Dre_hsbegin_c7,Dtr_c7_v2,Atr_c7_v2];
T=array2table(metrics,"VariableNames",["D_rec7","D_tr7","A_tr7"]);
writetable(T,'HSpeak1_pulsedPriming_Mmetrics_c14_test.xlsx');

    