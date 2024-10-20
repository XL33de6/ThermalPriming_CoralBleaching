% Compute the priming metrics in all the coral reefs cells globally
% followed the Millennium Coral Reef Mapping Project data (UNEP, 2021)
% Author: Xinru Li; Date: Aug. 2022

clear


% load data
DHD = ncread('DHD_MMMct5km_cc.nc','DHD_hsy');
DHD_coor = ncread('DHD_MMMct5km_cc.nc','coor_cc');
HS_coor = ncread('HS_1994_MMMct5km_cc_nmc.nc','coor_cc');

year = linspace(1986,2020,35);
yr_str = string(year);
HSY = linspace(1986,2019,34);

day_m_nlpy = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % the daily data in non leap year leave out the day 28 in Feb.
day_m_lpy = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

days_nlpy=365;
days_lpy=366;

x_sel=ncread('LMC_mk_sim_ct5km_v3.1_occ_nnan.nc','x_sel_ocean');


%%
tic

nc=size(DHD,1);              
Dtr_c=zeros(nc,34);
Atr_c=zeros(nc,34);
Dtr=zeros(nc,34);
Atr=zeros(nc,34);      
Dre=zeros(nc,34);
Asig=zeros(nc,34);  
Dsig=zeros(nc,34);
n_sig=zeros(34,1);
n_pr=zeros(34,1); 

for y=1:34
      filepath_p1 = 'F:/RC3/coraltemp_MMM_Hotspots for coral grids/HS_';   
      filepath_p2_1 = yr_str(y);
      filepath_p2_2 = yr_str(y+1);
      filepath_p3 = '_MMMct5km_cc_nmc.nc';       % change the value depending on which subset of dataset is being used
      filepath_cat1 = strcat(filepath_p1,filepath_p2_1,filepath_p3);
      filepath1 = char(filepath_cat1);  
      HS_cc_1 = ncread(filepath1,'hs');  % size(sst)=(7200,2400,#days:1), scale down to true value
      filepath_cat2 = strcat(filepath_p1,filepath_p2_2,filepath_p3);
      filepath2 = char(filepath_cat2);  
      HS_cc_2 = ncread(filepath2,'hs');  % size(sst)=(7200,2400,#days:1), scale down to true value  
      coor = ncread(filepath2,'coor_cc');
      e1_b=122;
      s2_b=123;
      e1_a=306;
      s2_a=307;
  if (mod(year(y),4) ~= 0)
    if (mod(year(y+1),4) ~= 0)   % determine the #days of heat stress year, use 'HSY_total' since Feb of a HSY is in the next calendar year
         nd=days_nlpy;
         m_day=day_m_nlpy;    % the #days of each month of a non leap year using the HSY definition b
         HS_hsy=zeros(nd,1);
         s1_b=244;      % definition b from Sept. - Aug.
         e2_b=365;
         e3_b=365;
         e4_b=243;
         s1_a=60;       % definition a from Mar. - Feb.
         e2_a=365;
         e3_a=365;
         e4_a=59;
    else
         nd=days_lpy;
         m_day=day_m_lpy;    % the #days of each month of a leap year using the HSY definition b
         HS_hsy=zeros(nd,1); 
         s1_b=244;
         e2_b=365;
         e3_b=366;
         e4_b=244;
         s1_a=60;
         e2_a=365;
         e3_a=366;
         e4_a=60;
    end
  else
       nd=days_nlpy;
       m_day=day_m_nlpy;    % the #days of each month of a non leap year using the HSY definition b
       HS_hsy=zeros(nd,1);
       s1_b=245;  % definition b from Sept. -Aug.
       e2_b=366;
       e3_b=365;
       e4_b=243;
       s1_a=61;
       e2_a=366;
       e3_a=365;
       e4_a=59;
  end

      n_s=0;
      n_prime=0;
  for n=1:nc
    if (x_sel(n) == 4)  % HSY definition b Sept.-Aug.
         HS_cc(1:e1_b) = squeeze(HS_cc_1(n,s1_b:e2_b));
         HS_cc(s2_b:e3_b) = squeeze(HS_cc_2(n,1:e4_b));
    elseif (x_sel(n) == 10) % HSY definition a: Mar. - Feb.
         HS_cc(1:e1_a) = squeeze(HS_cc_1(n,s1_a:e2_a));
         HS_cc(s2_a:e3_a) = squeeze(HS_cc_2(n,1:e4_a));
    end
  
    if (DHD(n,y)>=56)
          n_s=n_s+1;
          hs_max=0;   % the SSTpeak within a HSY
          dummy=0;    % a variable used to indicate the date of first sst>baseline
      for r=1:(nd-12)                
            hs=HS_cc(r);
            hs_1=HS_cc(r+1);
            hs_2=HS_cc(r+2);
            hs_3d=[hs, hs_1, hs_2]; 
            hs_3davg = mean(hs_3d);
            hs_max=max(hs_max,hs_3davg); 
        if (dummy == 0)
          if (hs > 0 && hs_1 > 0 && hs_2 > 0)
            if (r == 1)
                 hsbegin_day=r;    % if the first day of a HSY has HS>0
                 dummy=1;
            else
                 hsbegin_day=r-1;  % the day before the first day SST>baseline  
                 dummy=1;
            end
          end 
          if (r == (nd-12))
               hsbegin_day=0;
          end 
        end            
      end 
      
      % find the date of first SSTpeak within a HSY
          dummy_hp=0;   % indicate that reach the the date of first SSTpeak
      for r=1:(nd-2)
            hs=HS_cc(r);
            hs_1=HS_cc(r+1);
            hs_2=HS_cc(r+2);
            hs_3d=[hs, hs_1, hs_2]; 
            hs_3davg = mean(hs_3d);
        if (dummy_hp == 0)
          if (hs_3davg == hs_max)
               hsmax_day=r+1; % the date of first SSTpeak 
               dummy_hp=1;
          end
        end
      end 
    
      % find Asig
          dummy_e=0;
      for r=hsmax_day:nd
            hs=HS_cc(r);
        if (r <= (nd-2))
              hs_1=HS_cc(r+1);
              hs_2=HS_cc(r+2);            
          if (dummy_e == 0)
            if (hs == 0 && hs_1 == 0 && hs_2 == 0)
                  hslast_day_c=r;  % potential errors in singular 0 in a period of postive HS e.g. lon:224, lat=90, y=22
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
            if (r == nd)
                 hslast_day_c=r;
            end
          end   
        end         
      end 
      
      % find the date of first continuous positive HSp90 
          dummy_c=0;  % a variable used to indicate the date of continuous sst>baseline
      for r=hsmax_day:-1:9  % if hsmax_day = 0, the loop will only run one relization
            hs=HS_cc(r);
            hs_1=HS_cc(r-1);
            hs_2=HS_cc(r-2);
        if (dummy_c == 0)
          if (hs == 0 && hs_1 == 0 && hs_2 == 0)
                hsbegin_day_c=r;   % ? the date of first continuous sst>baseline
                dummy_c=1;
          end
          if (r == 9)
            if (hs ~= 0 | hs_1 ~=0 | hs_2 ~=0)
                 hsbegin_day_c=0;
            end
          end
        end
      end          
      
      % find the start of recovery period
      if (hsbegin_day_c == 0)
           HS_rec_s=0;
      else
            dummy_rec_s=0;
        for r=hsbegin_day_c:-1:5
              hs=HS_cc(r);
              hs_1=HS_cc(r-1);
              hs_2=HS_cc(r-2);
          if (dummy_rec_s == 0)
            if (hs>0 && hs_1>0 && hs_2>0)
                 HS_rec_s=r-1;
                 dummy_rec_s=1;
            end
            if (r == 5)
              if (hs == 0 | hs_1 ==0 | hs_2 ==0)
                   HS_rec_s=0;
              end
            end
          end
        end
      end
      
      % find the start of priming period
      if (HS_rec_s == 0)
           HS_tr_s=0;
      else
            dummy_tr_s=0;
        for r=HS_rec_s:-1:3                
              hs=HS_cc(r);
              hs_1=HS_cc(r-1);
              hs_2=HS_cc(r-2);    
          if (dummy_tr_s == 0)
            if (hs==0 && hs_1==0 && hs_2==0)
                 HS_tr_s=r;
                 dummy_tr_s=1;
                 n_prime=n_prime+1;
            end
            if (r == 3)
              if (hs ~= 0 | hs_1 ~=0 | hs_2 ~=0)
                   HS_tr_s=1;
                   n_prime=n_prime+1;
              end
            end
          end         
        end         
      end    
    
      % find continuous training period
      if (HS_tr_s~=0 && HS_rec_s~=0)
           Dtr_c(n,y)=HS_rec_s-HS_tr_s+1;
           Atr_c(n,y)=sum(HS_cc(HS_tr_s:HS_rec_s));
      else
           Dtr_c(n,y)=0;
           Atr_c(n,y)=0;
      end
      % find incotinuous training period
      if (hsbegin_day~=0 && HS_rec_s~=0)
           Dtr(n,y)=HS_rec_s-hsbegin_day+1;
           Atr(n,y)=sum(HS_cc(hsbegin_day:HS_rec_s));
      else
           Dtr(n,y)=0;
           Atr(n,y)=0;
      end    
      % find recovery period
      if (hsbegin_day_c~=0 && HS_rec_s~=0)
           Dre(n,y)=hsbegin_day_c-HS_rec_s+1;
      else
           Dre(n,y)=0;
      end  
      % find Asig
      if (hsbegin_day_c~=0 && hslast_day_c~=0)
           Asig(n,y)=sum(HS_cc(hsbegin_day_c:hslast_day_c));
           Dsig(n,y)=hslast_day_c-hsbegin_day_c+1;
      else
           Asig(n,y)=0;
           Dsig(n,y)=0;
      end   
    end       
  end  
      n_sig(y)=n_s;
      n_pr(y)=n_prime;      
end

toc



%%
% write out metrics involved in the priming phenomenon
 ncnc= netcdf.create('priming_metrics_c749_HSY8619_DHD56.nc','NC_WRITE');   

nID=netcdf.defDim(ncnc,'the number of coral cells',nc);
coorID=netcdf.defDim(ncnc,'two columns for coordinate',2);

tID=netcdf.defDim(ncnc,'time',34);
vtID=netcdf.defVar(ncnc,'time','float',tID);
netcdf.putAtt(ncnc,vtID,'axis','T');
long_namet = 'the heat stress year corresponding to each DHD';
netcdf.putAtt(ncnc,vtID,'long_name',long_namet);

varname = 'Dtrc';
long_name = 'the number of days druing the training period with continuous HS';
unit = 'day';
vm1ID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vm1ID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vm1ID,'units',unit);          % The unit

varname = 'Atrc';
long_name = 'the accumulated heat stress magnitude during the Dtr_c period';
unit = 'degree celcius*day';
vm3ID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vm3ID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vm3ID,'units',unit);          % The unit

varname = 'Ac';
long_name = 'the accumulated heat stress magnitude for the continuous period with HSpeak == Ac';
unit = 'degree celcius*day';
vm5ID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vm5ID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vm5ID,'units',unit);          % The unit

varname = 'Dc';
long_name = 'the number of days during the continuous duration with HSpeak == Dc';
unit = 'day';
vm6ID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vm6ID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vm6ID,'units',unit);          % The unit

varname = 'Drec';
long_name = 'the number of days druing the recovery period';
unit = 'day';
vm7ID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vm7ID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vm7ID,'units',unit);          % The unit
      
var = 'coor_cc';
long_name = 'coordinate of coral cells without NaNs';
unit = 'degree celcius';
vcoorID=netcdf.defVar(ncnc,var,'float',[nID,coorID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vcoorID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vcoorID,'units',unit);          % The unit         

% end define mode
netcdf.endDef(ncnc)
% input data
netcdf.putVar(ncnc,vm1ID,Dp);
netcdf.putVar(ncnc,vm3ID,Ap);
netcdf.putVar(ncnc,vm5ID,Ac);
netcdf.putVar(ncnc,vm6ID,Dc);
netcdf.putVar(ncnc,vm7ID,Dr);
netcdf.putVar(ncnc,vcoorID,coor);
netcdf.putVar(ncnc,vtID,HSY);
netcdf.close(ncnc)    


