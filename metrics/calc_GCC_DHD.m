% filter out the cells expericend DHD>56 and the years the significant heat
% stress events occurred
% Author: Xinru Li; Date: Feb. 2022

clear


% extract values in coral cells
coral_mask = ncread('F:/RC3/coraltemp_MMM_Hotspots for coral grids/WCMC_reef_cells_1_Resample_0.05.nc','WCMC_reef_cells_1_Resample',[1, 901],[7200,1800]);
ind=[3601:7200 1:3600]; % Move left side to right, and then lon: 179.975W-0-179.975E, consistent with MMM
coral_mask_cp=coral_mask(ind,:); % W to E

coral_mask_cp(coral_mask_cp~=1)=NaN;
[row, col] = find(~isnan(coral_mask_cp)); % check!!!

%
nc=length(row);
coor = cast(zeros(nc,2),'single'); 
coor(:,1)=Lat(col);
coor(:,2)=Long(row);
x_sel = x(~isnan(coral_mask_cp));

for n=1:nc
  if (isnan(x_sel(n)))
       x_sel(n)=x_sel(n-1);
  end
end

%%
%-------------------------write out simplified LMCM------------------------------------------
ncnc = netcdf.create('LMC_mk_sim_ct5km_v3.1_cc_nnan.nc','NC_WRITE');

nID=netcdf.defDim(ncnc,'the number of coral cells',nc);
coorID=netcdf.defDim(ncnc,'two columns for coordinate',2);
var = 'coor_cc';
long_name = 'coordinate of coral cells';
unit = 'degree celcius';
vcoorID=netcdf.defVar(ncnc,var,'float',[nID,coorID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vcoorID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vcoorID,'units',unit);          % The unit 

varname = 'x';
long_name = 'the number of months in a first calender year consisting of HSY';  
unit = '1';
vmmmID=netcdf.defVar(ncnc,varname,'float',nID); % we need to define axis of the field
netcdf.putAtt(ncnc,vmmmID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vmmmID,'units',unit);          % The unit        

% end define mode
netcdf.endDef(ncnc)
% input data
netcdf.putVar(ncnc,vmmmID,x_sel);
netcdf.putVar(ncnc,vcoorID,coor);
netcdf.close(ncnc)


%%
tic

% calculate DHD for all coral cells
year = (1985:2020);
HSY =(1986:2019);  % heat stress years
yr_str = string(year);
t=length(HSY);

day_m_nlpy = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % the daily data in non leap year leave out the day 28 in Feb.
day_m_lpy = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
days_nlpy=448;
days_lpy=449;

cc = ncread('HS_2020_MMMct5km_cc_REV.nc','coor_cc');
nc = size(cc,1);
% n_row_nan_total = zeros(10,34);
dhd = zeros(10,34);

for y=2%:35   
      filepath_p1 = 'F:/RC3/coraltemp_MMM_Hotspots for coral grids/HS_';
      filepath_p2_1 = yr_str(y-1);
      filepath_p2_2 = yr_str(y);
      filepath_p2_3 = yr_str(y+1);
      filepath_p3 = '_MMMct5km_cc_REV.nc';       % change the value depending on which subset of dataset is being used
      filepath_cat1 = strcat(filepath_p1,filepath_p2_1,filepath_p3);
      filepath1 = char(filepath_cat1);  
      HS_hsy_1 = ncread(filepath1,'hs');  % size(sst)=(7200,2400,#days:1), scale down to true value
      filepath_cat2 = strcat(filepath_p1,filepath_p2_2,filepath_p3);
      filepath2 = char(filepath_cat2);  
      HS_hsy_2 = ncread(filepath2,'hs');  % size(sst)=(7200,2400,#days:1), scale down to true value  
      filepath_cat3 = strcat(filepath_p1,filepath_p2_3,filepath_p3);
      filepath3 = char(filepath_cat3);  
      HS_hsy_3 = ncread(filepath3,'hs');
      
  if (mod(year(y-1),4) ~= 0) 
    if (mod(year(y),4) ~= 0)
      if (mod(year(y+1),4) ~= 0)   % determine the #days of heat stress year, use 'HSY_total' since Feb of a HSY is in the next calendar year
           days=days_nlpy;      % condition (1): 1985, 1986, 1987
%            m_day=day_m_nlpy;    % the #days of each month of a non leap year using the HSY definition b
           HS_hsy=zeros(nc,days);
           HS_hsy_filled=zeros(nc,days);
           e1_b=83;
           s1_b=161;  % Jun. 10
           e2_b=243;  % Aug.31
           s2_b=84; 
           e3_b=205;
           s3_b=244;  % Sept.1   
           e4_b=365;
           s4_b=206;
           e5_b=448;
           e6_b=243;
           e1_a=24;
           s1_a=342;  %Dec.8      
           e2_a=365;
           s2_a=25;
           e3_a=389;
           s3_a=390;
           e4_a=448;
           e5_a=59;
      else
         days=days_lpy;
%          m_day=day_m_lpy;    % condition (2): 1986, 1987, 1988
         HS_hsy=zeros(nc,days);
         HS_hsy_filled=zeros(nc,days);
         e1_b=83;
         s1_b=161;  % Jun. 10
         e2_b=243;  % Aug.31
         s2_b=84; 
         e3_b=205;
         s3_b=244;  % Sept.1   
         e4_b=365;
         s4_b=206;
         e5_b=449;
         e6_b=244;
         e1_a=24;
         s1_a=342;  %Dec.8      
         e2_a=365;
         s2_a=25;
         e3_a=389;
         s3_a=390;
         e4_a=449;
         e5_a=60;
      end
    else
         days=days_nlpy;
%          m_day=day_m_nlpy;    % condition (3): 1987, 1988, 1989
         HS_hsy=zeros(nc,days);
         HS_hsy_filled=zeros(nc,days);
         e1_b=83;
         s1_b=162;  % Jun. 10
         e2_b=244;  % Aug.31
         s2_b=84; 
         e3_b=205;
         s3_b=245;  % Sept.1   
         e4_b=366;
         s4_b=206;
         e5_b=448;
         e6_b=243;
         e1_a=23;
         s1_a=343;  %Dec.9      
         e2_a=365;
         s2_a=24;
         e3_a=389;
         s3_a=390;
         e4_a=448;
         e5_a=59;
    end
  else
       days=days_nlpy;
%        m_day=day_m_lpy;    % condition (4): 1988, 1989, 1990
       HS_hsy=zeros(nc,days);
       HS_hsy_filled=zeros(nc,days);
       e1_b=83;
       s1_b=161;  % Jun. 10
       e2_b=243;  % Aug.31
       s2_b=84; 
       e3_b=205;
       s3_b=244;  % Sept.1   
       e4_b=365;
       s4_b=206;
       e5_b=448;
       e6_b=243;
       e1_a=24;
       s1_a=343;  %Dec.8      
       e2_a=366;
       s2_a=25;
       e3_a=389;
       s3_a=390;
       e4_a=448;
       e5_a=59;                             
  end  
      
      coor_nnan=zeros(10,2);
      n_mcc=0;
  for n=1:nc      
    if (x_sel(n) == 4)  % HSY definition b: e.g. Jun.10,1986-Aug.31,1987
         HS_hsy(n,1:e1_b) = squeeze(HS_hsy_2(n,s1_b:e2_b));
         HS_hsy(n,s2_b:e3_b) = squeeze(HS_hsy_2(n,s3_b:e4_b));
         HS_hsy(n,s4_b:e5_b) = squeeze(HS_hsy_3(n,1:e6_b));
    elseif (x_sel(n) == 10) % HSY definition a: e.g. Dec.8/9,1985 - Feb.28/29,1987
             HS_hsy(n,1:e1_a) = squeeze(HS_hsy_1(n,s1_a:e2_a));
             HS_hsy(n,s2_a:e3_a) = squeeze(HS_hsy_2(n,:));
             HS_hsy(n,s3_a:e4_a) = squeeze(HS_hsy_3(n,1:e5_a));
    end
    
    % fill random gaps in a time series of HS
    if (sum(isnan(HS_hsy(n,:)))>0 && sum(isnan(HS_hsy(n,:)))~=days)
          ts=(1:days);
          HS_hsy_filled(n,:)=fillmissing(HS_hsy(n,:),'linear','SamplePoints',ts);
    else
          HS_hsy_filled(n,:)=HS_hsy(n,:);
    end
    
    if (sum(isnan(HS_hsy_filled(n,:)))~=days)           
          n_mcc=n_mcc+1;
          coor_nnan(n_mcc,:)=cc(n,:);
          DHDmax=0;          
      for i=e4_a:-1:84
            DHD=sum(HS_hsy_filled(n,(i-83):i));
            DHDmax=max(DHDmax,DHD);
      end
          dhd(n_mcc,(y-1))=DHDmax;       
    end
  end
end

toc




%%
% write out DHD
ncnc= netcdf.create('DHD_MMMct5km_cc.nc','NC_WRITE');   % Write netCDF file  

nID=netcdf.defDim(ncnc,'the number of coral cells',n_mcc);
coorID=netcdf.defDim(ncnc,'two columns for coordinate',2);

tID=netcdf.defDim(ncnc,'time',34);
vtID=netcdf.defVar(ncnc,'time','float',tID);
netcdf.putAtt(ncnc,vtID,'axis','T');
long_namet = 'the heat stress year corresponding to each DHD';
netcdf.putAtt(ncnc,vtID,'long_name',long_namet);

varname = 'DHD_hsy';
long_name = 'the DHD within each HSY from HS_ct5km_cc';
unit = 'degree celcius*day';
vmmmID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vmmmID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vmmmID,'units',unit);          % The unit
      
var = 'coor_cc';
long_name = 'coordinate of coral cells without NaNs';
unit = 'degree celcius';
vcoorID=netcdf.defVar(ncnc,var,'float',[nID,coorID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vcoorID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vcoorID,'units',unit);          % The unit         

% end define mode
netcdf.endDef(ncnc)
% input data
netcdf.putVar(ncnc,vmmmID,dhd);
netcdf.putVar(ncnc,vcoorID,coor_nnan);
netcdf.putVar(ncnc,vtID,HSY);
netcdf.close(ncnc)      




